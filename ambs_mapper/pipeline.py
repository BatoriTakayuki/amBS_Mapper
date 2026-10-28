"""Analysis pipeline: primers -> amplicons -> SNPs -> per-read calls -> reports."""
from __future__ import annotations

from collections import Counter
import csv
import gzip
import json
import multiprocessing
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import zipfile

from . import __version__, bsmap, igv_view, indexes, report
from .caller import DEFAULT_PARAMS, Caller, Genotype, Sample, Target
from .progress import emit, log
from .reference import find_amplicon, parse_region, reference_fasta, region_amplicon, revcomp
from .variants import query_region, resolve_strain, strain_list, strand_indels, strand_snps

READ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".bam")
TAG = re.compile(r"[ACGT]*")
NAME = re.compile(r"[^\\/:*?\"<>|\x00-\x1f]{1,60}")
CHUNK = 500


def is_read_file(path):
    # "._name" files are macOS metadata copied alongside the real files.
    return path.is_file() and not path.name.startswith(".") and path.name.lower().endswith(READ_SUFFIXES)


def list_inputs(input_root):
    root = Path(input_root)
    files = []
    if root.is_dir():
        for path in sorted(root.rglob("*")):
            if is_read_file(path):
                files.append(dict(path=path.relative_to(root).as_posix(), size=path.stat().st_size))
    return files


def expand_inputs(input_root, selected):
    root = Path(input_root).resolve()
    files = []
    for entry in selected:
        path = (root / entry).resolve()
        if not path.is_relative_to(root) or not path.exists():
            raise ValueError(f"入力ファイルが見つかりません: {entry}")
        candidates = sorted(path.rglob("*")) if path.is_dir() else [path]
        for candidate in candidates:
            if is_read_file(candidate):
                relative = candidate.relative_to(root).as_posix()
                if relative not in [item[0] for item in files]:
                    files.append((relative, candidate))
    if not files:
        raise ValueError("解析するFASTQが選択されていません。")
    return files


def validate_request(payload):
    """Normalise a submitted analysis request, raising ValueError with a user-facing message."""
    if not isinstance(payload, dict):
        raise ValueError("入力が不正です。")
    params = dict(DEFAULT_PARAMS)
    for key, value in (payload.get("params") or {}).items():
        if key not in DEFAULT_PARAMS:
            raise ValueError(f"不明な設定です: {key}")
        if key in ("exclude_primer_cpg", "single_index"):
            params[key] = bool(value)
        elif key == "max_unconverted" and value in (None, ""):
            params[key] = None
        else:
            try:
                number = float(value)
            except (TypeError, ValueError):
                raise ValueError(f"設定{key}は数値で指定してください。") from None
            if number < 0:
                raise ValueError(f"設定{key}は0以上にしてください。")
            params[key] = int(number) if key in ("max_unconverted", "tag_max_errors", "tag_slack", "snp_min_quality",
                                                  "photo_rows", "min_length", "auto_min_reads") else number
    layout = payload.get("layout") or "inline"
    if layout not in indexes.LAYOUTS:
        raise ValueError(f"不明なタグ形式です: {layout}")
    swapped = []
    amplicons = []
    for row in payload.get("amplicons") or []:
        name = str(row.get("name", "")).strip()
        if not NAME.fullmatch(name):
            raise ValueError("アンプリコン名を入力してください（60文字以内、/ \\ : * ? \" < > | 以外）。")
        if name in [item["name"] for item in amplicons]:
            raise ValueError(f"アンプリコン名が重複しています: {name}")
        max_size = int(float(row.get("max_size") or 5000))
        if not 50 <= max_size <= 50000:
            raise ValueError("最大産物長は50～50,000 bpで指定してください。")
        region = str(row.get("region") or "").strip()
        strand = str(row.get("strand") or "+").strip()
        forward, reverse = str(row.get("forward") or "").strip(), str(row.get("reverse") or "").strip()
        if region:
            parse_region(region)
            if strand not in ("+", "-"):
                raise ValueError(f"アンプリコン{name}: 鎖は + か - で指定してください。")
        elif not forward or not reverse:
            raise ValueError(f"アンプリコン{name}: プライマー（ForwardとReverse）か領域を入力してください。")
        amplicons.append(dict(name=name, forward=forward, reverse=reverse, max_size=max_size, region=region,
                              strand=strand))
    if not amplicons:
        raise ValueError("アンプリコン（プライマーか領域）を1つ以上入力してください。")
    samples = []
    for row in payload.get("samples") or []:
        name = str(row.get("name", "")).strip()
        if not NAME.fullmatch(name):
            raise ValueError("サンプル名を入力してください（60文字以内、/ \\ : * ? \" < > | 以外）。")
        amplicon = str(row.get("amplicon", "")).strip() or (amplicons[0]["name"] if len(amplicons) == 1 else "")
        if amplicon not in [item["name"] for item in amplicons]:
            raise ValueError(f"サンプル{name}のアンプリコン「{amplicon}」が定義されていません。")
        if layout == "nextera":
            # i7 entered in the i5 column (or the reverse) is recognised from the kit tables.
            kinds = (indexes.kind_of(row.get("tag_f")), indexes.kind_of(row.get("tag_r")))
            if "i7" in kinds[:1] or "i5" in kinds[1:]:
                row = {**row, "tag_f": row.get("tag_r"), "tag_r": row.get("tag_f")}
                swapped.append(name)
        tags = []
        for key, kind in (("tag_f", "i5"), ("tag_r", "i7")):
            # Kept as entered (Nextera names stay names) so a saved request restores the form.
            tag = re.sub(r"\s+", "", str(row.get(key) or "")).upper()
            if layout == "nextera":
                indexes.parse_index(tag, kind)
            elif not TAG.fullmatch(tag):
                raise ValueError(f"サンプル{name}のタグにACGT以外の文字があります。")
            tags.append(tag)
        mother, father = str(row.get("mother") or "").strip(), str(row.get("father") or "").strip()
        if not mother or not father:
            raise ValueError(f"サンプル{name}の母系統・父系統を入力してください。")
        samples.append(dict(name=name, group=str(row.get("group") or "").strip(), amplicon=amplicon,
                            tag_f=tags[0], tag_r=tags[1], mother=mother, father=father,
                            fastq=str(row.get("fastq") or "").strip().strip("/")))
    # Without a sample table, samples are made from the tag combinations found
    # (Nextera) or from all reads of each amplicon (inline).
    auto = None
    if not samples:
        mother = str(payload.get("auto_mother") or "").strip()
        father = str(payload.get("auto_father") or "").strip()
        if not mother or not father:
            raise ValueError("サンプルを入力するか、サンプル表を空にする場合は「自動サンプル」の母系統・父系統を入力してください。")
        auto = dict(mother=mother, father=father)
    keys = [(s["amplicon"], s["tag_f"], s["tag_r"], s["fastq"]) for s in samples]
    if len(set(keys)) != len(keys):
        raise ValueError("同じアンプリコン・タグ・FASTQの組み合わせのサンプルが重複しています。")
    if len({(s["name"], s["amplicon"]) for s in samples}) != len(samples):
        raise ValueError("同じアンプリコンで同じサンプル名が重複しています。")
    mode = payload.get("mode", "analyze")
    inputs = payload.get("inputs") or []
    if isinstance(inputs, str):
        inputs = [inputs]
    if not all(isinstance(item, str) for item in inputs) or (mode == "analyze" and not inputs):
        raise ValueError("解析するFASTQを選択してください。")
    threads = int(payload.get("threads") or os.cpu_count() or 1)
    title = str(payload.get("title") or "").strip()[:80]
    notes = []
    entered = [s[key] for s in samples for key in ("tag_f", "tag_r") if s[key]]
    if layout == "inline" and entered and sum(1 for tag in entered if indexes.kind_of(tag)) >= 0.8 * len(entered):
        notes.append("タグのほとんどがNextera XTのインデックス配列です。Nexteraライブラリなら、タグ形式を"
                     "「Nextera デュアルインデックス」にしてください（インラインではアダプター分ずれてタグが見つかりません）。")
    if swapped:
        notes.append(f"i5とi7が逆の列に入っていたため入れ替えました（{', '.join(swapped[:5])}"
                     f"{' ほか' if len(swapped) > 5 else ''}、{len(swapped)}サンプル）。")
    return dict(title=title, inputs=inputs, amplicons=amplicons, samples=samples, params=params,
                threads=max(1, min(threads, 256)), mode=mode, layout=layout, auto=auto, notes=notes)


def validate_explore(payload):
    """Request for the read-structure survey: only files are needed."""
    inputs = payload.get("inputs") or []
    if isinstance(inputs, str):
        inputs = [inputs]
    if not inputs or not all(isinstance(item, str) for item in inputs):
        raise ValueError("調べるFASTQを選択してください。")
    try:
        max_reads = int(payload.get("max_reads") or 300000)
    except (TypeError, ValueError):
        raise ValueError("最大リード数は数値で指定してください。") from None
    threads = int(payload.get("threads") or os.cpu_count() or 1)
    title = "リード構造の調査"
    return dict(mode="explore", title=title, inputs=inputs, max_reads=max(1000, max_reads),
                threads=max(1, min(threads, 256)))


def tag_warnings(tags, max_errors):
    """Pairs of tags that cannot be told apart with the allowed errors, as at most one note."""
    import edlib
    close = []
    tags = sorted(set(tag for tag in tags if tag))
    for i, first in enumerate(tags):
        for second in tags[i + 1:]:
            distance = edlib.align(first, second, mode="NW")["editDistance"]
            if distance <= 2 * max_errors:
                close.append((distance, first, second))
    if not close:
        return []
    close.sort()
    shown = "、".join(f"{a}/{b}（{d}）" for d, a, b in close[:3])
    return [f"距離が近いタグの組が{len(close)}組あります（{shown}{' など' if len(close) > 3 else ''}）。"
            f"許容誤差{max_errors}では誤判定の可能性があるので、許容誤差を下げることを検討してください。"]


def strain_sequence(sequence, product, variants, strain):
    bases = list(sequence)
    for variant in variants:
        allele = variant["alleles"].get(strain)
        if allele is None or len(allele) != 1 or len(variant["ref"]) != 1:
            continue
        if not product["start"] <= variant["pos"] < product["end"]:
            continue
        if product["strand"] == "+":
            bases[variant["pos"] - product["start"]] = allele
        else:
            bases[product["end"] - 1 - variant["pos"]] = revcomp(allele)
    return "".join(bases)


def build_targets(request, resources):
    available = strain_list(resources)
    layout = request.get("layout", "inline")
    samples = []
    for index, row in enumerate(request["samples"]):
        tag_f, tag_r = row["tag_f"], row["tag_r"]
        if layout == "nextera":
            tag_f, tag_r = indexes.parse_index(tag_f, "i5"), indexes.parse_index(tag_r, "i7")
        samples.append(Sample(index=index, name=row["name"], group=row["group"], amplicon=row["amplicon"],
                              tag_f=tag_f, tag_r=tag_r, fastq=row["fastq"],
                              mother=resolve_strain(row["mother"], available),
                              father=resolve_strain(row["father"], available),
                              labels=tag_labels(layout, tag_f, tag_r)))
    auto = request.get("auto")
    auto_pair = (resolve_strain(auto["mother"], available), resolve_strain(auto["father"], available)) if auto else None
    targets, details, warnings = [], [], list(request.get("notes") or [])
    for row in request["amplicons"]:
        if row.get("region"):
            product = region_amplicon(resources, row["name"], row["region"], row["strand"], row["forward"], row["reverse"])
        else:
            emit("amplicon", name=row["name"], message=f"{row['name']}: プライマー位置をゲノム上で検索しています。")
            product = find_amplicon(resources, row["name"], row["forward"], row["reverse"], row["max_size"])
        if product.get("swapped"):
            warnings.append(f"{row['name']}: Forward/Reverseを入れ替えると一致したため、入れ替えて解析しました。")
        members = [sample for sample in samples if sample.amplicon == row["name"]]
        pairs = sorted({sample.pair for sample in members} | ({auto_pair} if auto_pair else set()))
        strains = sorted({strain for pair in pairs for strain in pair})
        emit("variants", name=row["name"], message=f"{row['name']}: MGPから系統SNPを取得しています。")
        variants = query_region(resources, product["chrom"], product["start"], product["end"], strains)
        indels = {}
        for mother, father in pairs:
            if mother != father:
                try:
                    indels[(mother, father)] = strand_indels(resources, product, mother, father)
                except RuntimeError as error:
                    warnings.append(f"{row['name']}: indel情報を取得できませんでした（{error}）。")
        target, detail, notes = make_target(row["name"], product, variants, members, request["params"], indels,
                                            layout=layout, auto_pair=auto_pair)
        targets.append(target)
        details.append(detail)
        warnings += [f"{row['name']}: {text}" for text in notes]
    return samples, targets, details, list(dict.fromkeys(warnings))


def tag_labels(layout, tag_f, tag_r):
    if layout == "nextera":
        return indexes.label(tag_f, "i5"), indexes.label(tag_r, "i7")
    return tag_f, tag_r


def make_target(name, product, variants, members, params, indels=None, layout="inline", auto_pair=None):
    """Analysis target for one amplicon from its product, MGP variants and samples."""
    indels = indels or {}
    notes = []
    pairs_used = sorted({sample.pair for sample in members} | ({auto_pair} if auto_pair else set()))
    strains = sorted({strain for pair in pairs_used for strain in pair})
    sequence, length = product["sequence"], len(product["sequence"])
    # Region-defined amplicons without primers have no primer CpGs to exclude.
    lf, lr = product.get("primer_lengths") or (len(product["forward"]), len(product["reverse"]))
    strain_seqs = {strain: strain_sequence(sequence, product, variants, strain) for strain in strains}
    in_primer = (lambda x: x < lf or x + 1 >= length - lr) if params["exclude_primer_cpg"] else (lambda x: False)
    sites = sorted({x for seq in strain_seqs.values() for x in range(length - 1)
                    if seq[x] == "C" and seq[x + 1] == "G" and not in_primer(x)})
    variable = {}
    for variant in variants:
        if len(variant["ref"]) != 1 or not product["start"] <= variant["pos"] < product["end"]:
            continue
        x = variant["pos"] - product["start"] if product["strand"] == "+" else product["end"] - 1 - variant["pos"]
        variable[x] = {seq[x] for seq in strain_seqs.values()} | {sequence[x]}
    conversion = [x for x in range(length - 1) if sequence[x] == "C" and x not in variable
                  and lf <= x < length - lr and all(seq[x + 1] != "G" for seq in strain_seqs.values())]
    target = Target(name=name, sequence=sequence, forward=product["forward"], reverse=product["reverse"],
                    chrom=product["chrom"], start=product["start"], end=product["end"], strand=product["strand"],
                    sites=sites, conversion=conversion, variable=set(variable), variants=variable,
                    samples=members, tags_f=sorted({s.tag_f for s in members if s.tag_f}),
                    tags_r=sorted({s.tag_r for s in members if s.tag_r}), layout=layout, auto_pair=auto_pair)
    if layout != "nextera":
        # Nextera kit indexes are designed to be distinct, and every kit index is a
        # candidate (a unique best match is required), so only inline tags are checked.
        notes += tag_warnings(target.tags_f, params["tag_max_errors"])
        notes += tag_warnings(target.tags_r, params["tag_max_errors"])
    if layout == "nextera":
        # Every kit index is a candidate, so a read from an index left out of the
        # sample table is reported as such instead of being forced onto a neighbour.
        target.tags_f = sorted(set(target.tags_f) | set(indexes.whitelist("i5")))
        target.tags_r = sorted(set(target.tags_r) | set(indexes.whitelist("i7")))
    pairs = []
    to3 = str.maketrans("C", "T")
    for pair in pairs_used:
        mother, father = pair
        snps, uncertain = strand_snps(product, variants, mother, father)
        informative = [(s["x"], s["mother"], s["father"]) for s in snps if s["informative"]]
        seq_m, seq_f = strain_seqs[mother], strain_seqs[father]
        target.genotypes[pair] = Genotype(
            snps=informative,
            mother_sites=[seq_m[x] == "C" and seq_m[x + 1] == "G" for x in sites],
            father_sites=[seq_f[x] == "C" and seq_f[x + 1] == "G" for x in sites],
            mother3=seq_m.translate(to3), father3=seq_f.translate(to3))
        pair_indels = indels.get(pair, [])
        # Reciprocal crosses share one note.
        if mother != father and (father, mother) not in target.genotypes:
            if not informative:
                notes.append(f"{mother}と{father}を区別できるSNPがアンプリコン内にありません"
                             "（C/T型のSNPはバイサルファイト処理後に区別できません）。アレルは判定不能になります。")
            if pair_indels:
                notes.append(f"{mother}/{father}間のindelが{len(pair_indels)}か所あります（SNPのみでアレルを判定しています）。")
        pairs.append(dict(mother=mother, father=father, snps=snps, informative=len(informative),
                          uncertain=uncertain, indels=pair_indels))
    detail = dict(name=name, chrom=product["chrom"], start=product["start"] + 1, end=product["end"],
                  strand=product["strand"], size=length, forward=product["forward"],
                  reverse=product["reverse"], sequence=sequence, sites=sites,
                  # 1-based top-strand position of the CpG's C
                  site_positions=[(product["start"] + x + 1) if product["strand"] == "+"
                                  else (product["end"] - x - 1) for x in sites],
                  pairs=pairs, conversion_sites=len(conversion))
    return target, detail, notes


def read_records(path):
    """(name, sequence, quality) from FASTQ(.gz) or unaligned BAM."""
    if path.name.lower().endswith(".bam"):
        process = subprocess.Popen(["samtools", "fastq", "-0", "/dev/stdout", str(path)],
                                   stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        handle = process.stdout
    else:
        handle = (gzip.open(path, "rt", encoding="ascii", errors="replace") if path.name.lower().endswith(".gz")
                  else open(path, encoding="ascii", errors="replace"))
    with handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence, _, quality = handle.readline().strip(), handle.readline(), handle.readline().strip()
            if not header.startswith("@"):
                raise ValueError(f"FASTQ形式ではありません: {path.name}")
            yield header[1:].split()[0], sequence, quality


def _chunks(files):
    chunk = []
    for relative, path in files:
        for name, sequence, quality in read_records(path):
            chunk.append((name, sequence, quality, relative))
            if len(chunk) >= CHUNK:
                yield chunk
                chunk = []
    if chunk:
        yield chunk


_caller = None
_mapper = None                  # loaded once in the parent; forked workers share the genome index


def _init_worker(targets, params, index_path):
    global _caller
    _caller = Caller(targets, params, index_path, _mapper)


def _work(chunk):
    return [_caller.call(*item) for item in chunk]


FIELDS = ["read", "source", "length", "amplicon", "orientation", "sample", "tag_f", "tag_r", "tag_f_label",
          "tag_r_label", "tag_sides", "status", "coverage", "identity", "mismatch", "conversion", "unconverted", "votes_m", "votes_f", "allele", "pattern",
          "chrom", "gstart", "gend", "bs_strand", "mapq"]


HAPLOTYPE = {"mother": 1, "father": 2}


def sam_line(record, sample, read_groups, layout="inline"):
    """The worker's SAM line plus the final sample, allele, status and tags found as SAM tags."""
    tags = [f"YT:Z:{record['status']}"]
    if record["amplicon"]:
        tags.append(f"YA:Z:{record['amplicon']}")
    # Index (or inline tag) found on each side, or why none was: i5/i7 for Nextera, tf/tr inline.
    for key, value in zip(("i5", "i7") if layout == "nextera" else ("tf", "tr"),
                          (record["tag_f_label"], record["tag_r_label"])):
        if value:
            tags.append(f"{key}:Z:{value}")
    if sample:
        if sample not in read_groups:
            # read-group IDs stay ASCII; the sample name itself goes to SM
            base = report.safe(sample)
            used = set(read_groups.values())
            read_groups[sample] = base if base not in used else f"{base}_{len(read_groups) + 1}"
        tags.append(f"RG:Z:{read_groups[sample]}")
    if record["allele"] in HAPLOTYPE:
        tags.append(f"HP:i:{HAPLOTYPE[record['allele']]}")   # IGV: Group alignments by phase
    return record["sam"] + "\t" + "\t".join(tags)


def write_bam(path, body_path, read_groups, resources, threads):
    """Coordinate-sorted, indexed BAM of every mapped read, for IGV."""
    fai = Path(str(reference_fasta(resources)) + ".fai")
    header = ["@HD\tVN:1.6"]
    header += [f"@SQ\tSN:{name}\tLN:{length}"
               for name, length, *_ in (line.split("\t") for line in fai.read_text().splitlines() if line)]
    header += [f"@RG\tID:{rid}\tSM:{name}" for name, rid in read_groups.items()]
    header.append(f"@PG\tID:ambs-mapper\tPN:ambs-mapper\tVN:{__version__}")
    command = ["samtools", "sort", "-@", str(max(1, min(threads, 8))), "-m", "768M",
               "-T", str(body_path.with_suffix(".tmp")), "-o", str(path), "-"]
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        process.stdin.write(("\n".join(header) + "\n").encode("utf-8"))
        with open(body_path, "rb") as body:
            shutil.copyfileobj(body, process.stdin, length=1024 * 1024)
        process.stdin.close()
    except BrokenPipeError:
        pass
    error = process.stderr.read().decode("utf-8", "replace")
    if process.wait() != 0:
        raise RuntimeError(f"BAMを作成できませんでした: {error.strip()[-300:]}")
    subprocess.run(["samtools", "index", str(path)], check=True, capture_output=True)
    body_path.unlink()


def file_name(text):
    """A file name that keeps Japanese but drops characters Windows does not allow."""
    return re.sub(r'[\\/:*?"<>|\s]+', "_", text).strip("._") or "x"


def split_bam(path, folder, read_groups, samples, amplicons=()):
    """One BAM per sample holding its passed reads, so IGV shows each sample as a track,
    plus one "タグ不明" BAM per amplicon.

    Files are named after the sample (IGV uses the file name as the track name)
    and zipped for a single download. Returns the BAMs written as dicts
    (path, name, amplicon, unknown), samples first.
    """
    folder.mkdir(exist_ok=True)
    wanted = {read_groups[s.name]: s for s in samples if s.name in read_groups}
    passed = folder / "passed.tmp.bam"
    subprocess.run(["samtools", "view", "-b", "-d", "YT:pass", "-o", str(passed), str(path)],
                   check=True, capture_output=True)
    subprocess.run(["samtools", "split", "-u", str(folder / "other.tmp.bam"), "-f", str(folder / "%!.tmp.%."),
                    str(passed)], check=True, capture_output=True)
    written = []
    for part in folder.glob("*.tmp.bam"):
        rid = part.name[:-len(".tmp.bam")]
        if rid in wanted:
            sample = wanted[rid]
            target = folder / f"{file_name(sample.name)}_{file_name(sample.amplicon)}.bam"
            part.replace(target)
            subprocess.run(["samtools", "index", str(target)], check=True, capture_output=True)
            written.append(dict(path=target, name=sample.name, amplicon=sample.amplicon, unknown=False))
        else:
            part.unlink()            # reads outside the reported samples (or none at all)
    written.sort(key=lambda entry: (entry["amplicon"], entry["name"]))
    # Reads on an amplicon whose tag could not be read or fits several samples, for
    # checking what the unassigned reads look like; YT tells the two apart.
    for amplicon in amplicons:
        target = folder / f"タグ不明_{file_name(amplicon)}.bam"
        expression = f'[YA]=="{amplicon}" && ([YT]=="tag_unknown" || [YT]=="tag_ambiguous")'
        subprocess.run(["samtools", "view", "-b", "-e", expression, "-o", str(target), str(path)],
                       check=True, capture_output=True)
        count = subprocess.run(["samtools", "view", "-c", str(target)], check=True, capture_output=True, text=True)
        if int(count.stdout.strip() or 0):
            subprocess.run(["samtools", "index", str(target)], check=True, capture_output=True)
            written.append(dict(path=target, name=f"タグ不明（{amplicon}）", amplicon=amplicon, unknown=True))
        else:
            target.unlink()
    with zipfile.ZipFile(folder.with_suffix(".zip"), "w", zipfile.ZIP_STORED) as archive:
        for entry in written:
            target = entry["path"]
            archive.write(target, f"{folder.name}/{target.name}")
            archive.write(str(target) + ".bai", f"{folder.name}/{target.name}.bai")
    return written


def run(job_dir, input_root, resources):
    job_dir = Path(job_dir)
    request = json.loads((job_dir / "request.json").read_text(encoding="utf-8"))
    out = job_dir / "output"
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()
    if request.get("mode") != "check":
        files = expand_inputs(input_root, request["inputs"])
    if request.get("mode") == "explore":
        from . import explore
        explore.run(out, request, files, resources, read_records)
        emit("done", message="リード構造の調査が完了しました。")
        return
    samples, targets, details, warnings = build_targets(request, resources)
    for text in warnings:
        log("注意: " + text)
    (out / "amplicons.json").write_text(json.dumps(dict(amplicons=details, warnings=warnings), ensure_ascii=False,
                                                   indent=1), encoding="utf-8")
    index_path = job_dir / "work" / "amplicons_ct.fa"
    index_path.parent.mkdir(exist_ok=True)
    index_path.write_text("".join(f">{t.name}\n{t.sequence.translate(str.maketrans('C', 'T'))}\n" for t in targets),
                          encoding="ascii")
    if request.get("mode") == "check":
        report.write_check(out, request, samples, details, warnings)
        emit("done", message="プライマーとSNPの確認が完了しました。")
        return

    global _mapper
    if not bsmap.index_ready(resources):
        raise RuntimeError("ゲノム索引がありません。画面の「参照ゲノムを準備」を実行してください。")
    emit("reads", message="ゲノム索引を読み込んでいます。")
    _mapper = bsmap.BisulfiteMapper(bsmap.index_path(resources))
    total_bytes = sum(path.stat().st_size for _, path in files)
    emit("reads", message=f"{len(files)}ファイルを解析しています。", done=0, files=len(files))
    collector = report.Collector(samples, targets)
    auto_samples = {}
    by_name = {target.name: target for target in targets}

    def auto_sample(amplicon, key):
        if (amplicon, key) not in auto_samples:
            target = by_name[amplicon]
            labels = tag_labels(target.layout, *key)
            name = "-".join(label for label in labels if label) or "全リード"
            # Indexes stay unique after rare combinations are dropped from the list.
            sample = Sample(index=len(request["samples"]) + len(auto_samples), name=name, group="", amplicon=amplicon,
                            tag_f=key[0], tag_r=key[1],
                            mother=target.auto_pair[0], father=target.auto_pair[1], labels=labels)
            samples.append(sample)
            auto_samples[(amplicon, key)] = sample.index
        return auto_samples[(amplicon, key)]

    processed = 0
    last = 0.0
    pending = []                # automatic samples: single-index reads, resolved at the end
    body_path = job_dir / "work" / "reads.body.sam"
    read_groups = {}            # sample name -> read-group ID in the BAM
    with gzip.open(out / "reads.tsv.gz", "wt", encoding="utf-8", newline="") as handle, \
            open(body_path, "w", encoding="utf-8") as body:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        def keep(record):
            collector.add(record)
            row = dict(record)
            row["sample"] = by_index[record["sample"]].name if record["sample"] >= 0 else ""
            writer.writerow(row)
            if record.get("sam"):
                body.write(sam_line(record, row["sample"], read_groups, request.get("layout", "inline")) + "\n")

        by_index = {}
        context = multiprocessing.get_context("fork")
        with context.Pool(request["threads"], _init_worker, (targets, request["params"], index_path)) as pool:
            for records in pool.imap(_work, _chunks(files), chunksize=1):
                for record in records:
                    key = record.pop("auto_key", None)
                    if key is not None and by_name[record["amplicon"]].layout == "nextera" and not all(key):
                        pending.append((key, record))
                        continue
                    if key is not None:
                        record["sample"] = auto_sample(record["amplicon"], key)
                    by_index.update((s.index, s) for s in samples[len(by_index):])
                    keep(record)
                processed += len(records)
                if time.time() - last > 2:
                    last = time.time()
                    emit("reads", done=processed, passed=collector.passed, total_bytes=total_bytes,
                         message=f"{processed:,}リードを処理（合格 {collector.passed:,}）")
        if request.get("auto"):
            # Rare combinations are index errors, chimeras or carry-over, not samples.
            minimum = request["params"]["auto_min_reads"]
            created = len(auto_samples)
            kept = [s for s in samples if collector.sample_status[s.index]["pass"] >= minimum]
            for sample in samples:
                if sample not in kept:
                    collector.drop(sample)
            samples[:] = kept
            full = {s.amplicon for s in kept}
            # Amplicons without any dual-indexed sample can only be split by the one
            # index each read carries; those groups are reported as such.
            single_pass = Counter((record["amplicon"], key) for key, record in pending
                                  if record["amplicon"] not in full and record["status"] == "pass")
            single_samples = 0
            for (amplicon, key), count in single_pass.items():
                if count >= minimum:
                    index = auto_sample(amplicon, key)
                    sample = next(s for s in samples if s.index == index)
                    sample.name += "（片側のみ）"
                    by_index[index] = sample
                    single_samples += 1
            rescued = 0
            for key, record in pending:
                matches = [s for s in samples if s.amplicon == record["amplicon"]
                           and (not key[0] or s.tag_f == key[0]) and (not key[1] or s.tag_r == key[1])]
                if record["amplicon"] not in full:
                    matches = [s for s in matches if (s.tag_f, s.tag_r) == key]
                if len(matches) == 1:
                    record["sample"] = matches[0].index
                    rescued += 1
                else:
                    record["status"] = "tag_ambiguous" if matches else "tag_unknown"
                keep(record)
            samples.sort(key=lambda s: (s.amplicon, s.name))
            warnings.append(f"サンプル表が空のため、タグの組み合わせごとに自動でサンプルを作成しました"
                            f"（合格リード{minimum}未満の{created - len(kept)}組み合わせは除外）。")
            if pending:
                warnings.append(f"片側のインデックスしか読めないリード{len(pending):,}本のうち、"
                                f"{rescued:,}本をサンプルに割り当てました。")
            if single_samples:
                names = ", ".join(sorted(amplicon for amplicon in {a for a, _ in single_pass} if amplicon not in full))
                warnings.append(f"{names}: 両側のインデックスがそろったリードがないため、片側のインデックスごとにまとめました"
                                "（「片側のみ」のサンプルは複数サンプルの混合の可能性があります）。")
    emit("reads", done=processed, passed=collector.passed, message=f"{processed:,}リードを処理しました。")
    emit("bam", message="IGV用のBAMを作成しています。")
    write_bam(out / "reads.bam", body_path, read_groups, resources, request["threads"])
    split = split_bam(out / "reads.bam", out / "bam_by_sample", read_groups, samples,
                      [target.name for target in targets])
    sample_bams = sum(not entry["unknown"] for entry in split)
    emit("bam", message="画面のIGV用データを作成しています。")
    igv_tracks = igv_view.build(out, details, split)
    emit("report", message="図と表を作成しています。")
    report.write_all(out, request, samples, targets, details, warnings, collector,
                     dict(files=[relative for relative, _ in files], reads=processed,
                          seconds=round(time.time() - started, 1), bam=True, sample_bams=sample_bams,
                          igv=igv_tracks))
    emit("done", message="解析が完了しました。")
