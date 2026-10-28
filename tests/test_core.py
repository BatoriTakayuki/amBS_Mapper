"""Offline tests on a synthetic genome: in-silico PCR, SNP logic, per-read calls."""
import random

import pytest

from ambs_mapper import indexes
from ambs_mapper.caller import DEFAULT_PARAMS, Caller, Sample, anchored_tag, best_tag
from ambs_mapper.pipeline import make_target, strain_sequence, validate_request
from ambs_mapper.reference import Genome, in_silico_pcr, revcomp
from ambs_mapper.report import quma_order
from ambs_mapper.simulate import simulate
from ambs_mapper.variants import resolve_strain, strand_snps

STRAINS = ["C57BL_6J", "JF1_MsJ", "129S1_SvImJ", "C57BL_6NJ", "CAST_EiJ"]


def bisulfite_primer(sequence):
    """Primer as designed on the converted strand: C->T, CpG C -> Y."""
    out = []
    for i, base in enumerate(sequence):
        if base == "C":
            out.append("Y" if i + 1 < len(sequence) and sequence[i + 1] == "G" else "T")
        else:
            out.append(base)
    return "".join(out)


@pytest.fixture(scope="module")
def genome(tmp_path_factory):
    rng = random.Random(7)
    directory = tmp_path_factory.mktemp("genome")
    chroms = {"1": "".join(rng.choice("ACGT") for _ in range(30000)),
              "2": "".join(rng.choice("ACGT") for _ in range(20000))}
    fasta = directory / "reference.fa"
    lines, fai, offset = [], [], 0
    for name, sequence in chroms.items():
        header = f">{name}\n"
        offset += len(header)
        body = "".join(sequence[i:i + 60] + "\n" for i in range(0, len(sequence), 60))
        fai.append(f"{name}\t{len(sequence)}\t{offset}\t60\t61\n")
        lines.append(header + body)
        offset += len(body)
    fasta.write_text("".join(lines))
    (directory / "reference.fa.fai").write_text("".join(fai))
    return Genome(fasta), chroms


def test_fetch(genome):
    g, chroms = genome
    assert g.fetch("1", 59, 125) == chroms["1"][59:125]
    assert g.chromosome("2") == chroms["2"].encode()


def test_in_silico_pcr_top_strand(genome):
    g, chroms = genome
    top = chroms["1"][10000:10400]
    forward = bisulfite_primer(top[:24])
    reverse = revcomp(bisulfite_primer(top[-24:]))
    products = in_silico_pcr(g, forward, reverse, 1000)
    assert [(p["chrom"], p["start"], p["end"], p["strand"]) for p in products] == [("1", 10000, 10400, "+")]


def test_in_silico_pcr_bottom_strand(genome):
    g, chroms = genome
    x = revcomp(chroms["2"][5000:5600])        # bisulfite strand is the bottom strand
    forward = bisulfite_primer(x[:25])
    reverse = revcomp(bisulfite_primer(x[-25:]))
    products = in_silico_pcr(g, forward, reverse, 1000)
    assert [(p["chrom"], p["start"], p["end"], p["strand"]) for p in products] == [("2", 5000, 5600, "-")]
    assert in_silico_pcr(g, forward, reverse, 500) == []


def test_resolve_strain():
    assert resolve_strain("B6", STRAINS) == "C57BL_6J"
    assert resolve_strain("jf1", STRAINS) == "JF1_MsJ"
    assert resolve_strain("CAST", STRAINS) == "CAST_EiJ"
    with pytest.raises(ValueError):
        resolve_strain("C57BL", STRAINS)         # ambiguous prefix
    with pytest.raises(ValueError):
        resolve_strain("XYZ", STRAINS)


def test_strand_snps_bottom_and_uninformative():
    product = dict(chrom="1", start=100, end=200, strand="-")
    variants = [dict(pos=150, ref="G", alleles={"C57BL_6J": "G", "JF1_MsJ": "A"}),   # G/A top = C/T on X
                dict(pos=160, ref="A", alleles={"C57BL_6J": "A", "JF1_MsJ": "C"}),
                dict(pos=170, ref="A", alleles={"C57BL_6J": "A", "JF1_MsJ": None})]
    snps, uncertain = strand_snps(product, variants, "C57BL_6J", "JF1_MsJ")
    assert [(s["x"], s["mother"], s["father"], s["informative"]) for s in snps] == [
        (39, "T", "G", True), (49, "C", "T", False)]
    assert uncertain == 1


def test_best_tag_anchored():
    params = dict(DEFAULT_PARAMS)
    primer = "ATTTTTAAGAATTGGTGATTGTTG"
    read = "GGTTAACC" + "ACGTACGA" + primer + "TTTGGGAAATTT"
    tags = ["ACGTACGA", "TGCATGTC"]
    assert best_tag(read, primer, 16, tags, params)[0] == "ACGTACGA"
    # one substitution and one deletion still resolve
    assert best_tag("GGTTAACC" + "ACGAACG" + primer + "TTT", primer, 15, tags, params)[0] == "ACGTACGA"
    assert best_tag("GGTTAACCGGTTAACC" + primer, primer, 16, tags, params)[0] is None


def test_quma_order():
    assert quma_order(["MUM", "UUU", "MMU", "UMU"]) == ["UUU", "UMU", "MMU", "MUM"]


def test_validate_request():
    base = dict(inputs=["x.fastq"], amplicons=[dict(name="A", forward="ACGTACGTACGTAC", reverse="ACGTACGTACGTAC")],
                samples=[dict(name="s1", tag_f="acgt", mother="B6", father="JF1")])
    request = validate_request(base)
    assert request["samples"][0]["amplicon"] == "A" and request["samples"][0]["tag_f"] == "ACGT"
    with pytest.raises(ValueError):
        validate_request({**base, "samples": [dict(name="s1", tag_f="ACGN", mother="B6", father="JF1")]})
    with pytest.raises(ValueError):
        validate_request({**base, "samples": base["samples"] * 2})
    with pytest.raises(ValueError):
        validate_request({**base, "inputs": []})
    assert validate_request({**base, "inputs": [], "mode": "check"})["mode"] == "check"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_caller_on_simulated_reads(genome, tmp_path, strand):
    g, chroms = genome
    start, end = 20000, 20450
    top = chroms["1"][start:end]
    x = top if strand == "+" else revcomp(top)
    product = dict(chrom="1", start=start, end=end, strand=strand, sequence=x,
                   forward=bisulfite_primer(x[:24]), reverse=revcomp(bisulfite_primer(x[-24:])))
    # SNPs inside the amplicon (top-strand alleles); one creates a JF1-only CpG
    rng = random.Random(3)
    variants = []
    for pos in (start + 80, start + 150, start + 230, start + 310):
        ref = top[pos - start]
        alt = rng.choice([b for b in "ACGT" if b != ref and {b, ref} not in ({"C", "T"}, {"G", "A"})])
        variants.append(dict(pos=pos, ref=ref, alleles={"C57BL_6J": ref, "JF1_MsJ": alt}))
    samples = [Sample(0, "WT", "ctrl", "A", "ACGTACGA", "TTGCAGCA", "C57BL_6J", "JF1_MsJ"),
               Sample(1, "REC", "ctrl", "A", "TGCATGTC", "TTGCAGCA", "JF1_MsJ", "C57BL_6J"),
               Sample(2, "B6", "ctrl", "A", "GATCCTAG", "", "C57BL_6J", "C57BL_6J")]
    params = dict(DEFAULT_PARAMS)
    target, detail, _ = make_target("A", product, variants, samples, params)
    b6 = strain_sequence(x, product, variants, "C57BL_6J")
    jf1 = strain_sequence(x, product, variants, "JF1_MsJ")
    spec = [dict(name="WT", tag_f="ACGTACGA", tag_r="TTGCAGCA", methylation=dict(mother=0.05, father=0.95),
                 alleles=dict(mother=b6, father=jf1)),
            dict(name="REC", tag_f="TGCATGTC", tag_r="TTGCAGCA", methylation=dict(mother=0.05, father=0.95),
                 alleles=dict(mother=jf1, father=b6)),
            dict(name="B6", tag_f="GATCCTAG", tag_r="", methylation=dict(mother=0.5, father=0.5),
                 alleles=dict(mother=b6, father=b6))]
    simulate(tmp_path / "r.fastq", tmp_path / "t.jsonl", product, None, spec, 150, error_rate=0.04)
    index = tmp_path / "idx.fa"
    index.write_text(f">A\n{x.replace('C', 'T')}\n")
    caller = Caller([target], params, index)
    import json
    truth = {json.loads(line)["read"]: json.loads(line) for line in (tmp_path / "t.jsonl").read_text().splitlines()}
    lines = (tmp_path / "r.fastq").read_text().splitlines()
    names = {s.index: s.name for s in samples}
    ok_sample = wrong_sample = ok_allele = wrong_allele = calls = ok_calls = passed = 0
    for i in range(0, len(lines), 4):
        record = caller.call(lines[i][1:], lines[i + 1], lines[i + 3], "r.fastq")
        t = truth[record["read"]]
        if record["sample"] < 0:
            continue
        if names[record["sample"]] != t["sample"]:
            wrong_sample += 1
            continue
        ok_sample += 1
        if record["status"] != "pass":
            continue
        passed += 1
        if record["allele"] in ("mother", "father"):
            ok_allele += record["allele"] == t["allele"]
            wrong_allele += record["allele"] != t["allele"]
        methylated = set(t["methylated"])
        for site, value in zip(target.sites, record["pattern"]):
            if value in "MU":
                calls += 1
                ok_calls += (value == "M") == (site in methylated)
    assert ok_sample >= 0.9 * 450 and wrong_sample <= 2
    assert passed >= 0.85 * ok_sample
    assert ok_allele >= 0.8 * 300 and wrong_allele <= 1
    assert ok_calls / calls > 0.95
    # CpG positions are reported on the top strand at the C of the CpG
    for x_site, position in zip(detail["sites"], detail["site_positions"]):
        top_c = position - 1
        assert chroms["1"][top_c:top_c + 2] == "CG" or any(v["pos"] in (top_c, top_c + 1) for v in variants)


def test_parse_index():
    assert indexes.parse_index("S520", "i5") == "AAGGCTAT"
    assert indexes.parse_index("n712", "i7") == revcomp("GTAGAGGA")      # oligo carries the complement
    assert indexes.parse_index("712", "i7") == revcomp("GTAGAGGA")
    assert indexes.parse_index("acgtacgt", "i5") == "ACGTACGT"
    assert indexes.label(revcomp("GTAGAGGA"), "i7") == "N712"
    with pytest.raises(ValueError):
        indexes.parse_index("N799X", "i7")


def test_nextera_tag_through_adapter():
    params = dict(DEFAULT_PARAMS)
    primer = "ATTTTTAAGAATTGGTGATTGTTG"
    a1 = indexes.A1
    noisy_a1 = a1[:5] + a1[6:20] + "G" + a1[21:]                  # adapter with a deletion and a substitution
    i5 = indexes.I5["S520"]
    read = "AATGATACGGCGACCACCGAGATCTACAC" + i5[:3] + "G" + i5[4:] + noisy_a1 + primer + "TTTGGG"
    start = read.index(primer)
    assert anchored_tag(read, start, indexes.whitelist("i5"), params, a1)[0] == i5
    assert anchored_tag(read, start, indexes.whitelist("i5"), params, a1, indexes.P5_END)[0] == i5
    # Nanopore adapter directly on A1 (no P5/i5): its end resembles an index, the P5 check rejects it
    bare = "GTTCAGTTACGTATTGCT" + a1 + primer
    assert anchored_tag(bare, bare.index(primer), indexes.whitelist("i5"), params, a1, indexes.P5_END) \
        == (None, "no_index")
    assert anchored_tag(read, start, indexes.whitelist("i5"), params, a1, indexes.P5_END) == (i5, "found")
    # adapter missing: no tag rather than a guess
    assert anchored_tag("ACGT" * 15 + primer, 60, indexes.whitelist("i5"), params, a1) == (None, "no_adapter")


def test_nested_amplicons_are_not_chimeras(genome, tmp_path):
    g, chroms = genome
    x = chroms["1"][25000:25800]
    full = dict(chrom="1", start=25000, end=25800, strand="+", sequence=x,
                forward=bisulfite_primer(x[:24]), reverse=revcomp(bisulfite_primer(x[-24:])))
    short = dict(full, end=25300, sequence=x[:300], reverse=revcomp(bisulfite_primer(x[276:300])))
    samples = [Sample(0, "L", "", "LONG", "", "", "C57BL_6J", "C57BL_6J"),
               Sample(1, "S", "", "SHORT", "", "", "C57BL_6J", "C57BL_6J")]
    params = dict(DEFAULT_PARAMS)
    targets = [make_target("LONG", full, [], samples[:1], params)[0],
               make_target("SHORT", short, [], samples[1:], params)[0]]
    index = tmp_path / "idx.fa"
    index.write_text("".join(f">{t.name}\n{t.sequence.replace('C', 'T')}\n" for t in targets))
    caller = Caller(targets, params, index)

    def converted(sequence):
        return "".join("C" if b == "C" and sequence[i + 1:i + 2] == "G" else ("T" if b == "C" else b)
                       for i, b in enumerate(sequence))
    for sequence, expected in ((x, "LONG"), (x[:300], "SHORT")):
        for read in (converted(sequence), revcomp(converted(sequence))):
            record = caller.call("r", "ACGTTGCA" + read + "TTGACA", "", "f")
            assert (record["amplicon"], record["status"]) == (expected, "pass")


def test_parse_region_and_region_amplicon(genome, tmp_path):
    from ambs_mapper.reference import parse_region, region_amplicon
    assert parse_region("chr12:109,494,123-109,494,444") == ("12", 109494122, 109494444)
    with pytest.raises(ValueError):
        parse_region("chr12:500-100")
    g, chroms = genome
    product = region_amplicon(g.fasta.parent, "R", "chr2:5001-5600", "-")
    assert product["sequence"] == revcomp(chroms["2"][5000:5600]) and product["primer_lengths"] == (0, 0)


def test_genome_mapping_mode(genome, tmp_path, monkeypatch):
    from ambs_mapper import bsmap
    g, chroms = genome
    monkeypatch.setenv("AMBS_INDEX_DIR", str(tmp_path))
    mapper = bsmap.BisulfiteMapper(bsmap.build_index(g.fasta.parent, threads=2))

    def converted(sequence):
        return "".join("C" if b == "C" and sequence[i + 1:i + 2] == "G" else ("T" if b == "C" else b)
                       for i, b in enumerate(sequence))
    top = chroms["1"][20000:20450]
    product = dict(chrom="1", start=20000, end=20450, strand="+", sequence=top,
                   forward=bisulfite_primer(top[:24]), reverse=revcomp(bisulfite_primer(top[-24:])))
    samples = [Sample(0, "S", "", "A", "", "", "C57BL_6J", "C57BL_6J")]
    params = dict(DEFAULT_PARAMS)
    target = make_target("A", product, [], samples, params)[0]
    index = tmp_path / "amp.fa"
    index.write_text(f">A\n{top.replace('C', 'T')}\n")
    caller = Caller([target], params, index, mapper)
    on_target = converted(top)
    cases = {"OT": on_target, "CTOT": revcomp(on_target),
             "off": converted(chroms["2"][5000:5400]),
             "wrong": converted(revcomp(top))}          # the other strand of the same region
    record = caller.call("r", cases["OT"], "", "f")
    assert (record["status"], record["chrom"], record["bs_strand"], record["gstart"]) == ("pass", "1", "OT", 20001)
    assert caller.call("r", cases["CTOT"], "", "f")["bs_strand"] == "CTOT"
    assert caller.call("r", cases["off"], "", "f")["status"] == "off_target"
    assert caller.call("r", cases["wrong"], "", "f")["status"] == "wrong_strand"

    # BAM for IGV: a PCR copy (CTOT) is written like its original strand (OT):
    # forward flag, top-strand bases, so bisulfite colouring reads the top-strand Cs.
    import subprocess
    from ambs_mapper.pipeline import sam_line, write_bam
    body = tmp_path / "work" / "reads.body.sam"
    body.parent.mkdir()
    groups = {}
    with open(body, "w") as handle:
        for name, read in (("ot", cases["OT"]), ("ctot", cases["CTOT"]), ("wrong", cases["wrong"])):
            record = caller.call(name, read, "I" * len(read), "f")
            handle.write(sam_line(record, "S" if record["sample"] >= 0 else "", groups) + "\n")
    write_bam(tmp_path / "reads.bam", body, groups, g.fasta.parent, 2)
    assert (tmp_path / "reads.bam.bai").exists() and not body.exists()
    lines = subprocess.run(["samtools", "view", str(tmp_path / "reads.bam")], capture_output=True, text=True,
                           check=True).stdout.splitlines()
    fields = {line.split("\t")[0]: line.split("\t") for line in lines}
    assert fields["ot"][1:4] == fields["ctot"][1:4] == ["0", "1", "20001"]
    assert fields["ot"][9] == fields["ctot"][9] == on_target
    assert "YB:Z:CTOT" in fields["ctot"] and "RG:Z:S" in fields["ot"] and "YT:Z:pass" in fields["ot"]
    assert fields["wrong"][1] == "16" and "YT:Z:wrong_strand" in fields["wrong"]

    # one BAM per sample with its passed reads, named after the sample
    from ambs_mapper.pipeline import split_bam
    split = split_bam(tmp_path / "reads.bam", tmp_path / "bam_by_sample", groups, samples, ["A"])
    assert [(e["name"], e["amplicon"], e["unknown"]) for e in split] == [("S", "A", False)]
    assert not (tmp_path / "bam_by_sample" / "タグ不明_A.bam").exists()      # no unassigned reads here
    per_sample = tmp_path / "bam_by_sample" / "S_A.bam"
    names = subprocess.run(["samtools", "view", str(per_sample)], capture_output=True, text=True,
                           check=True).stdout.split()
    assert {"ot", "ctot"} <= set(names) and "wrong" not in names
    assert (tmp_path / "bam_by_sample.zip").exists() and not list((tmp_path / "bam_by_sample").glob("*.tmp.*"))


def test_clip_alignment():
    from ambs_mapper.bsmap import clip_alignment
    # 4S 10M 2I 5M 3D 20M 6S at 100, restricted to [105, 120)
    cigar = [(10, 0), (2, 1), (5, 0), (3, 2), (20, 0)]
    start, ops, left, right = clip_alignment(100, cigar, 4, 47, 105, 120)
    assert (start, ops, left, right) == (105, [(5, 0), (2, 1), (5, 0), (3, 2), (2, 0)], 9, 24)
    assert left + sum(n for n, op in ops if op in (0, 1)) + right == 47
    # an insertion just before the range start becomes part of the clip
    assert clip_alignment(100, [(5, 0), (3, 1), (10, 0)], 0, 18, 105, 115) == (105, [(10, 0)], 8, 0)
    assert clip_alignment(100, cigar, 4, 47, 500, 600) is None


def test_nextera_entry_tolerance():
    assert indexes.parse_index("CAGCCTCG", "i7") == indexes.parse_index("N710", "i7")   # written reverse-complemented
    base = dict(inputs=["x"], layout="nextera", amplicons=[dict(name="A", forward="ACGTACGTACGTAC", reverse="ACGTACGTACGTAC")],
                samples=[dict(name="s1", tag_f="GCTCATGA", tag_r="CTATTAAG", mother="B6", father="JF1")])
    request = validate_request(base)
    assert (request["samples"][0]["tag_f"], request["samples"][0]["tag_r"]) == ("CTATTAAG", "GCTCATGA")
    assert request["notes"]
    inline = validate_request({**base, "layout": "inline"})
    assert any("Nextera" in note for note in inline["notes"])


def test_methylation_tags_for_viewer():
    from ambs_mapper.igv_view import methylation_tags
    # top strand: CpG Cs at reference 100 (methylated, C) and 104 (unmethylated, T)
    fields = ["r", "0", "1", "101", "60", "2S8M", "*", "0", "0", "GGCATTTACA", "*"]   # read[2] at 100
    out = methylation_tags(list(fields), [100, 104], True)
    assert out[9] == "GGCATTCACA"                       # unmethylated T shown as C
    assert "MM:Z:C+m?,0,0;" in out and "ML:B:C,255,0" in out
    # bottom strand (flag 16): measured G/A on the top strand; MM counts Cs of the reverse-complemented read
    fields = ["r", "16", "1", "100", "60", "6M", "*", "0", "0", "AGTTGA", "*"]       # G at 100, A at 104
    out = methylation_tags(list(fields), [100, 104], False)
    assert out[9] == "AGTTGG"
    # reversed read CCAACT: the A->G site is its first C, the methylated G its third
    assert "MM:Z:C+m?,0,1;" in out and "ML:B:C,0,255" in out
