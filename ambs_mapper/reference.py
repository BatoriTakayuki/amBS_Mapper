"""GRCm39 access and bisulfite-aware in-silico PCR.

Coordinates are 0-based half-open on the top strand unless stated otherwise.
An amplicon is described on its "bisulfite strand" X: the original genomic
strand that the forward primer anneals to after conversion, written 5'->3'.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import urllib.request

from .progress import emit

FASTA_URL = ("https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/"
             "Mus_musculus.GRCm39.dna.primary_assembly.fa.gz")
ASSEMBLY = "GRCm39"

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT", "S": "CG",
         "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
_COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def revcomp(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def clean_primer(primer: str) -> str:
    primer = re.sub(r"\s+", "", str(primer or "")).upper().replace("U", "T")
    if not primer:
        raise ValueError("プライマー配列が空です。")
    if any(base not in IUPAC for base in primer):
        raise ValueError(f"プライマーに使えない文字があります: {primer}")
    if len(primer) < 12:
        raise ValueError(f"プライマーが短すぎます（12塩基以上）: {primer}")
    return primer


class Genome:
    """Random access to an uncompressed FASTA through its .fai index."""

    def __init__(self, fasta):
        self.fasta = Path(fasta)
        self.index = {}
        with open(str(self.fasta) + ".fai", encoding="utf-8") as handle:
            for line in handle:
                name, length, offset, bases, width = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(bases), int(width))

    def length(self, chrom):
        return self.index[chrom][0]

    def _read(self, chrom, start, end) -> bytes:
        length, offset, bases, width = self.index[chrom]
        start, end = max(0, start), min(length, end)
        if start >= end:
            return b""
        first = offset + (start // bases) * width + start % bases
        last = offset + ((end - 1) // bases) * width + (end - 1) % bases
        with open(self.fasta, "rb") as handle:
            handle.seek(first)
            raw = handle.read(last - first + 1)
        return raw.replace(b"\n", b"").replace(b"\r", b"").upper()

    def fetch(self, chrom, start, end) -> str:
        return self._read(chrom, start, end).decode("ascii")

    def chromosome(self, chrom) -> bytes:
        return self._read(chrom, 0, self.length(chrom))


def _pattern(primer: str, source: str, target: str) -> bytes:
    """Regex for a primer in a fully converted (3-letter) sequence."""
    parts = []
    for base in primer:
        options = sorted(set(IUPAC[base].replace(source, target)))
        parts.append(options[0] if len(options) == 1 else "[" + "".join(options) + "]")
    return "".join(parts).encode("ascii")


def _finditer(pattern: bytes, sequence: bytes):
    if re.fullmatch(rb"[ACGT]+", pattern):
        position = sequence.find(pattern)
        while position >= 0:
            yield position
            position = sequence.find(pattern, position + 1)
    else:
        for match in re.finditer(b"(?=" + pattern + b")", sequence):
            yield match.start()


_CT = bytes.maketrans(b"C", b"T")
_GA = bytes.maketrans(b"G", b"A")


def in_silico_pcr(genome: Genome, forward: str, reverse: str, max_size: int, chromosomes=None):
    """Find bisulfite PCR products of forward/reverse primers on both converted strands.

    On the top strand the forward primer matches C->T converted sequence and the
    reverse complement of the reverse primer matches downstream. On the bottom
    strand the same relation holds for the G->A converted top sequence read in
    reverse.
    """
    f_top, rc_r_top = _pattern(forward, "C", "T"), _pattern(revcomp(reverse), "C", "T")
    rc_f_bottom, r_bottom = _pattern(revcomp(forward), "G", "A"), _pattern(reverse, "G", "A")
    lf, lr = len(forward), len(reverse)
    products = []
    names = chromosomes or list(genome.index)
    for number, chrom in enumerate(names, 1):
        emit("primer_search", chrom=chrom, done=number - 1, total=len(names))
        sequence = genome.chromosome(chrom)
        top = sequence.translate(_CT)
        starts = list(_finditer(f_top, top))
        if starts:
            ends = sorted(position + lr for position in _finditer(rc_r_top, top))
            for start in starts:
                for end in ends:
                    if end - start > max_size:
                        break
                    if end - start >= lf + lr:
                        products.append(dict(chrom=chrom, start=start, end=end, strand="+"))
        del top
        bottom = sequence.translate(_GA)
        ends = list(_finditer(rc_f_bottom, bottom))
        if ends:
            starts = sorted(_finditer(r_bottom, bottom))
            for f_site in ends:
                end = f_site + lf
                for start in reversed(starts):
                    if end - start > max_size:
                        break
                    if end - start >= lf + lr and start <= f_site:
                        products.append(dict(chrom=chrom, start=start, end=end, strand="-"))
        del bottom, sequence
    for product in products:
        product["size"] = product["end"] - product["start"]
    return products


def amplicon_sequence(genome: Genome, product) -> str:
    """Unconverted reference sequence of the bisulfite strand X, 5'->3'."""
    sequence = genome.fetch(product["chrom"], product["start"], product["end"])
    return sequence if product["strand"] == "+" else revcomp(sequence)


def x_to_genome(product, x: int) -> int:
    """0-based top-strand coordinate of position x on the bisulfite strand."""
    return product["start"] + x if product["strand"] == "+" else product["end"] - 1 - x


def find_amplicon(resources, name, forward, reverse, max_size):
    """Resolve primers to a unique product, caching results under resources/amplicons."""
    forward, reverse = clean_primer(forward), clean_primer(reverse)
    max_size = int(max_size)
    key = hashlib.sha1(f"{ASSEMBLY}|{forward}|{reverse}|{max_size}".encode()).hexdigest()[:16]
    cache = Path(resources) / "amplicons" / f"{key}.json"
    if cache.is_file():
        products = json.loads(cache.read_text(encoding="utf-8"))
    else:
        genome = Genome(reference_fasta(resources))
        products = in_silico_pcr(genome, forward, reverse, max_size)
        for product in products:
            product["swapped"] = False
        if not products:
            # Forward and reverse columns are a common mix-up; accept the swapped pair with a note.
            products = in_silico_pcr(genome, reverse, forward, max_size)
            for product in products:
                product["swapped"] = True
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(products), encoding="utf-8")
    if not products:
        raise ValueError(f"アンプリコン{name}: プライマーに一致する{max_size:,} bp以下のPCR産物が"
                         f"{ASSEMBLY}の変換配列上に見つかりません。配列と最大産物長を確認してください。")
    if len(products) > 1:
        listed = ", ".join(f"{p['chrom']}:{p['start'] + 1}-{p['end']}({p['strand']})" for p in products[:5])
        raise ValueError(f"アンプリコン{name}: PCR産物が{len(products)}か所で予測されました（{listed}）。"
                         "一意になるようにプライマーか最大産物長を見直してください。")
    product = dict(products[0])
    genome = Genome(reference_fasta(resources))
    if product["swapped"]:
        forward, reverse = reverse, forward
    product.update(name=name, forward=forward, reverse=reverse, assembly=ASSEMBLY,
                   sequence=amplicon_sequence(genome, product))
    return product


REGION = re.compile(r"(?:chr)?([0-9A-Za-z_.]+):([0-9,]+)-([0-9,]+)")
PSEUDO_PRIMER = 20


def parse_region(text):
    """'chr12:109,494,123-109,494,444' -> ('12', 109494122, 109494444) (0-based half-open)."""
    match = REGION.fullmatch(re.sub(r"\s+", "", str(text or "")))
    if not match:
        raise ValueError(f"領域「{text}」は chr12:109494123-109494444 の形式で入力してください。")
    chrom, start, end = match.group(1), int(match.group(2).replace(",", "")), int(match.group(3).replace(",", ""))
    if not 1 <= start < end:
        raise ValueError(f"領域「{text}」の開始と終了を確認してください。")
    return chrom, start - 1, end


def _bisulfite_primer(sequence):
    # C->T except CpG C, which may stay C (Y), as a primer for strand X would be designed
    return "".join("Y" if base == "C" and sequence[i + 1:i + 2] == "G" else "T" if base == "C" else base
                   for i, base in enumerate(sequence))


def region_amplicon(resources, name, region, strand, forward="", reverse=""):
    """Product for a region given directly (no primer search). Missing primers are
    replaced by the converted ends of the region for read anchoring only."""
    chrom, start, end = parse_region(region)
    if strand not in ("+", "-"):
        raise ValueError(f"アンプリコン{name}: 鎖は + か - で指定してください。")
    genome = Genome(reference_fasta(resources))
    if chrom not in genome.index:
        raise ValueError(f"アンプリコン{name}: 染色体「{chrom}」が{ASSEMBLY}にありません（例: 12, X, MT）。")
    if end > genome.length(chrom) or end - start > 200000:
        raise ValueError(f"アンプリコン{name}: 領域が染色体の範囲外か、200 kbを超えています。")
    product = dict(chrom=chrom, start=start, end=end, strand=strand, size=end - start, swapped=False)
    sequence = amplicon_sequence(genome, product)
    forward, reverse = (clean_primer(p) if str(p or "").strip() else "" for p in (forward, reverse))
    lengths = (len(forward), len(reverse))
    forward = forward or _bisulfite_primer(sequence[:PSEUDO_PRIMER])
    reverse = reverse or revcomp(_bisulfite_primer(sequence[-PSEUDO_PRIMER:]))
    product.update(name=name, forward=forward, reverse=reverse, assembly=ASSEMBLY, sequence=sequence,
                   primer_lengths=lengths)
    return product


def reference_fasta(resources) -> Path:
    path = Path(resources) / "reference.fa"
    if not path.is_file() or not Path(str(path) + ".fai").is_file():
        raise ValueError("参照ゲノム(GRCm39)が未準備です。画面の「参照ゲノムを準備」を実行してください。")
    return path


def reference_ready(resources) -> bool:
    try:
        reference_fasta(resources)
        return True
    except ValueError:
        return False


def setup_reference(resources, url=FASTA_URL):
    """Download GRCm39 primary assembly and build the FASTA index."""
    resources = Path(resources)
    resources.mkdir(parents=True, exist_ok=True)
    target = resources / "reference.fa"
    if reference_ready(resources):
        emit("reference", message="参照ゲノムは準備済みです。")
        return target
    archive = resources / "reference.fa.gz.partial"
    emit("reference", message="GRCm39をダウンロードしています（約800 MB）。")
    with urllib.request.urlopen(url, timeout=60) as response, open(archive, "wb") as handle:
        total = int(response.headers.get("Content-Length") or 0)
        done = 0
        while chunk := response.read(8 * 1024 * 1024):
            handle.write(chunk)
            done += len(chunk)
            emit("download", done=done, total=total)
    emit("reference", message="展開しています。")
    partial = resources / "reference.fa.partial"
    with gzip.open(archive, "rb") as source, open(partial, "wb") as handle:
        shutil.copyfileobj(source, handle, length=16 * 1024 * 1024)
    subprocess.run(["samtools", "faidx", str(partial)], check=True)
    Path(str(partial) + ".fai").replace(str(target) + ".fai")
    partial.replace(target)
    archive.unlink()
    emit("reference", message="参照ゲノムの準備が完了しました。")
    return target
