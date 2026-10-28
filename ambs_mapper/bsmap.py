"""Genome-wide bisulfite read mapping with minimap2 on a 3-letter genome.

The index holds two converted copies of GRCm39: "CT|chrom" (C->T, the top
strand after conversion) and "GA|chrom" (G->A, the bottom strand after
conversion, kept in top-strand coordinates). Each read is mapped twice, after
C->T and after G->A conversion; the four valid combinations give the strand
the read derives from:

    read C->T on CT+  : OT   (original top)
    read G->A on CT-  : CTOT (complement of the original top, PCR copy)
    read C->T on GA-  : OB   (original bottom)
    read G->A on GA+  : CTOB (complement of the original bottom)

OT/CTOT report methylation of top-strand CpGs, OB/CTOB of bottom-strand CpGs.
"""
from __future__ import annotations

import os
from pathlib import Path
import time

import mappy

from .progress import emit, log
from .reference import Genome, reference_fasta, revcomp

CT = str.maketrans("C", "T")
GA = str.maketrans("G", "A")
KMER = 19                       # 3-letter sequence needs longer k-mers than plain DNA
VALID = {("CT", 1, "ct"): "OT", ("CT", -1, "ga"): "CTOT", ("GA", -1, "ct"): "OB", ("GA", 1, "ga"): "CTOB"}
TOP = {"OT", "CTOT"}


def index_dir(resources):
    # Docker Desktop reads bind-mounted Windows folders slowly; a named volume
    # (AMBS_INDEX_DIR) loads the ~13 GB index several times faster.
    return Path(os.environ.get("AMBS_INDEX_DIR") or resources)


def index_path(resources, k=KMER):
    return index_dir(resources) / f"bisulfite_k{k}.mmi"


def index_ready(resources, k=KMER):
    return index_path(resources, k).exists()


def converted_fasta(resources):
    """Write the CT and GA genomes once; returns the path."""
    path = index_dir(resources) / "bisulfite.fa"
    if path.exists():
        return path
    genome = Genome(reference_fasta(resources))
    partial = path.with_suffix(".fa.part")
    with open(partial, "w", encoding="ascii") as handle:
        for label, table in (("CT", CT), ("GA", GA)):
            for n, chrom in enumerate(genome.index, 1):
                emit("bisulfite_index", message=f"変換ゲノムを作成しています（{label} {chrom}）。",
                     done=n, total=len(genome.index))
                sequence = genome.chromosome(chrom).decode("ascii").upper().translate(table)
                handle.write(f">{label}|{chrom}\n")
                for i in range(0, len(sequence), 80):
                    handle.write(sequence[i:i + 80] + "\n")
    partial.rename(path)
    return path


def build_index(resources, k=KMER, threads=8):
    """minimap2 index of the converted genome (about 5 min with 8 threads, ~13 GB)."""
    target = index_path(resources, k)
    if target.exists():
        return target
    fasta = converted_fasta(resources)
    partial = target.with_suffix(".mmi.part")
    emit("bisulfite_index", message="バイサルファイト用ゲノム索引を作成しています（初回のみ、5分程度）。")
    started = time.time()
    aligner = mappy.Aligner(str(fasta), preset="map-ont", k=k, n_threads=threads, fn_idx_out=str(partial))
    if not aligner:
        raise RuntimeError("ゲノム索引を作成できませんでした。")
    del aligner
    partial.rename(target)
    fasta.unlink()              # only needed to build the index
    log(f"ゲノム索引を作成しました（{time.time() - started:.0f}秒）。")
    return target


class Hit:
    """One read placement on the genome (top-strand coordinates, 0-based half-open)."""
    __slots__ = ("chrom", "start", "end", "strand", "mapq", "q_st", "q_en", "mlen", "blen", "primary",
                 "cigar", "orientation")

    def __init__(self, chrom, start, end, strand, mapq, q_st, q_en, mlen, blen, primary, cigar=(), orientation=1):
        self.chrom, self.start, self.end, self.strand = chrom, start, end, strand
        self.mapq, self.q_st, self.q_en, self.mlen, self.blen, self.primary = mapq, q_st, q_en, mlen, blen, primary
        self.cigar, self.orientation = cigar, orientation   # minimap2 CIGAR (reference order), read strand

    @property
    def top(self):
        return self.strand in TOP


class BisulfiteMapper:
    def __init__(self, index, threads=1):
        self.aligner = mappy.Aligner(str(index), preset="map-ont", best_n=5, n_threads=threads)
        if not self.aligner:
            raise RuntimeError("ゲノム索引を読み込めませんでした。")

    def map(self, sequence):
        """Valid placements, best first. Query coordinates refer to the read as given."""
        sequence = sequence.upper()
        hits = []
        for kind, query in (("ct", sequence.translate(CT)), ("ga", sequence.translate(GA))):
            for h in self.aligner.map(query, cs=False):
                label, chrom = h.ctg.split("|", 1)
                strand = VALID.get((label, h.strand, kind))
                if strand:
                    hits.append(Hit(chrom, h.r_st, h.r_en, strand, h.mapq, h.q_st, h.q_en, h.mlen, h.blen,
                                    h.is_primary, h.cigar, h.strand))
        hits.sort(key=lambda h: -h.mlen)
        return hits


CIGAR_OPS = "MIDNSHP=X"
CONVERSIONS = {"OT": ("CT", "CT"), "CTOT": ("GA", "CT"), "OB": ("CT", "GA"), "CTOB": ("GA", "GA")}


def clip_alignment(start, cigar, left, length, lo, hi):
    """Restrict an alignment to reference [lo, hi): (start, ops, left clip, right clip).

    cigar is in reference order without clips; left is the number of read bases
    (in SEQ order) before the alignment. Bases outside the range become soft clips.
    """
    ops, r, q = [], start, left
    first = None
    for n, op in cigar:
        if op in (0, 7, 8):                         # M, =, X
            a, b = max(r, lo), min(r + n, hi)
            if a < b:
                if first is None:
                    first, new_start = q + (a - r), a
                ops.append((b - a, op))
            r, q = r + n, q + n
        elif op in (2, 3):                          # D, N: only between kept bases
            if first is not None and r + n <= hi:
                ops.append((n, op))
            r += n
        else:                                       # I
            if first is not None and r < hi:
                ops.append((n, op))
            q += n
    while ops and ops[-1][1] not in (0, 7, 8):      # no indel next to a clip
        ops.pop()
    if first is None or not ops:
        return None
    end = first + sum(n for n, op in ops if op in (0, 1, 7, 8))
    return new_start, ops, first, length - end


def sam_record(name, sequence, quality, hit, clip=None):
    """SAM line (without sample tags) for a read placed on the unconverted reference.

    SEQ holds the original read bases, so C/T at CpGs show methylation. The
    reverse flag follows the bisulfite strand (top: forward, bottom: reverse)
    rather than the read direction, so PCR copies (CTOT/CTOB) are coloured
    correctly by IGV's bisulfite mode and read as OT/OB by tools that assume a
    directional library. YB keeps the original strand; XR/XG are as in Bismark.
    clip=(lo, hi) soft-clips everything aligned outside that reference range.
    """
    length = len(sequence)
    if hit.orientation == 1:
        seq, qual, left, right = sequence, quality, hit.q_st, length - hit.q_en
    else:
        seq, qual, left, right = revcomp(sequence), quality[::-1], length - hit.q_en, hit.q_st
    start, ops = hit.start, list(hit.cigar)
    clipped = clip_alignment(start, ops, left, length, *clip) if clip else None
    if clipped:
        start, ops, left, right = clipped
    cigar = ((f"{left}S" if left else "") + "".join(f"{n}{CIGAR_OPS[op]}" for n, op in ops)
             + (f"{right}S" if right else ""))
    flag = 0 if hit.top else 16
    read_conversion, genome_conversion = CONVERSIONS[hit.strand]
    return "\t".join([name, str(flag), hit.chrom, str(start + 1), str(hit.mapq), cigar, "*", "0", "0", seq,
                      qual if qual and len(qual) == length else "*",
                      f"YB:Z:{hit.strand}", f"XR:Z:{read_conversion}", f"XG:Z:{genome_conversion}"])
