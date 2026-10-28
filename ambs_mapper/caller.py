"""Per-read analysis: 3-letter alignment, inline tags, QC, allele and CpG calls.

Reads are aligned with minimap2 (map-ont) after converting C->T in both the
read (each orientation) and the amplicon's bisulfite strand X, as in Bismark.
The base-level calls then use the original read bases on that alignment.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import edlib
import mappy

from . import indexes
from .bsmap import sam_record
from .reference import IUPAC, revcomp

CT = str.maketrans("C", "T")

STATUS_LABELS = {
    "pass": "合格",
    "unmapped": "ゲノムにマッピングされない",
    "off_target": "アンプリコン外にマッピング",
    "wrong_strand": "反対鎖の産物（ミスプライミング等）",
    "chimera": "キメラ（複数アンプリコン）",
    "tag_unknown": "タグ不明",
    "tag_ambiguous": "タグが複数サンプルに一致",
    "tag_minor": "リードが少ないインデックス組み合わせ",
    "low_coverage": "アンプリコン被覆率が低い",
    "low_identity": "一致率が低い",
    "high_mismatch": "ミスマッチが多い",
    "low_conversion": "変換率が低い",
    "unconverted_cph": "未変換CpHが多い",
}
ALLELE_LABELS = {"mother": "母アレル", "father": "父アレル", "unassigned": "判定不能", "same": "同一系統"}

DEFAULT_PARAMS = dict(
    min_coverage=90.0,       # % of the amplicon spanned by the alignment
    min_identity=90.0,       # QUMA -p: lower limit of percent identity
    max_mismatch=10.0,       # QUMA -m: upper limit of mismatches (%)
    min_conversion=90.0,     # QUMA -c: lower limit of converted CpH (%)
    max_unconverted=6,       # QUMA -u: upper limit of unconverted CpH (count); None disables
    tag_max_errors=2,        # edit distance allowed for an inline tag
    tag_slack=3,             # bases allowed between tag and primer (spacer or indels)
    snp_min_quality=10,      # base quality required at a SNP for allele votes
    allele_min_fraction=80.0,  # % of SNP votes that must agree with the assigned allele
    exclude_primer_cpg=True,
    photo_rows=200,
    min_length=0,
    auto_min_reads=50,       # automatic samples: index combinations with fewer passed reads are not reported
    single_index=True,       # Nextera: assign reads carrying only i5 or only i7 when that index is unique
)


@dataclass
class Sample:
    index: int
    name: str
    group: str
    amplicon: str
    tag_f: str
    tag_r: str
    mother: str
    father: str
    fastq: str = ""
    labels: tuple = ("", "")        # display names of the tags (Nextera index names)

    @property
    def pair(self):
        return (self.mother, self.father)


@dataclass
class Genotype:
    """Allele-specific information for one amplicon and one mother/father pair."""
    snps: list                      # informative SNPs: (x, mother_base, father_base)
    mother_sites: list              # bool per CpG site
    father_sites: list
    mother3: str = ""               # C->T converted allele sequences on X
    father3: str = ""


@dataclass
class Target:
    name: str
    sequence: str                   # bisulfite strand X, unconverted, reference (C57BL/6J)
    forward: str
    reverse: str
    sites: list                     # x of the C of every CpG analysed (union over strains)
    conversion: list                # x of non-CpG C used for conversion rate
    variable: set                   # x of any SNP between analysed strains
    variants: dict                  # x -> set of bases seen in analysed strains
    samples: list = field(default_factory=list)
    genotypes: dict = field(default_factory=dict)
    tags_f: list = field(default_factory=list)
    tags_r: list = field(default_factory=list)
    layout: str = "inline"          # "inline" or "nextera"
    auto_pair: tuple = None         # (mother, father) when samples are created from the tags found
    chrom: str = ""                 # genomic region (0-based half-open, top-strand coordinates)
    start: int = 0
    end: int = 0
    strand: str = "+"               # "+": X is the top strand, "-": the bottom strand


def _equalities(sequence):
    pairs = []
    for letter in set(sequence):
        if len(IUPAC.get(letter, letter)) > 1:
            pairs.extend((letter, base) for base in IUPAC[letter])
    return pairs


def find_primer(read, primer, approx_start, margin=40, limit=None):
    """Start of the primer in the read near approx_start, or None.

    approx_start comes from the alignment and is off when the read does not span
    the amplicon (e.g. a primer-to-primer by-product); then the primer is looked
    for anywhere in read[:limit], the part 5' of the aligned stretch.
    """
    if approx_start is None:
        return None
    windows = [max(0, approx_start - margin), max(0, approx_start) + len(primer) + margin]
    for lo, hi in (windows, (0, limit)) if limit else (windows,):
        window = read[lo:hi]
        if len(window) < len(primer) // 2:
            continue
        hit = edlib.align(primer, window, mode="HW", task="locations", k=max(2, len(primer) // 4),
                          additionalEqualities=_equalities(primer))
        if hit["editDistance"] >= 0:
            start, end = hit["locations"][0]
            return lo + (start if start is not None else max(0, end + 1 - len(primer)))
    return None


def _upstream(read, primer_start, length):
    # The tag's 3' end is anchored at the primer (or adapter); reversing both sequences
    # lets the prefix mode (SHW) keep that end fixed while the 5' end stays free.
    return read[max(0, primer_start - length):primer_start][::-1]


def spacer_errors(read, primer_start, spacer):
    """Edit distance of an adapter placed directly 5' of the primer (large if absent)."""
    if primer_start is None:
        return 999
    region = _upstream(read, primer_start, len(spacer) + 10)
    return edlib.align(spacer[::-1], region, mode="SHW")["editDistance"] if region else 999


def anchored_tag(read, primer_start, tags, params, spacer="", outer=""):
    """Tag whose 3' end meets the primer, optionally through a fixed adapter (spacer).

    Returns (tag, reason): tag is None when nothing matches and "" when two tags
    tie; reason then says why (see TAG_REASONS).
    Every candidate carries the same spacer, so the spacer's own errors cancel out
    in the comparison and are subtracted before applying the tag error limit.
    outer is constant sequence 5' of the tag (P5/P7); requiring it rejects reads
    whose index is missing, where a sequencing adapter can mimic a short tag.
    """
    if not tags or primer_start is None:
        return None, "no_primer"
    longest = max(len(tag) for tag in tags)
    region = _upstream(read, primer_start, longest + len(spacer) + len(outer) + params["tag_slack"])
    if not region:
        return None, "no_adapter" if spacer else "unreadable"
    base = 0
    if spacer:
        base = edlib.align(spacer[::-1], region, mode="SHW")["editDistance"]
        if base > max(3, len(spacer) // 5):
            return None, "no_adapter"
    if outer and edlib.align(outer[::-1], region, mode="HW", k=len(outer) // 4)["editDistance"] < 0:
        return None, "no_index"         # the molecule ends at the adapter: no index was attached
    scores = sorted((edlib.align((tag + spacer)[::-1], region, mode="SHW")["editDistance"] - base, tag)
                    for tag in tags)
    best, tag = scores[0]
    if best > params["tag_max_errors"]:
        return None, "unreadable"
    if outer:
        inner = best + base
        if edlib.align((outer + tag + spacer)[::-1], region, mode="SHW")["editDistance"] - inner > len(outer) // 4:
            return None, "unreadable"
    if len(scores) > 1 and scores[1][0] == best:
        return "", "ambiguous"
    return tag, "found"


def best_tag(read, primer, approx_start, tags, params, spacer=""):
    """Tag 5' of the primer, or None. approx_start is the expected primer start."""
    if not tags:
        return None, None
    return anchored_tag(read, find_primer(read, primer, approx_start), tags, params, spacer)


SNP_WINDOW = 8


def snp_vote(read3, quality, mapping, x, genotype, hit, min_quality):
    """Which parental sequence the read matches better around SNP x.

    A single alignment column is unreliable next to homopolymers, where Nanopore
    indels let the aligner shift the SNP base. Comparing a short window against
    both converted allele sequences is independent of gap placement.
    """
    lo, hi = x - SNP_WINDOW, x + SNP_WINDOW + 1
    if lo < hit.r_st or hi > hit.r_en:
        return None
    if quality and mapping[x] >= 0 and ord(quality[mapping[x]]) - 33 < min_quality:
        return None
    left = next((mapping[i] for i in range(lo, x + 1) if mapping[i] >= 0), -1)
    right = next((mapping[i] for i in range(hi - 1, x - 1, -1) if mapping[i] >= 0), -1)
    if left < 0 or right < left:
        return None
    segment = read3[max(0, left - 3):right + 4]
    d_m = edlib.align(genotype.mother3[lo:hi], segment, mode="HW")["editDistance"]
    d_f = edlib.align(genotype.father3[lo:hi], segment, mode="HW")["editDistance"]
    if d_m == d_f:
        return None
    return "mother" if d_m < d_f else "father"


def reference_map(hit, length):
    """ref position -> read position (or -1) and counts of inserted/deleted bases."""
    mapping = [-1] * length
    qi, ri, insertions, deletions = hit.q_st, hit.r_st, 0, 0
    for size, op in hit.cigar:
        if op in (0, 7, 8):
            for offset in range(size):
                mapping[ri + offset] = qi + offset
            qi += size
            ri += size
        elif op == 1:
            qi += size
            insertions += size
        elif op in (2, 3):
            ri += size
            deletions += size
    return mapping, insertions, deletions


class Caller:
    def __init__(self, targets, params, index_path, mapper=None):
        self.targets = {target.name: target for target in targets}
        self.params = params
        self.mapper = mapper            # bsmap.BisulfiteMapper; genome-wide placement of every read
        self.aligner = mappy.Aligner(str(index_path), preset="map-ont", best_n=5)
        if not self.aligner:
            raise RuntimeError("アンプリコン参照のインデックスを作成できませんでした。")

    def _align(self, sequence, allowed=None):
        forward3, reverse_sequence = sequence.translate(CT), revcomp(sequence)
        candidates = []
        for orientation, query in ((1, forward3), (-1, reverse_sequence.translate(CT))):
            for hit in self.aligner.map(query):
                if hit.strand == 1 and (allowed is None or hit.ctg in allowed):
                    candidates.append((hit.mlen, orientation, hit))
        return candidates, reverse_sequence

    def _place(self, sequence, record, placed):
        """Genome-wide placement: names of the targets the read overlaps, or None with the status set.
        The chosen genome hit is left in placed["hit"] for the BAM record."""
        hits = self.mapper.map(sequence)
        if not hits:
            return None

        def overlapping(hit):
            return [t for t in self.targets.values()
                    if t.chrom == hit.chrom and min(t.end, hit.end) > max(t.start, hit.start)]
        # Nanopore reads may join two library molecules (e.g. an amplicon and an
        # unrelated product). Tags are read next to the amplicon's own primers, so
        # the amplicon part stays usable: any placement on a target counts.
        def shadowed(hit):
            # an alternative placement of the same read segment that a better hit explains
            return any(other.mlen > hit.mlen and min(other.q_en, hit.q_en) - max(other.q_st, hit.q_st)
                       > 0.5 * (hit.q_en - hit.q_st) for other in hits)
        on_target = [hit for hit in hits if overlapping(hit) and not shadowed(hit)]
        hit = next((h for h in on_target if any((t.strand == "+") == h.top for t in overlapping(h))),
                   on_target[0] if on_target else hits[0])
        record.update(chrom=hit.chrom, gstart=hit.start + 1, gend=hit.end, bs_strand=hit.strand, mapq=hit.mapq)
        placed["hit"] = hit
        if not on_target:
            record["status"] = "off_target"
            return None
        same_strand = [t.name for t in overlapping(hit) if (t.strand == "+") == hit.top]
        if not same_strand:
            record["status"] = "wrong_strand"
            return None
        return set(same_strand)

    def _tags(self, target, read, hit, length):
        """((tag_f, reason), (tag_r, reason)) on the read oriented to strand X; see anchored_tag."""
        p = self.params
        if not target.tags_f and not target.tags_r:
            return (None, ""), (None, "")
        rc_read = revcomp(read)
        # The fallback windows end where the aligned part starts (plus the primer length).
        f_start = find_primer(read, target.forward, hit.q_st - hit.r_st, limit=hit.q_st + len(target.forward))
        approx = (len(read) - hit.q_en) - (length - hit.r_en)
        r_start = find_primer(rc_read, target.reverse, approx, limit=len(read) - hit.q_en + len(target.reverse))
        if target.layout != "nextera":
            # Inline tags sit 5' of each primer: on the oriented read for the forward
            # primer, and on its reverse complement for the reverse primer.
            return (anchored_tag(read, f_start, target.tags_f, p) if target.tags_f else (None, ""),
                    anchored_tag(rc_read, r_start, target.tags_r, p) if target.tags_r else (None, ""))
        # Nextera: i5 sits behind adapter A1 and i7 behind A2. Which gene-specific
        # primer carries A1 is read from the adapters themselves, so swapped
        # Forward/Reverse entries still give the right indexes.
        a1_forward = (spacer_errors(read, f_start, indexes.A1) + spacer_errors(rc_read, r_start, indexes.A2)
                      <= spacer_errors(read, f_start, indexes.A2) + spacer_errors(rc_read, r_start, indexes.A1))
        i5_side, i7_side = ((read, f_start), (rc_read, r_start)) if a1_forward else ((rc_read, r_start), (read, f_start))
        return (anchored_tag(*i5_side, target.tags_f, p, indexes.A1, indexes.P5_END),
                anchored_tag(*i7_side, target.tags_r, p, indexes.A2, indexes.P7_END))

    def call(self, name, sequence, quality, source):
        placed = {}
        sequence = sequence.upper()
        record = self._call(name, sequence, quality, source, placed)
        hit = placed.get("hit")
        if hit is not None:
            # On an amplicon, bases aligned beyond its ends are adapter or tag
            # sequence that the aligner stretched over the genome; IGV shows them clipped.
            target = self.targets.get(record["amplicon"])
            clip = (target.start, target.end) if target is not None and target.chrom == hit.chrom else None
            record["sam"] = sam_record(name, sequence, quality, hit, clip)
        return record

    def _call(self, name, sequence, quality, source, placed):
        p = self.params
        record = dict(read=name, source=source, length=len(sequence), amplicon="", orientation="",
                      sample=-1, tag_f="", tag_r="", tag_sides=0, status="unmapped", coverage=0.0, identity=0.0,
                      mismatch=0.0, conversion=0.0, unconverted=0, votes_m=0, votes_f=0, allele="", pattern="",
                      chrom="", gstart=0, gend=0, bs_strand="", mapq=0, tag_f_label="", tag_r_label="")
        if len(sequence) < max(30, p["min_length"]):
            return record
        allowed = None
        if self.mapper is not None:
            allowed = self._place(sequence, record, placed)
            if allowed is None:
                return record
        candidates, reverse_sequence = self._align(sequence, allowed)
        if not candidates:
            if allowed:
                # Placed in the amplicon's region but too little of it to align locally.
                record.update(amplicon=sorted(allowed)[0], status="low_coverage")
            return record
        def fraction(item):
            return (item[2].r_en - item[2].r_st) / len(self.targets[item[2].ctg].sequence)

        def read_span(item):
            hit = item[2]
            return (hit.q_st, hit.q_en) if item[1] == 1 else (len(sequence) - hit.q_en, len(sequence) - hit.q_st)

        # Nested amplicons (a short product inside a long one, sharing a primer) give
        # near-equal alignments; prefer the amplicon the read covers end to end.
        best_mlen = max(item[0] for item in candidates)
        _, orientation, hit = max((item for item in candidates if item[0] >= 0.9 * best_mlen), key=fraction)
        target = self.targets[hit.ctg]
        length = len(target.sequence)
        read = sequence if orientation == 1 else reverse_sequence
        if quality:
            quality = quality if orientation == 1 else quality[::-1]
        record.update(amplicon=target.name, orientation="+" if orientation == 1 else "-",
                      coverage=round(100.0 * (hit.r_en - hit.r_st) / length, 1))
        # A chimera holds two amplicons in different parts of the read; overlapping
        # hits are the same molecule matching nested or overlapping amplicons.
        long_hits = [read_span(item) for item in candidates if fraction(item) >= 0.5]
        for i, (s1, e1) in enumerate(long_hits):
            for s2, e2 in long_hits[i + 1:]:
                if min(e1, e2) - max(s1, s2) < 0.5 * min(e1 - s1, e2 - s2):
                    record["status"] = "chimera"
                    return record

        (tag_f, why_f), (tag_r, why_r) = self._tags(target, read, hit, length)
        # Nanopore also reads molecules indexed on one side only (one-sided
        # extension in the index PCR); the missing index (None) may then be inferred.
        single = p["single_index"] and target.layout == "nextera" and (tag_f is None) != (tag_r is None)

        def label(tag, why, kind):
            if not tag:
                return why                  # no_index, unreadable, ambiguous, no_adapter, no_primer or ""
            return indexes.label(tag, kind) if target.layout == "nextera" else tag
        record.update(tag_f=tag_f or "", tag_r=tag_r or "", tag_sides=(tag_f is not None) + (tag_r is not None),
                      tag_f_label=label(tag_f, why_f, "i5"), tag_r_label=label(tag_r, why_r, "i7"))
        if tag_f == "" or tag_r == "":
            record["status"] = "tag_ambiguous"
            return record
        if target.auto_pair:
            # Samples are created from the tag combinations found (in the main process);
            # single-index reads are resolved there once the combinations are known.
            if target.layout == "nextera" and not (tag_f and tag_r) and not single:
                record["status"] = "tag_unknown"
                return record
            record["auto_key"] = (tag_f or "", tag_r or "")
            pair = target.auto_pair
        else:
            matches = [sample for sample in target.samples
                       if (not sample.fastq or source == sample.fastq
                           or source.startswith(sample.fastq.rstrip("/") + "/"))
                       and (not sample.tag_f or sample.tag_f == tag_f or (single and tag_f is None))
                       and (not sample.tag_r or sample.tag_r == tag_r or (single and tag_r is None))]
            if len(matches) != 1:
                record["status"] = "tag_ambiguous" if len(matches) > 1 else "tag_unknown"
                return record
            record["sample"] = matches[0].index
            pair = matches[0].pair

        mapping, insertions, deletions = reference_map(hit, length)
        reference = target.sequence
        matched = mismatched = 0
        for x in range(hit.r_st, hit.r_en):
            position = mapping[x]
            if position < 0:
                continue
            base, expected = read[position], reference[x]
            if base == expected or (expected == "C" and base == "T"):
                matched += 1
            elif x in target.variants and (base in target.variants[x]
                                           or ("C" in target.variants[x] and base == "T")):
                matched += 1
            else:
                mismatched += 1
        aligned = matched + mismatched + insertions + deletions
        span = hit.r_en - hit.r_st
        converted = unconverted = 0
        for x in target.conversion:
            if hit.r_st <= x < hit.r_en and mapping[x] >= 0:
                base = read[mapping[x]]
                if base == "T":
                    converted += 1
                elif base == "C":
                    unconverted += 1
        conversion = 100.0 * converted / (converted + unconverted) if converted + unconverted else 0.0
        record.update(identity=round(100.0 * matched / aligned, 1) if aligned else 0.0,
                      mismatch=round(100.0 * mismatched / span, 1) if span else 0.0,
                      conversion=round(conversion, 1), unconverted=unconverted)

        genotype = target.genotypes[pair]
        votes_m = votes_f = 0
        read3 = read.translate(CT)
        for x, _, _ in genotype.snps:
            vote = snp_vote(read3, quality, mapping, x, genotype, hit, p["snp_min_quality"])
            votes_m += vote == "mother"
            votes_f += vote == "father"
        if pair[0] == pair[1]:
            allele = "same"
        elif votes_m + votes_f and max(votes_m, votes_f) * 100.0 >= p["allele_min_fraction"] * (votes_m + votes_f):
            allele = "mother" if votes_m > votes_f else "father"
        else:
            allele = "unassigned"
        present = (genotype.mother_sites if allele == "mother" else genotype.father_sites if allele == "father"
                   else [m and f for m, f in zip(genotype.mother_sites, genotype.father_sites)])
        pattern = []
        for x, applies in zip(target.sites, present):
            position = mapping[x]
            if not applies:
                pattern.append(".")
            elif position < 0:
                pattern.append("-")
            else:
                base = read[position]
                pattern.append("M" if base == "C" else "U" if base == "T" else "N")
        record.update(votes_m=votes_m, votes_f=votes_f, allele=allele, pattern="".join(pattern))

        if record["coverage"] < p["min_coverage"]:
            record["status"] = "low_coverage"
        elif record["identity"] < p["min_identity"]:
            record["status"] = "low_identity"
        elif record["mismatch"] > p["max_mismatch"]:
            record["status"] = "high_mismatch"
        elif record["conversion"] < p["min_conversion"]:
            record["status"] = "low_conversion"
        elif p["max_unconverted"] is not None and unconverted > p["max_unconverted"]:
            record["status"] = "unconverted_cph"
        else:
            record["status"] = "pass"
        return record
