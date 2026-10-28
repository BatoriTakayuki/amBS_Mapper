"""Illumina Nextera dual-index layout: [P5][i5][A1][primer] ... [primer][A2][i7][P7].

Index sequences are written as in a MiSeq sample sheet (i5 forward-strand
workflow, i7 as read). The 14 i7 and 12 i5 indexes seen in the first shared
Nanopore run matched these tables exactly.
"""
from __future__ import annotations

import re

from .reference import revcomp

A1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"   # Nextera read 1 adapter (5' tail of one primer)
A2 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"  # Nextera read 2 adapter (5' tail of the other)
P5_END = "CGAGATCTACAC"                        # 3' end of P5, directly 5' of i5
P7_END = "GGCATACGAGAT"                        # 3' end of P7, directly 5' of i7 on the oligo

I7 = {
    "N701": "TAAGGCGA", "N702": "CGTACTAG", "N703": "AGGCAGAA", "N704": "TCCTGAGC", "N705": "GGACTCCT",
    "N706": "TAGGCATG", "N707": "CTCTCTAC", "N708": "CAGAGAGG", "N709": "GCTACGCT", "N710": "CGAGGCTG",
    "N711": "AAGAGGCA", "N712": "GTAGAGGA", "N714": "GCTCATGA", "N715": "ATCTCAGG", "N716": "ACTCGCTA",
    "N718": "GGAGCTAC", "N719": "GCGTAGTA", "N720": "CGGAGCCT", "N721": "TACGCTGC", "N722": "ATGCGCAG",
    "N723": "TAGCGCTC", "N724": "ACTGAGCG", "N726": "CCTAAGAC", "N727": "CGATCAGT", "N728": "TGCAGCTA",
    "N729": "TCGACGTC",
}
I5 = {
    "S502": "CTCTCTAT", "S503": "TATCCTCT", "S505": "GTAAGGAG", "S506": "ACTGCATA", "S507": "AAGGAGTA",
    "S508": "CTAAGCCT", "S510": "CGTCTAAT", "S511": "TCTCTCCG", "S513": "TCGACTAG", "S515": "TTCTAGCT",
    "S516": "CCTAGAGT", "S517": "GCGTAAGA", "S518": "CTATTAAG", "S520": "AAGGCTAT", "S521": "GAGCCTTA",
    "S522": "TTATGCGA",
}
LAYOUTS = {"inline": "インライン（プライマー直前のタグ）", "nextera": "Nextera デュアルインデックス（i5/i7）"}
_SEQUENCE = re.compile(r"[ACGT]+")


def parse_index(value, kind):
    """Sample-sheet index name or sequence -> sequence as it appears on the oligo, 5' of the adapter.

    kind is "i5" or "i7"; i7 is stored reverse-complemented because the oligo
    carries the complement of the index read.
    """
    text = re.sub(r"\s+", "", str(value or "")).upper()
    if not text:
        return ""
    table = I5 if kind == "i5" else I7
    # Accept "N701", "701", "S520", "520", "[N7]01"-style names.
    name = text if text in table else next((key for key in table if key[1:] == text), None)
    if name is None:
        # Sample sheets differ in index orientation (e.g. i5 on forward- vs reverse-
        # complement workflows); a kit index written either way is recognised.
        name = next((key for key, value in table.items() if text in (value, revcomp(value))), None)
    sequence = table[name] if name else text
    if not _SEQUENCE.fullmatch(sequence):
        raise ValueError(f"{kind}インデックス「{value}」はNextera XTの名前（例 {'S520' if kind == 'i5' else 'N712'}）"
                         "かACGTの配列で入力してください。")
    return sequence if kind == "i5" else revcomp(sequence)


def kind_of(value):
    """"i5" or "i7" if the entry names or spells a kit index of that kind, else None."""
    text = re.sub(r"\s+", "", str(value or "")).upper()
    for kind, table in (("i5", I5), ("i7", I7)):
        if text in table or any(text in (seq, revcomp(seq)) for seq in table.values()):
            return kind
    return None


def label(oligo_sequence, kind):
    """Index name (or sample-sheet sequence) for an oligo-orientation sequence."""
    if not oligo_sequence:
        return ""
    sequence = oligo_sequence if kind == "i5" else revcomp(oligo_sequence)
    table = I5 if kind == "i5" else I7
    return next((name for name, value in table.items() if value == sequence), sequence)


def whitelist(kind):
    return sorted(I5.values()) if kind == "i5" else sorted(revcomp(value) for value in I7.values())
