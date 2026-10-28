"""Simulated Nanopore amplicon bisulfite reads with known truth, for validation."""
from __future__ import annotations

import json
from pathlib import Path
import random

from .reference import IUPAC, revcomp

BASES = "ACGT"


def bisulfite(sequence, methylated_sites, conversion, rng):
    bases = list(sequence)
    for x, base in enumerate(bases):
        if base == "C" and x not in methylated_sites and rng.random() < conversion:
            bases[x] = "T"
    return "".join(bases)


def nanopore_errors(sequence, rate, rng):
    """Substitutions, insertions and deletions in a 2:1:1 ratio."""
    out = []
    for base in sequence:
        roll = rng.random()
        if roll < rate / 2:
            out.append(rng.choice([b for b in BASES if b != base]))
        elif roll < rate * 3 / 4:
            out.append(base + rng.choice(BASES))
        elif roll < rate:
            continue
        else:
            out.append(base)
    return "".join(out)


def simulate(path, truth_path, amplicon, alleles, samples, reads_per_sample, error_rate=0.05,
             conversion=0.99, seed=1, unconverted_fraction=0.0, junk_fraction=0.0):
    """Write FASTQ and a truth TSV.

    amplicon: dict(sequence, forward, reverse) on the bisulfite strand X.
    alleles: {"mother": X sequence, "father": X sequence}.
    samples: list of dict(name, tag_f, tag_r, methylation={"mother": p, "father": p}).
    """
    rng = random.Random(seed)
    sequence, forward, reverse = amplicon["sequence"], amplicon["forward"], amplicon["reverse"]
    lf, lr = len(forward), len(reverse)
    records, truth = [], []
    number = 0
    for sample in samples:
        for _ in range(reads_per_sample):
            number += 1
            allele = rng.choice(["mother", "father"])
            template = (sample.get("alleles") or alleles)[allele]
            sites = [x for x in range(len(template) - 1) if template[x] == "C" and template[x + 1] == "G"]
            methylated = {x for x in sites if rng.random() < sample["methylation"][allele]}
            rate = 0.2 if rng.random() < unconverted_fraction else conversion
            converted = bisulfite(template, methylated, rate, rng)
            primer_f = "".join(rng.choice(IUPAC[b]) for b in forward)
            primer_r = "".join(rng.choice(IUPAC[b]) for b in reverse)
            product = primer_f + converted[lf:len(converted) - lr] + revcomp(primer_r)
            molecule = sample["tag_f"] + product + revcomp(sample["tag_r"])
            if rng.random() < junk_fraction:
                molecule = "".join(rng.choice(BASES) for _ in range(len(molecule)))
                allele = "junk"
            flank = lambda: "".join(rng.choice(BASES) for _ in range(rng.randint(10, 60)))  # noqa: E731
            read = flank() + molecule + flank()
            if rng.random() < 0.5:
                read = revcomp(read)
            read = nanopore_errors(read, error_rate, rng)
            name = f"sim{number}"
            records.append(f"@{name}\n{read}\n+\n{'5' * len(read)}\n")
            truth.append(dict(read=name, sample=sample["name"], allele=allele,
                              methylated=sorted(methylated), sites=sites))
    Path(path).write_text("".join(records), encoding="ascii")
    Path(truth_path).write_text("\n".join(json.dumps(t) for t in truth) + "\n", encoding="utf-8")
