"""Strain-specific SNPs from the Mouse Genomes Project (REL-2112 v8, GRCm39).

Variants are read remotely through the tabix index, so only the amplicon
regions are transferred. C57BL/6J is the reference genome itself.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess

from .reference import revcomp, x_to_genome

SNP_URL = os.environ.get("AMBS_SNP_VCF", "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/"
                         "REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz")
INDEL_URL = os.environ.get("AMBS_INDEL_VCF", "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/"
                           "REL-2112-v8-SNPs_Indels/mgp_REL2021_indels.vcf.gz")
REFERENCE_STRAIN = "C57BL_6J"
ALIASES = {"B6": REFERENCE_STRAIN, "B6J": REFERENCE_STRAIN, "C57BL/6": REFERENCE_STRAIN,
           "C57BL/6J": REFERENCE_STRAIN, "BL6": REFERENCE_STRAIN, "JF1": "JF1_MsJ",
           "B6N": "C57BL_6NJ", "MSM": "MSM_Ms"}


def _bcftools(arguments, cache):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    # bcftools stores the downloaded remote index in the working directory.
    result = subprocess.run(["bcftools", *arguments], cwd=cache, capture_output=True, text=True, timeout=600)
    if result.returncode:
        raise RuntimeError("MGPのVCFを取得できませんでした。インターネット接続を確認してください: "
                           + result.stderr.strip()[-400:])
    return result.stdout


def strain_list(resources):
    cache = Path(resources) / "strains.json"
    if cache.is_file():
        return json.loads(cache.read_text(encoding="utf-8"))
    names = [REFERENCE_STRAIN] + _bcftools(["query", "-l", SNP_URL], Path(resources) / "vcf_cache").split()
    cache.write_text(json.dumps(names), encoding="utf-8")
    return names


def resolve_strain(name, available):
    """Accept full MGP names, common aliases and unique prefixes (JF1 -> JF1_MsJ)."""
    text = str(name or "").strip()
    if not text:
        raise ValueError("系統名が空です。")
    key = text.upper().replace(" ", "")
    for alias, target in ALIASES.items():
        if key == alias.upper() and target in available:
            return target
    normalized = {strain.upper().replace("/", "_"): strain for strain in available}
    if key.replace("/", "_") in normalized:
        return normalized[key.replace("/", "_")]
    matches = [strain for upper, strain in normalized.items() if upper.startswith(key.replace("/", "_"))]
    if len(matches) == 1:
        return matches[0]
    hint = "、".join(matches[:6]) if matches else "一覧から選んでください"
    raise ValueError(f"系統「{text}」を特定できません（{hint}）。")


def _genotype_base(gt, fi, ref, alts):
    """Homozygous high-confidence base, otherwise None."""
    if fi not in ("1", "."):
        return None
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or alleles[0] != alleles[1] or alleles[0] == ".":
        return None
    index = int(alleles[0])
    return ref if index == 0 else alts[index - 1]


def query_region(resources, chrom, start, end, strains, url=SNP_URL):
    """Variants in [start, end) (0-based) as dicts with per-strain alleles on the top strand."""
    strains = sorted(set(strains))
    others = [strain for strain in strains if strain != REFERENCE_STRAIN]
    key = hashlib.sha1(f"{url}|{chrom}|{start}|{end}|{','.join(strains)}".encode()).hexdigest()[:20]
    cache = Path(resources) / "variant_cache" / f"{key}.json"
    if cache.is_file():
        return json.loads(cache.read_text(encoding="utf-8"))
    variants = []
    if others:
        output = _bcftools(["query", "-r", f"{chrom}:{start + 1}-{end}", "-s", ",".join(others),
                            "-f", "%POS\t%REF\t%ALT\t%FILTER[\t%GT:%FI]\n", url], Path(resources) / "vcf_cache")
        for line in output.splitlines():
            pos, ref, alt, flt, *calls = line.split("\t")
            if flt not in ("PASS", "."):
                continue
            alts = alt.split(",")
            alleles = {REFERENCE_STRAIN: ref}
            for strain, call in zip(others, calls):
                gt, _, fi = call.partition(":")
                alleles[strain] = _genotype_base(gt, fi or ".", ref, alts)
            variants.append(dict(pos=int(pos) - 1, ref=ref, alleles=alleles))
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(variants), encoding="utf-8")
    return variants


def strand_snps(product, variants, mother, father):
    """SNPs differing between two strains, expressed on the bisulfite strand X.

    Returns (snps, skipped) where each SNP is (x, mother_base, father_base,
    informative). Under bisulfite conversion only C on X can change (to T), so a
    C/T difference on X cannot separate alleles.
    """
    snps, skipped = [], 0
    for variant in variants:
        a, b = variant["alleles"].get(mother), variant["alleles"].get(father)
        if len(variant["ref"]) != 1 or a is None or b is None:
            if a != b:
                skipped += 1
            continue
        if a == b or len(a) != 1 or len(b) != 1:
            continue
        pos = variant["pos"]
        if not product["start"] <= pos < product["end"]:
            continue
        if product["strand"] == "+":
            x = pos - product["start"]
        else:
            x = product["end"] - 1 - pos
            a, b = revcomp(a), revcomp(b)
        snps.append(dict(x=x, genome=pos + 1, mother=a, father=b, informative={a, b} != {"C", "T"}))
    snps.sort(key=lambda snp: snp["x"])
    return snps, skipped


def strand_indels(resources, product, mother, father):
    """Indel positions that differ between the strains (reported as a warning)."""
    variants = query_region(resources, product["chrom"], product["start"], product["end"],
                            [mother, father], INDEL_URL)
    result = []
    for variant in variants:
        a, b = variant["alleles"].get(mother), variant["alleles"].get(father)
        if a != b and product["start"] <= variant["pos"] < product["end"]:
            result.append(dict(genome=variant["pos"] + 1, mother=a, father=b))
    return result


__all__ = ["strain_list", "resolve_strain", "query_region", "strand_snps", "strand_indels",
           "REFERENCE_STRAIN", "x_to_genome"]
