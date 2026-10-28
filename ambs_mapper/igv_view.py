"""Embedded IGV (igv.js): viewer BAMs with CpG calls as base modifications, and the browser configuration.

igv.js has no bisulfite colouring, but it colours base modifications from MM/ML
tags ("basemod2": modified red, unmodified blue). For the viewer, each read's
CpG calls are written as 5mC probabilities (255 methylated, 0 unmethylated);
an unmethylated CpG reads T (or A on the bottom strand), which MM cannot point
at, so that base is shown as C (G) in the viewer copy only. Soft clips keep
the tag/adapter sequence outside the alignment visible. The downloadable BAMs
are left unchanged.
"""
from __future__ import annotations

from bisect import bisect_left
import json
from pathlib import Path
import re
import subprocess

from .reference import reference_fasta, revcomp

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def methylation_tags(fields, positions, top):
    """Add MM/ML to one SAM record (list of fields); positions are 0-based reference
    positions of the measured base of each CpG (top C, or top G for the bottom strand)."""
    start, seq = int(fields[3]) - 1, list(fields[9])
    index, r, q = {}, start, 0
    for n, op in CIGAR.findall(fields[5]):
        n = int(n)
        if op in "M=X":
            for k in range(n):
                index[r + k] = q + k
            r, q = r + n, q + n
        elif op in "DN":
            r += n
        elif op in "IS":
            q += n
    methylated_base, converted_base = ("C", "T") if top else ("G", "A")
    calls = []
    for p in positions:
        i = index.get(p)
        if i is None or seq[i] not in (methylated_base, converted_base):
            continue
        calls.append((i, seq[i] == methylated_base))
        seq[i] = methylated_base            # viewer copy only: MM can only point at a C
    if not calls:
        return fields
    reverse = int(fields[1]) & 16
    length = len(seq)
    original = revcomp("".join(seq)) if reverse else "".join(seq)
    cs = [k for k, base in enumerate(original) if base == "C"]
    ranked = sorted((bisect_left(cs, length - 1 - i if reverse else i), m) for i, m in calls)
    skips, previous = [], -1
    for rank, _ in ranked:
        skips.append(rank - previous - 1)
        previous = rank
    fields = fields[:9] + ["".join(seq)] + fields[10:]
    fields = [f for f in fields if not f.startswith(("MM:", "ML:", "Mm:", "Ml:"))]
    fields.append("MM:Z:C+m?," + ",".join(map(str, skips)) + ";")
    fields.append("ML:B:C," + ",".join("255" if m else "0" for _, m in ranked))
    return fields


def view_bam(source, target, sites):
    """Copy a BAM for the viewer, adding MM/ML from the read's amplicon (YA tag)."""
    reader = subprocess.Popen(["samtools", "view", "-h", str(source)], stdout=subprocess.PIPE, text=True)
    writer = subprocess.Popen(["samtools", "view", "-b", "-o", str(target), "-"], stdin=subprocess.PIPE, text=True)
    for line in reader.stdout:
        if not line.startswith("@"):
            fields = line.rstrip("\n").split("\t")
            amplicon = next((f[5:] for f in fields[11:] if f.startswith("YA:Z:")), None)
            if amplicon in sites:
                fields = methylation_tags(fields, *sites[amplicon])
            line = "\t".join(fields) + "\n"
        writer.stdin.write(line)
    writer.stdin.close()
    if reader.wait() or writer.wait():
        raise RuntimeError(f"IGV表示用のBAMを作成できませんでした: {source.name}")
    subprocess.run(["samtools", "index", str(target)], check=True, capture_output=True)


def build(out, details, entries):
    """Viewer BAMs for the per-sample and tag-unknown BAMs (entries from pipeline.split_bam),
    CpG/amplicon tracks, and the track list. Returns the number of tracks."""
    out = Path(out)
    folder = out / "igv"
    folder.mkdir(exist_ok=True)
    sites = {}
    beds = {"amplicons.bed": [], "cpg.bed": []}
    for d in details:
        top = d["strand"] == "+"
        # measured base: the CpG's top-strand C, or its top-strand G (= bottom-strand C)
        sites[d["name"]] = ([p - 1 if top else p for p in d["site_positions"]], top)
        beds["amplicons.bed"].append(f"{d['chrom']}\t{d['start'] - 1}\t{d['end']}\t{d['name']}\t0\t{d['strand']}")
        for n, p in enumerate(d["site_positions"], 1):
            beds["cpg.bed"].append(f"{d['chrom']}\t{p - 1}\t{p + 1}\tCpG{n}")
    for name, lines in beds.items():
        (folder / name).write_text("\n".join(lines) + "\n", encoding="utf-8")
    tracks = []
    for n, entry in enumerate(entries):
        asset = f"track{n}.bam"
        view_bam(entry["path"], folder / asset, sites)
        tracks.append(dict(asset=asset, name=entry["name"], amplicon=entry["amplicon"], unknown=entry["unknown"]))
    amplicons = [dict(name=d["name"], locus=f"{d['chrom']}:{max(1, d['start'] - 60)}-{d['end'] + 60}") for d in details]
    (folder / "tracks.json").write_text(json.dumps(dict(tracks=tracks, amplicons=amplicons), ensure_ascii=False),
                                        encoding="utf-8")
    return len(tracks)


def assets(out, resources):
    """Files the viewer may fetch: the reference, the BED tracks and the viewer BAMs."""
    folder = Path(out) / "igv"
    fasta = reference_fasta(resources)
    files = {"reference.fa": fasta, "reference.fa.fai": Path(str(fasta) + ".fai"),
             "amplicons.bed": folder / "amplicons.bed", "cpg.bed": folder / "cpg.bed"}
    for track in json.loads((folder / "tracks.json").read_text(encoding="utf-8"))["tracks"]:
        files[track["asset"]] = folder / track["asset"]
        files[track["asset"] + ".bai"] = folder / (track["asset"] + ".bai")
    return files


def configuration(out, job_id):
    listing = json.loads((Path(out) / "igv" / "tracks.json").read_text(encoding="utf-8"))
    base = f"/api/jobs/{job_id}/igv/"
    first = listing["amplicons"][0]["locus"] if listing["amplicons"] else "1"
    tracks = [dict(name="アンプリコン", type="annotation", format="bed", url=base + "amplicons.bed", indexed=False,
                   displayMode="EXPANDED", order=0),
              dict(name="CpG番号（写真の列と同じ）", type="annotation", format="bed", url=base + "cpg.bed",
                   indexed=False, displayMode="EXPANDED", color="#c0392b", order=1)]
    alignments = [dict(name=t["name"], type="alignment", format="bam", url=base + t["asset"],
                       indexURL=base + t["asset"] + ".bai", showSoftClips=True, showMismatches=False,
                       colorBy="basemod2", displayMode="SQUISHED", showCoverage=True, height=220,
                       amplicon=t["amplicon"], unknown=t["unknown"]) for t in listing["tracks"]]
    return dict(reference=dict(id="GRCm39", name="GRCm39", fastaURL=base + "reference.fa",
                               indexURL=base + "reference.fa.fai"),
                locus=first, tracks=tracks, showUCSCLink=False,
                loadDefaultGenomes=False), alignments, listing["amplicons"]   # no fetch from igv.org


def byte_range(header, size):
    """Inclusive single byte range of a Range header (as in KO_Checker)."""
    match = re.fullmatch(r"bytes=(\d*)-(\d*)", header or "")
    if not match or not any(match.groups()) or size == 0:
        raise ValueError("Invalid byte range")
    first, last = match.groups()
    if first:
        start, end = int(first), min(int(last), size - 1) if last else size - 1
    else:
        if int(last) == 0:
            raise ValueError("Invalid suffix range")
        start, end = max(0, size - int(last)), size - 1
    if start >= size or start > end:
        raise ValueError("Unsatisfiable byte range")
    return start, end
