#!/usr/bin/env python3
"""Benchmark step 1: pick and download one NCBI RefSeq assembly per panel species (genome FASTA + GFF3).
Rule: RefSeq reference genome if one exists; otherwise the RefSeq assembly with the highest contig N50.
Output: bench/genomes/<accession>/{genome.fna,genomic.gff}; results/bench/genome_accessions.tsv
Usage (VM): python bench_fetch.py <datasets binary>"""
import json, pathlib, shutil, subprocess, sys, zipfile

DS = sys.argv[1]
B = pathlib.Path(__file__).resolve().parent; G = B / "genomes"; G.mkdir(exist_ok=True)
R = B.parent / "results/bench"; R.mkdir(parents=True, exist_ok=True)
panel = [l.rstrip("\n").split("\t") for l in open(B / "genome_panel.tsv")][1:]


def summary(sp, ref):
    cmd = [DS, "summary", "genome", "taxon", sp, "--assembly-source", "RefSeq", "--as-json-lines"] + (["--reference"] if ref else [])
    p = subprocess.run(cmd, capture_output=True, text=True)
    return [json.loads(l) for l in p.stdout.splitlines() if l.strip()]


rows = []
for sp, grp in panel:
    rep = summary(sp, True); basis = "reference"
    if not rep:
        rep = summary(sp, False); basis = "max_contig_N50"
    rep = [r for r in rep if r.get("organism", {}).get("organism_name", "").startswith(sp)] or rep
    if not rep:
        rows.append((sp, grp, "", "", "NOT_FOUND")); print("NOT FOUND", sp); continue
    best = max(rep, key=lambda r: int(r.get("assembly_stats", {}).get("contig_n50", 0)))
    acc = best["accession"]; org = best["organism"]["organism_name"]
    d = G / acc
    if not (d / "genome.fna").exists():
        z = G / f"{acc}.zip"
        subprocess.run([DS, "download", "genome", "accession", acc, "--include", "genome,gff3", "--filename", str(z)], check=True,
                       capture_output=True)
        with zipfile.ZipFile(z) as zf:
            d.mkdir(exist_ok=True)
            for n in zf.namelist():
                if n.endswith(".fna"):
                    with zf.open(n) as s, open(d / "genome.fna", "ab") as t: shutil.copyfileobj(s, t)
                elif n.endswith(".gff"):
                    with zf.open(n) as s, open(d / "genomic.gff", "wb") as t: shutil.copyfileobj(s, t)
        z.unlink()
    rows.append((sp, grp, acc, org, basis)); print(sp, acc, org, basis, flush=True)
with open(R / "genome_accessions.tsv", "w") as o:
    o.write("species\tgroup\taccession\torganism_name\tselection\n")
    for r in rows: o.write("\t".join(r) + "\n")
