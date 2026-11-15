#!/usr/bin/env python3
"""Convert DIAMOND hits into per-guild read counts normalised to genome equivalents (GE).

GE = median coverage across SingleM single-copy marker genes (sum of OTU coverage per marker).
Marker abundance (copies per GE) = reads_hit * read_len / (mean_target_len_nt) / GE_coverage_in_reads,
approximated here as reads per kb of target per GE, the convention used for cross-sample comparison.
Guild-specific identity thresholds come from the labels file (column min_id).
"""
import argparse
import pandas as pd


def genome_equivalents(otu_path: str) -> float:
    otu = pd.read_csv(otu_path, sep="\t")
    cov_col = "coverage" if "coverage" in otu.columns else otu.columns[-2]
    per_gene = otu.groupby("gene")[cov_col].sum()
    return float(per_gene.median())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m8", required=True); ap.add_argument("--labels", required=True)
    ap.add_argument("--singlem", required=True); ap.add_argument("--run", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cols = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"]
    hits = pd.read_csv(a.m8, sep="\t", names=cols)
    lab = pd.read_csv(a.labels, sep="\t")  # columns: sseqid, marker, guild, min_id, confidence
    hits = hits.merge(lab, on="sseqid", how="inner")
    hits = hits[hits["pident"] >= hits["min_id"]]
    ge = genome_equivalents(a.singlem)
    g = hits.groupby(["marker", "guild"]).agg(reads=("qseqid", "nunique"), mean_slen_aa=("slen", "mean")).reset_index()
    g["rpk"] = g["reads"] / (g["mean_slen_aa"] * 3 / 1000)
    g["copies_per_GE"] = g["rpk"] / ge if ge > 0 else float("nan")
    g.insert(0, "run_accession", a.run)
    g["genome_equivalent_cov"] = ge
    g.to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
