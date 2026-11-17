#!/usr/bin/env python3
"""acsB attribution, descriptive (POST-HOC): share of acsB reads by NCBI phylum / family of the best-hit reference,
per analysis group (P1 chamber; P2 fermentation site, modelled set). Outputs in analysis/results/acsb/."""
import pathlib
import pandas as pd
A = pathlib.Path(__file__).resolve().parents[2] / "analysis"; R = A / "results"
c = pd.read_csv(R / "acsb/acsb_ref_counts.tsv", sep="\t").merge(pd.read_csv(R / "acsb/acsb_ref_taxonomy.tsv", sep="\t").fillna(""), on="seq_id")
p1 = pd.read_csv(R / "p1_matrix.tsv", sep="\t"); p1 = p1[p1.qc_pass]; p2 = pd.read_csv(R / "p2_analysis_set.tsv", sep="\t")
grp = {**{r: "P1 " + ch for r, ch in zip(p1.run_accession, p1.chamber)}, **{r: "P2 " + s for r, s in zip(p2.run_accession, p2.fermentation_site)}}
c["group"] = c.run.map(grp); c = c[c.group.notna()]
for lvl, fn in (("phylum", "acsb_reads_by_phylum.tsv"), ("family", "acsb_reads_by_family.tsv")):
    key = c.domain + " | " + c.phylum + (" | " + c.family if lvl == "family" else "")
    t = c.assign(taxon=key).pivot_table(index="taxon", columns="group", values="reads", aggfunc="sum", fill_value=0)
    pct = (t / t.sum() * 100).round(2); pct.loc["TOTAL_READS"] = t.sum()
    pct.sort_values("P2 hindgut", ascending=False).to_csv(R / "acsb" / fn, sep="\t")
print(pd.read_csv(R / "acsb/acsb_reads_by_phylum.tsv", sep="\t").head(6).to_string(index=False))
