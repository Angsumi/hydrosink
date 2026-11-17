#!/usr/bin/env python3
"""Assemble per-sample marker tables + metadata into analysis-ready matrices.

Outputs:
  analysis/results/p1_matrix.tsv  one row per P1 run (study, host_species, animal, compartment, chamber, markers)
  analysis/results/p2_matrix.tsv  one row per P2 run (host traits, study, markers)
Marker columns: <gene>_cpg (copies per genome) and <gene>_n (reads) for the H2 markers below,
plus derived sums: resp_cpg (dsrA+asrA+nrfA+narG+napA+frdA), NiFe1_cpg (group-1 uptake [NiFe]),
FeFeB_cpg, rpoB_n. Pseudocounts and log-ratios are computed in R (05_meta_analysis.R) so that
c = half-min-nonzero is defined on the analysed set.
QC flag: rpoB_n >= 200 (enough single-copy reads for a stable normaliser); others kept but flagged.
"""
import pandas as pd, numpy as np, pathlib, os

ROOT = pathlib.Path(__file__).resolve().parents[2]
PS = pathlib.Path(os.environ.get("HS_PS", ROOT / "analysis/results/per_sample"))   # override: sensitivity re-aggregations only
OUTD = pathlib.Path(os.environ.get("HS_OUTDIR", ROOT / "analysis/results"))
GENES = ["mcrA", "acsB", "dsrA", "asrA", "nrfA", "narG", "napA", "frdA", "sdhA", "cooS", "fdhA", "NiFe", "FeFe", "rpoB"]
RESP = ["dsrA", "asrA", "nrfA", "narG", "napA", "frdA"]
MIN_RPOB = 200


def load(run):
    f = PS / f"{run}.markers.tsv"
    if not f.exists():
        return None
    m = pd.read_csv(f, sep="\t")
    g = m[m.level == "gene"].set_index("feature"); s = m[m.level == "subgroup"].set_index("feature")
    row = {"run_accession": run, "qc_reads": int(m.qc_reads.iloc[0]), "n_streamed": int(m.n_streamed.iloc[0])}
    for k in GENES:
        row[f"{k}_cpg"] = float(g.loc[k, "cpg"]) if k in g.index else 0.0
        row[f"{k}_n"] = int(g.loc[k, "reads"]) if k in g.index else 0
    row["resp_cpg"] = sum(row[f"{k}_cpg"] for k in RESP)
    row["resp_n"] = sum(row[f"{k}_n"] for k in RESP)
    row["NiFe1_cpg"] = float(s[s.index.str.match(r"NiFe_1")].cpg.sum()) if len(s) else 0.0
    row["FeFeB_cpg"] = float(s.loc["FeFe_B", "cpg"]) if "FeFe_B" in s.index else 0.0
    return row


def build(inv_path, out_path, keep_cols):
    inv = pd.read_csv(inv_path, sep="\t", dtype=str)
    rows = [r for r in (load(x) for x in inv.run_accession) if r]
    if not rows:
        print(f"{out_path.name}: no quantified runs yet"); return
    M = pd.DataFrame(rows)
    d = inv[keep_cols].merge(M, on="run_accession", how="inner")
    d["qc_pass"] = d.rpoB_n >= MIN_RPOB
    cur = ROOT / "analysis/inventory/p2_curation.tsv"
    if "captivity" in d.columns and cur.exists():   # [prereg-v1.1] outcome-blind curated metadata (01d_curate_metadata.py)
        c = pd.read_csv(cur, sep="\t", dtype=str)[["run_accession", "captivity_v11", "captivity_basis", "species_v11", "species_basis"]]
        ott = set(pd.read_csv(ROOT / "analysis/inventory/host_ott_map.tsv", sep="\t", dtype=str).species)
        d = d.merge(c, on="run_accession", how="left")
        d["captivity_v1"] = d.captivity; d["captivity"] = d.captivity_v11.fillna(d.captivity)
        swap = d.species_v11.notna() & (d.species_v11 != d.species) & d.species_v11.isin(ott)
        d["species_v1"] = d.species; d.loc[swap, "species"] = d.loc[swap, "species_v11"]
        d["genus_only"] = d.species_basis.eq("genus_only_unresolved") | (~swap & d.species_basis.str.startswith("manual", na=False))
        d = d.drop(columns=["captivity_v11", "species_v11"])
    d.to_csv(out_path, sep="\t", index=False)
    print(f"{out_path.name}: {len(d)} quantified / {len(inv)} planned; qc_pass {d.qc_pass.sum()}")


build(ROOT / "analysis/inventory/p1_runs.tsv", OUTD / "p1_matrix.tsv",
      ["run_accession", "study", "host_species", "animal", "compartment", "chamber", "sample_alias"])
build(ROOT / "analysis/inventory/p2_runs.tsv", OUTD / "p2_matrix.tsv",
      ["run_accession", "study_accession", "sample_accession", "species", "host_tax_id", "vclass",
       "fermentation_site", "diet_class", "plant_pct", "captivity", "sample_type", "read_len",
       "instrument_model", "country"])
