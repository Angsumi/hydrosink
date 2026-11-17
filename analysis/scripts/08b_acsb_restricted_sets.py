#!/usr/bin/env python3
"""acsB attribution, step 2 (POST-HOC sensitivity; the pre-planned read placement was not possible: reads deleted).

Builds per-sample marker tables in which acsB counts only reads whose best-hit acsB reference comes from a given
taxon set (NCBI taxonomy of the reference protein, results/acsb/acsb_ref_taxonomy.tsv):
  bacterial : all bacterial acsB references (archaeal ACS/CODH beta homologues excluded)
  bacillota : Bacillota references only (the phylum receiving 81-95% of acsB reads in every group)
acsB RPKM uses the mean length of the retained references; all other genes and qc_reads are unchanged.
Output: analysis/results/acsb_<set>/per_sample/*.markers.tsv, then 05a_assemble.py builds p1/p2 matrices there.
"""
import os, pathlib, statistics, subprocess, sys
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]; A = ROOT / "analysis"
tax = pd.read_csv(A / "results/acsb/acsb_ref_taxonomy.tsv", sep="\t").fillna("")
lab = pd.read_csv(A / "db/hydrosink.labels.tsv", sep="\t")
L = dict(zip(lab.seq_id, lab.len_aa)); rpoB_kb = statistics.mean(lab[lab.gene == "rpoB"].len_aa) * 3 / 1000
counts = pd.read_csv(A / "results/acsb/acsb_ref_counts.tsv", sep="\t")
SETS = {"bacterial": tax.domain.eq("Bacteria"), "bacillota": tax.phylum.eq("Bacillota")}
for name, mask in SETS.items():
    keep = set(tax.seq_id[mask]); kb = statistics.mean(L[s] for s in keep) * 3 / 1000
    n = counts[counts.seq_id.isin(keep)].groupby("run").reads.sum()
    out = A / f"results/acsb_{name}"; ps = out / "per_sample"; ps.mkdir(parents=True, exist_ok=True)
    k = 0
    for f in (A / "results/per_sample").glob("*.markers.tsv"):
        m = pd.read_csv(f, sep="\t")
        qc = m.qc_reads.iloc[0]; rpob = m[(m.level == "gene") & (m.feature == "rpoB")].reads.iloc[0]
        rr = rpob / rpoB_kb / (qc / 1e6)
        a = int(n.get(m.run.iloc[0], 0)); rp = a / kb / (qc / 1e6)
        i = m.index[(m.level == "gene") & (m.feature == "acsB")][0]
        m.loc[i, ["reads", "rpkm", "cpg"]] = [a, rp, rp / rr if rr else float("nan")]
        m.to_csv(ps / f.name, sep="\t", index=False, float_format="%.6g"); k += 1
    print(f"{name}: {len(keep)} references, {k} runs, acsB reads kept {int(n.sum())}")
    subprocess.run([sys.executable, "-I", str(A / "scripts/05a_assemble.py")], check=True,
                   env=dict(os.environ, HS_PS=str(ps), HS_OUTDIR=str(out)))
