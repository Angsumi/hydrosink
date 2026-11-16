#!/usr/bin/env python3
"""Gate 1b decision: head-of-file vs random subsample, per gene, copies per genome (cpg).

Pre-specified rule (fixed 2026-11-16, before any random-subsample output existed):
  Genes evaluated: mcrA, acsB, rpoB-normalised respiratory sum (dsrA asrA nrfA narG napA frdA), NiFe, FeFe.
  For a gene with >= 100 reads in BOTH subsamples: PASS if |cpg_head/cpg_random - 1| <= 0.10.
  For a gene with < 100 reads in either: PASS if the read counts are consistent with Poisson
    sampling (two-sided exact conditional binomial test on the two counts, scaled by rpoB reads, p >= 0.05).
  Gate passes if every evaluated gene passes in both runs. On failure: switch to multi-offset sampling.
Output: analysis/results/gate1b_head_vs_random.tsv
"""
import pandas as pd, pathlib
from scipy.stats import binomtest

ROOT = pathlib.Path(__file__).resolve().parents[2]
PS, G = ROOT / "analysis/results/per_sample", ROOT / "analysis/results/gate1b"
RUNS = {"SRR28204439": "Ce14", "SRR28204452": "Rum14"}
RESP = ["dsrA", "asrA", "nrfA", "narG", "napA", "frdA"]
GENES = ["mcrA", "acsB", "resp", "NiFe", "FeFe"]


def genes(f):
    m = pd.read_csv(f, sep="\t"); m = m[m.level == "gene"].set_index("feature")
    get = lambda k, c: float(m.loc[k, c]) if k in m.index else 0.0
    d = {k: (get(k, "reads"), get(k, "cpg")) for k in GENES if k != "resp"}
    d["resp"] = (sum(get(k, "reads") for k in RESP), sum(get(k, "cpg") for k in RESP))
    return d, get("rpoB", "reads")


rows = []
for run, alias in RUNS.items():
    h, hr = genes(PS / f"{run}.markers.tsv"); r, rr = genes(G / f"{run}.random.markers.tsv")
    for k in GENES:
        (nh, ch), (nr, cr) = h[k], r[k]
        if min(nh, nr) >= 100:
            rel = ch / cr - 1 if cr else float("inf"); rule = "rel<=10%"; ok = abs(rel) <= 0.10; p = None
        else:
            rel = ch / cr - 1 if cr else float("nan"); rule = "poisson"
            frac = hr / (hr + rr)  # expected share of head reads under equal per-genome abundance
            p = binomtest(int(nh), int(nh + nr), frac).pvalue if nh + nr > 0 else 1.0; ok = p >= 0.05
        rows.append(dict(run=run, sample=alias, gene=k, reads_head=nh, reads_random=nr, cpg_head=ch,
                         cpg_random=cr, rel_diff=rel, rule=rule, poisson_p=p, pass_=ok))
d = pd.DataFrame(rows).rename(columns={"pass_": "pass"})
d.to_csv(ROOT / "analysis/results/gate1b_head_vs_random.tsv", sep="\t", index=False)
print(d.to_string(index=False))
print("GATE 1b:", "PASS" if d["pass"].all() else "FAIL")
