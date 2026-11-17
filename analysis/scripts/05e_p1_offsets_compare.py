#!/usr/bin/env python3
"""Pre-specified descriptive sensitivity (amendment v1.2 §2): 5 seeded P1 runs, head (first 5 M R1 reads, used in P1)
vs multi-offset sampling (P2 protocol). Per gene group: copies per genome under both, ratio offsets/head, and an
approximate 95% CI of the ratio from Poisson counts (log-ratio SE = sqrt(1/a + 1/b + 1/ra + 1/rb), a/b marker reads,
ra/rb rpoB reads). Reported whatever it shows; does not alter the P1 analysis.
Output: analysis/results/p1_offsets_vs_head.tsv
"""
import math, pathlib
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]; A = ROOT / "analysis"
OFF = A / "sens/p1_offsets/root/analysis/results/per_sample"
RESP = ["dsrA", "asrA", "nrfA", "narG", "napA", "frdA"]
GROUPS = {"mcrA": ["mcrA"], "acsB": ["acsB"], "respiratory_sum": RESP, "NiFe": ["NiFe"], "FeFe": ["FeFe"]}


def load(p):
    m = pd.read_csv(p, sep="\t"); g = m[m.level == "gene"].set_index("feature")
    return g


rows = []
for run in [l.strip() for l in open(A / "sens/p1_offsets/runs.txt") if l.strip()]:
    h, o = load(A / "results/per_sample" / f"{run}.markers.tsv"), load(OFF / f"{run}.markers.tsv")
    for grp, genes in GROUPS.items():
        get = lambda t, c: sum(float(t.loc[x, c]) if x in t.index else 0.0 for x in genes)
        ch, co, nh, no = get(h, "cpg"), get(o, "cpg"), get(h, "reads"), get(o, "reads")
        rh, ro = float(h.loc["rpoB", "reads"]), float(o.loc["rpoB", "reads"])
        r = co / ch if ch > 0 else float("nan")
        se = math.sqrt(1 / nh + 1 / no + 1 / rh + 1 / ro) if nh > 0 and no > 0 else float("nan")
        rows.append(dict(run=run, group=grp, reads_head=int(nh), reads_offsets=int(no), cpg_head=ch, cpg_offsets=co,
                         ratio_offsets_head=r, ci_lb=r * math.exp(-1.96 * se) if se == se else float("nan"),
                         ci_ub=r * math.exp(1.96 * se) if se == se else float("nan")))
d = pd.DataFrame(rows)
d.to_csv(A / "results/p1_offsets_vs_head.tsv", sep="\t", index=False)
print(d.round(3).to_string(index=False))
