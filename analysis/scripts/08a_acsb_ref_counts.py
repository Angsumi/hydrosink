#!/usr/bin/env python3
"""acsB attribution, step 1 (POST-HOC; replaces the pre-planned read placement, which needs the reads - deleted).

Counts, per run, the reads assigned to acsB (best hit = acsB reference, identity >= the acsB threshold, exactly as
in aggregate_markers.py), split by the reference sequence that was the best hit.
Output: analysis/results/acsb/acsb_ref_counts.tsv (run, seq_id, reads)
"""
import collections, gzip, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]; A = ROOT / "analysis"
OUT = A / "results/acsb"; OUT.mkdir(parents=True, exist_ok=True)
thr = {}
for line in open(A / "db/hydrosink.labels.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[1] == "acsB": thr[f[0]] = float(f[4])
alts = [A / f"sens/p1_rerun_hits/{d}/analysis/results/hits" for d in ("root", "root2")]
runs = [l.strip() for q in ("p1_queue.txt", "p2_queue.txt") for l in open(A / "inventory" / q) if l.strip()]
rows, missing = [], []
for run in runs:
    if not (A / "results/per_sample" / f"{run}.markers.tsv").exists(): continue
    h = A / "results/hits" / f"{run}.m8.gz"
    for alt in alts:
        if not h.exists(): h = alt / f"{run}.m8.gz"
    if not h.exists(): missing.append(run); continue
    c = collections.Counter()
    with gzip.open(h, "rt") as f:
        for line in f:
            q, s, pid = line.split("\t")[:3]
            if s in thr and float(pid) >= thr[s]: c[s] += 1
    rows += [(run, s, n) for s, n in c.items()]
with open(OUT / "acsb_ref_counts.tsv", "w") as o:
    o.write("run\tseq_id\treads\n")
    for r in rows: o.write("%s\t%s\t%d\n" % r)
print(f"{len(runs) - len(missing)} runs; missing hits: {missing}; {len(rows)} run x ref rows; "
      f"{sum(r[2] for r in rows)} acsB reads; {len({r[1] for r in rows})} / {len(thr)} references hit")
