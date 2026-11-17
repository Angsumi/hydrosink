#!/usr/bin/env python3
"""Pre-specified sensitivity (config.yaml diamond.sensitivity_identity_shift: [-10, +10]).

Re-aggregates the saved DIAMOND best hits (results/hits/<run>.m8.gz) with every per-gene identity threshold shifted
by SHIFT percentage points (floored at 40, the DIAMOND --id used for the search, below which no hits were kept).
Counting, gene lengths, RPKM and copies-per-genome follow aggregate_markers.py exactly; only the threshold changes.
Output: analysis/results/sens_id_<m10|p10>/per_sample/<run>.markers.tsv (same format), then 05a_assemble.py builds
the matrices in that directory. Hits for runs re-run in analysis/sens/p1_rerun_hits are taken from there.
Usage: 05c_identity_sensitivity.py SHIFT
"""
import collections, gzip, json, os, pathlib, statistics, subprocess, sys

SHIFT = int(sys.argv[1]); FLOOR = 40.0
ROOT = pathlib.Path(__file__).resolve().parents[2]; A = ROOT / "analysis"
OUT = A / "results" / f"sens_id_{'m' if SHIFT < 0 else 'p'}{abs(SHIFT)}"; PSO = OUT / "per_sample"; PSO.mkdir(parents=True, exist_ok=True)
lab, rg, rs = {}, collections.defaultdict(list), collections.defaultdict(list)
with open(A / "db/hydrosink.labels.tsv") as f:
    next(f)
    for line in f:
        sid, gene, guild, sub, mid, L, src = line.rstrip("\n").split("\t")
        lab[sid] = (gene, sub, max(FLOOR, float(mid) + SHIFT)); rg[gene].append(int(L))
        if sub: rs[f"{gene}_{sub}"].append(int(L))
kb = lambda l: statistics.mean(l) * 3 / 1000
runs = [l.strip() for q in ("p1_queue.txt", "p2_queue.txt") for l in open(A / "inventory" / q) if l.strip()]
alt = A / "sens/p1_rerun_hits/root/analysis/results/hits"
missing = []
for run in runs:
    ps = A / "results/per_sample" / f"{run}.markers.tsv"
    if not ps.exists(): continue
    h = A / "results/hits" / f"{run}.m8.gz"
    if not h.exists(): h = alt / f"{run}.m8.gz"
    if not h.exists(): missing.append(run); continue
    first = open(ps).readlines()[1].split("\t"); qc, ns = int(first[6]), int(first[7])   # same qc_reads as original
    g, s = collections.Counter(), collections.Counter()
    with gzip.open(h, "rt") as f:
        for line in f:
            q, sid, pid = line.split("\t")[:3]; gene, sub, mid = lab[sid]
            if float(pid) < mid: continue
            g[gene] += 1
            if sub: s[f"{gene}_{sub}"] += 1
    rp = lambda n, L: n / kb(L) / (qc / 1e6) if qc else float("nan")
    rb = rp(g["rpoB"], rg["rpoB"])
    with open(PSO / f"{run}.markers.tsv", "w") as o:
        o.write("run\tfeature\tlevel\treads\trpkm\tcpg\tqc_reads\tn_streamed\n")
        for k in sorted(set(g) | {"mcrA", "acsB", "rpoB"}):
            r = rp(g[k], rg[k]); o.write(f"{run}\t{k}\tgene\t{g[k]}\t{r:.6g}\t{(r / rb if rb else float('nan')):.6g}\t{qc}\t{ns}\n")
        for k in sorted(s):
            r = rp(s[k], rs[k]); o.write(f"{run}\t{k}\tsubgroup\t{s[k]}\t{r:.6g}\t{(r / rb if rb else float('nan')):.6g}\t{qc}\t{ns}\n")
print(f"shift {SHIFT:+d}: {len(list(PSO.glob('*.markers.tsv')))} runs re-aggregated; no hits: {len(missing)} {missing[:20]}")
env = dict(os.environ, HS_PS=str(PSO), HS_OUTDIR=str(OUT))
subprocess.run([sys.executable, "-I", str(A / "scripts/05a_assemble.py")], env=env, check=True)
