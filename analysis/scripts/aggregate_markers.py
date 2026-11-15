#!/usr/bin/env python3
"""Aggregate DIAMOND best hits to per-sample marker abundances.

A read is counted for gene g if its best hit is a reference of g and pident >= min_id(g)
(Greening-lab per-gene thresholds; query cover >= 80 % enforced in DIAMOND).
RPKM_g = reads_g / (mean_ref_len_g[aa] * 3 / 1000) / (qc_reads / 1e6)
copies_per_genome_g = RPKM_g / RPKM_rpoB
Output (long format): run, feature, level, reads, rpkm, cpg, qc_reads, n_streamed
  level = gene (e.g. mcrA) or subgroup (e.g. NiFe_1b, FeFe_B)
"""
import argparse, json, collections, statistics

ap = argparse.ArgumentParser()
ap.add_argument("--m8", required=True); ap.add_argument("--labels", required=True)
ap.add_argument("--fastp", required=True); ap.add_argument("--run", required=True)
ap.add_argument("--n_streamed", type=int, required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()

lab = {}
reflen_gene, reflen_sub = collections.defaultdict(list), collections.defaultdict(list)
with open(a.labels) as f:
    next(f)
    for line in f:
        sid, gene, guild, sub, mid, L, src = line.rstrip("\n").split("\t")
        lab[sid] = (gene, sub, float(mid))
        reflen_gene[gene].append(int(L))
        if sub:
            reflen_sub[f"{gene}_{sub}"].append(int(L))
reads_gene, reads_sub = collections.Counter(), collections.Counter()
for line in open(a.m8):
    q, s, pid, alen, qlen, slen, ev, bs = line.split("\t")[:8]
    gene, sub, mid = lab[s]
    if float(pid) < mid:
        continue
    reads_gene[gene] += 1
    if sub:
        key = f"{gene}_{sub}"
        reads_sub[key] += 1

qc = json.load(open(a.fastp))["summary"]["after_filtering"]["total_reads"]

# gene length = mean length (aa) of all reference sequences of the family / subgroup
def kb(lst):
    return statistics.mean(lst) * 3 / 1000

rp = lambda n, L: n / kb(L) / (qc / 1e6) if qc else float("nan")
rpkm_rpob = rp(reads_gene["rpoB"], reflen_gene["rpoB"])
with open(a.out, "w") as o:
    o.write("run\tfeature\tlevel\treads\trpkm\tcpg\tqc_reads\tn_streamed\n")
    for g in sorted(set(reads_gene) | {"mcrA", "acsB", "rpoB"}):
        r = rp(reads_gene[g], reflen_gene[g])
        cpg = r / rpkm_rpob if rpkm_rpob else float("nan")
        o.write(f"{a.run}\t{g}\tgene\t{reads_gene[g]}\t{r:.6g}\t{cpg:.6g}\t{qc}\t{a.n_streamed}\n")
    for k in sorted(reads_sub):
        r = rp(reads_sub[k], reflen_sub[k])
        cpg = r / rpkm_rpob if rpkm_rpob else float("nan")
        o.write(f"{a.run}\t{k}\tsubgroup\t{reads_sub[k]}\t{r:.6g}\t{cpg:.6g}\t{qc}\t{a.n_streamed}\n")
