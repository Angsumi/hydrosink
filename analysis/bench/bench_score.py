#!/usr/bin/env python3
"""Benchmark step 3: score marker assignment against read origin (see BENCHMARK_PROTOCOL.md).
Inputs: results/bench/genome_accessions.tsv, results/bench/runs/<acc>.{m8.gz,names.gz}, bench/genomes/<acc>/genomic.gff
Outputs (results/bench/): bench_metrics.tsv, bench_per_genome.tsv, bench_assigned_products.tsv, bench_host_decoy.tsv
"""
import bisect, collections, gzip, pathlib, re, urllib.parse
import pandas as pd

B = pathlib.Path(__file__).resolve().parent; A = B.parent; R = A / "results/bench"
GENES = ["mcrA", "acsB", "dsrA", "asrA", "nrfA", "narG", "napA", "frdA", "rpoB"]
TRUTH = {
    "mcrA": r"methyl-coenzyme M reductase subunit alpha",
    "acsB": r"acetyl-CoA decarbonylase/synthase complex subunit (alpha/)?beta|acetyl-CoA synthase|CO dehydrogenase/acetyl-CoA synthase subunit beta",
    "dsrA": r"dissimilatory-type sulfite reductase subunit alpha|sulfite reductase, dissimilatory-type subunit alpha",
    "asrA": r"anaerobic sulfite reductase subunit A",
    "nrfA": r"ammonia-forming nitrite reductase cytochrome c552 subunit|cytochrome c nitrite reductase|nitrite reductase \(cytochrome; ammonia-forming\) c552 subunit",
    "narG": r"nitrate reductase subunit alpha",
    "napA": r"periplasmic nitrate reductase subunit alpha|nitrate reductase catalytic subunit NapA",
    "frdA": r"fumarate reductase flavoprotein subunit",
    "rpoB": r"DNA-directed RNA polymerase subunit beta(?!')",
}
TRUTH = {g: re.compile(p, re.I) for g, p in TRUTH.items()}
_memo = {}


def is_true(g, product):
    k = (g, product)
    if k not in _memo: _memo[k] = bool(TRUTH[g].search(product))
    return _memo[k]
lab = pd.read_csv(A / "db/hydrosink.labels.tsv", sep="\t")
GENE_OF = dict(zip(lab.seq_id, lab.gene)); THR = dict(zip(lab.seq_id, lab.min_id))
acc = pd.read_csv(R / "genome_accessions.tsv", sep="\t").dropna(subset=["accession"])


def load_cds(gff):
    cds = collections.defaultdict(list)
    for line in open(gff):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        prod = re.search(r"product=([^;]+)", f[8]); prod = urllib.parse.unquote(prod.group(1)) if prod else "unnamed"
        cds[f[0]].append((int(f[3]), int(f[4]), prod))
    out = {}
    for c, v in cds.items():
        v.sort(); out[c] = ([x[0] for x in v], v)
    return out


def origin(name, cds):
    contig, s, e, *_ = name.rsplit("_", 5)
    s, e = int(s), int(e)
    if contig not in cds: return "intergenic"
    starts, v = cds[contig]; i = bisect.bisect_right(starts, e)
    best, ov = "intergenic", 0
    for j in range(max(0, i - 50), i):          # CDS starting before read end (genes < ~50 overlap)
        a, b, p = v[j]; o = min(b, e) - max(a, s) + 1
        if o >= 75 and o > ov: best, ov = p, o
    return best


per_genome, products, host = [], collections.Counter(), []
tot = {g: collections.Counter() for g in GENES}   # tp, assigned, fp_noncarrier, truth_reads, truth_detected
for _, r in acc.iterrows():
    a = r.accession; m8 = R / f"runs/{a}.m8.gz"
    if not m8.exists(): print("missing run", a); continue
    host_decoy = r.group == "host_decoy"
    cds = {} if host_decoy else load_cds(B / "genomes" / a / "genomic.gff")
    assigned = {}
    with gzip.open(m8, "rt") as f:
        for line in f:
            q, s, pid = line.split("\t")[:3]
            g = GENE_OF[s]
            if g in GENES and float(pid) >= THR[s]: assigned[q.split()[0]] = g
    n_sim = 0; truth_names = collections.defaultdict(set); carrier = collections.Counter()
    if not host_decoy:
        for _, (_, v) in cds.items():
            for _, _, p in v:
                for g in GENES:
                    if TRUTH[g].search(p): carrier[g] += 1
        with gzip.open(R / f"runs/{a}.names.gz", "rt") as f:
            for name in f:
                name = name.strip(); n_sim += 1; o = origin(name, cds)
                for g in GENES:
                    if is_true(g, o): truth_names[g].add(name)
    else:
        with gzip.open(R / f"runs/{a}.names.gz", "rt") as f: n_sim = sum(1 for _ in f)
    row = dict(accession=a, species=r.species, group=r.group, simulated=n_sim)
    for g in GENES:
        qs = [q for q, x in assigned.items() if x == g]
        if host_decoy:
            host.append(dict(gene=g, assigned=len(qs), per_million=len(qs) / n_sim * 1e6)); continue
        tp = 0
        for q in qs:
            o = origin(q, cds); ok = is_true(g, o); tp += ok
            products[(g, o, ok)] += 1
        row[f"{g}_carrier_cds"] = carrier[g]; row[f"{g}_assigned"] = len(qs); row[f"{g}_tp"] = tp
        row[f"{g}_truth_reads"] = len(truth_names[g]); row[f"{g}_truth_detected"] = sum(1 for q in truth_names[g] if assigned.get(q) == g)
        t = tot[g]; t["assigned"] += len(qs); t["tp"] += tp; t["truth"] += len(truth_names[g]); t["detected"] += row[f"{g}_truth_detected"]
        if carrier[g] == 0: t["fp_noncarrier"] += len(qs); t["noncarrier_reads"] += n_sim
        else: t["carrier_genomes"] += 1
    if not host_decoy: per_genome.append(row)
    print(a, r.species, "done", flush=True)

pd.DataFrame(per_genome).to_csv(R / "bench_per_genome.tsv", sep="\t", index=False)
met = pd.DataFrame([dict(gene=g, carrier_genomes=t["carrier_genomes"], assigned=t["assigned"], true_assigned=t["tp"],
                         precision=t["tp"] / t["assigned"] if t["assigned"] else float("nan"),
                         truth_reads=t["truth"], sensitivity=t["detected"] / t["truth"] if t["truth"] else float("nan"),
                         fp_noncarrier_reads=t["fp_noncarrier"],
                         fp_per_million_noncarrier=t["fp_noncarrier"] / t["noncarrier_reads"] * 1e6 if t["noncarrier_reads"] else float("nan"))
                    for g, t in tot.items()])
met.to_csv(R / "bench_metrics.tsv", sep="\t", index=False)
pd.DataFrame([dict(gene=g, source_product=p, counted_true=ok, reads=n) for (g, p, ok), n in products.items()]) \
    .sort_values(["gene", "reads"], ascending=[True, False]).to_csv(R / "bench_assigned_products.tsv", sep="\t", index=False)
pd.DataFrame(host).to_csv(R / "bench_host_decoy.tsv", sep="\t", index=False)
print(met.round(4).to_string(index=False)); print(pd.DataFrame(host).to_string(index=False))
