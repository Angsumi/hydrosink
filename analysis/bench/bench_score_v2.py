#!/usr/bin/env python3
"""Benchmark scoring v2 (POST-HOC revision of the truth labels, written after the v1 audit table showed that the
pre-written product-name patterns missed RefSeq wording, e.g. McrA = "coenzyme-B sulfoethylthiotransferase subunit
alpha"). Changes vs v1 (bench_score.py, kept and reported):
  (1) a read is TRUE for gene g if ANY CDS it overlaps by >= 75 bp is a g source (v1: only the best-overlapping CDS);
  (3) sensitivity counts only reads with >= 120 bp inside a source CDS (the 80% query-cover filter needs >= 120 bp);
  (2) a CDS is a g source if its GFF gene symbol is in SYM[g], or (no symbol) its product matches PROD[g].
Misassignments to homologous but different enzymes (SdhA, L-aspartate oxidase, formate dehydrogenase, other
molybdopterin enzymes, McrG, homing endonucleases, ...) remain FALSE. Outputs: results/bench/*_v2.tsv"""
import bisect, collections, gzip, pathlib, re, urllib.parse
import pandas as pd

B = pathlib.Path(__file__).resolve().parent; A = B.parent; R = A / "results/bench"
GENES = ["mcrA", "acsB", "dsrA", "asrA", "nrfA", "narG", "napA", "frdA", "rpoB"]
SYM = {"mcrA": {"mcrA"}, "acsB": {"acsB", "cdhC"}, "dsrA": {"dsrA"}, "asrA": {"asrA"}, "nrfA": {"nrfA"},
       "narG": {"narG", "narZ"}, "napA": {"napA"}, "frdA": {"frdA"}, "rpoB": {"rpoB"}}
PROD = {k: re.compile(v, re.I) for k, v in {
    "mcrA": r"coenzyme-B sulfoethylthiotransferase subunit alpha|methyl-coenzyme M reductase subunit alpha",
    "acsB": r"acetyl-CoA decarbonylase/synthase complex subunit (alpha/)?beta|acetyl-CoA synthase complex subunit beta|CO dehydrogenase/acetyl-CoA synthase subunit beta",
    "dsrA": r"dissimilatory-type sulfite reductase subunit alpha|sulfite reductase, dissimilatory-type subunit alpha",
    "asrA": r"anaerobic sulfite reductase subunit (A|AsrA)\b",
    "nrfA": r"cytochrome c nitrite reductase( subunit c552| cytochrome c552 subunit)?$|ammonia-forming cytochrome c nitrite reductase( subunit c552)?$",
    "narG": r"nitrate reductase (\d )?(subunit )?alpha|nitrate reductase \d subunit alpha|nitrate reductase \d alpha subunit",
    "napA": r"periplasmic nitrate reductase( subunit alpha| NapA)?$|nitrate reductase catalytic subunit NapA",
    "frdA": r"fumarate reductase( flavoprotein subunit)?$|fumarate reductase catalytic and NAD/flavoprotein subunit",
    "rpoB": r"DNA-directed RNA polymerase subunit (beta(?!')|beta/beta'|B'{0,2}$)|RNA polymerase beta subunit"}.items()}
lab = pd.read_csv(A / "db/hydrosink.labels.tsv", sep="\t"); GENE_OF = dict(zip(lab.seq_id, lab.gene)); THR = dict(zip(lab.seq_id, lab.min_id))
acc = pd.read_csv(R / "genome_accessions.tsv", sep="\t").dropna(subset=["accession"])
_m = {}
def src(g, sym, prod):
    k = (g, sym, prod)
    if k not in _m: _m[k] = (sym in SYM[g]) if sym else bool(PROD[g].search(prod))
    return _m[k]

def load_cds(gff):
    cds = collections.defaultdict(list)
    for line in open(gff):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        p = re.search(r"product=([^;]+)", f[8]); s = re.search(r"(?:^|;)gene=([^;]+)", f[8])
        cds[f[0]].append((int(f[3]), int(f[4]), s.group(1) if s else "", urllib.parse.unquote(p.group(1)) if p else "unnamed"))
    return {c: ([x[0] for x in sorted(v)], sorted(v)) for c, v in cds.items()}

def overlaps(name, cds, min_ov=75):
    contig, s, e, *_ = name.rsplit("_", 5); s, e = int(s), int(e)
    if contig not in cds: return []
    starts, v = cds[contig]; i = bisect.bisect_right(starts, e)
    return [(sym, p) for a, b, sym, p in v[max(0, i - 50):i] if min(b, e) - max(a, s) + 1 >= min_ov]

tot = {g: collections.Counter() for g in GENES}; prods = collections.Counter(); host = []
for _, r in acc.iterrows():
    a = r.accession; hd = r.group == "host_decoy"
    assigned = {}
    with gzip.open(R / f"runs/{a}.m8.gz", "rt") as f:
        for line in f:
            q, s, pid = line.split("\t")[:3]; g = GENE_OF[s]
            if g in GENES and float(pid) >= THR[s]: assigned[q.split()[0]] = g
    with gzip.open(R / f"runs/{a}.names.gz", "rt") as f: names = [x.strip() for x in f]
    if hd:
        for g in GENES:
            n = sum(1 for x in assigned.values() if x == g); host.append(dict(gene=g, assigned=n, per_million=n / len(names) * 1e6))
        continue
    cds = load_cds(B / "genomes" / a / "genomic.gff")
    carrier = {g: any(src(g, sym, p) for _, v in cds.values() for _, _, sym, p in v) for g in GENES}
    for name in names:   # sensitivity denominator: >= 120 bp inside a source CDS (80% query cover needs >= 120 bp)
        ov = overlaps(name, cds, 120)
        for g in GENES:
            if any(src(g, sym, p) for sym, p in ov):
                tot[g]["truth"] += 1; tot[g]["detected"] += assigned.get(name) == g
    for q, g in assigned.items():
        ov = overlaps(q, cds); ok = any(src(g, sym, p) for sym, p in ov)
        tot[g]["assigned"] += 1; tot[g]["tp"] += ok
        if not carrier[g]: tot[g]["fp_noncarrier"] += 1
        top = max(ov, key=lambda x: 0, default=("", "intergenic"))
        prods[(g, top[0], top[1], ok)] += 1
    for g in GENES:
        if carrier[g]: tot[g]["carrier_genomes"] += 1
        else: tot[g]["noncarrier_reads"] += len(names)
    print(a, r.species, "done", flush=True)

met = pd.DataFrame([dict(gene=g, carrier_genomes=t["carrier_genomes"], assigned=t["assigned"], true_assigned=t["tp"],
                         precision=t["tp"] / t["assigned"] if t["assigned"] else float("nan"), truth_reads=t["truth"],
                         sensitivity=t["detected"] / t["truth"] if t["truth"] else float("nan"), fp_noncarrier_reads=t["fp_noncarrier"],
                         fp_per_million_noncarrier=t["fp_noncarrier"] / t["noncarrier_reads"] * 1e6 if t["noncarrier_reads"] else float("nan"))
                    for g, t in tot.items()])
met.to_csv(R / "bench_metrics_v2.tsv", sep="\t", index=False)
pd.DataFrame([dict(gene=g, gene_symbol=s, source_product=p, counted_true=ok, reads=n) for (g, s, p, ok), n in prods.items()]) \
    .sort_values(["gene", "reads"], ascending=[True, False]).to_csv(R / "bench_assigned_products_v2.tsv", sep="\t", index=False)
pd.DataFrame(host).to_csv(R / "bench_host_decoy_v2.tsv", sep="\t", index=False)
print(met.round(4).to_string(index=False))
