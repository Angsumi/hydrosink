# In-silico marker specificity benchmark — protocol (written 2026-11-17, before any benchmark read was simulated)

This benchmark was planned in amendment v1.1 §5 and executed after the main analyses. It estimates marker assignment performance on the selected RefSeq panel.

## Genome panel
NCBI RefSeq assemblies (reference genome if available, otherwise the RefSeq assembly with the highest contig N50) for
the species in `bench/genome_panel.tsv`: gut-associated bacteria and archaea spanning the H2-sink guilds
(methanogens, Wood–Ljungdahl carriers, sulfate/sulfite reducers, nitrate/nitrite/fumarate reducers) and common
fermenters without these pathways, plus one host genome (*Bos taurus*) as a decoy. The panel was chosen for
diversity; truth labels never depend on the panel list (see below). Accessions are recorded in
`results/bench/genome_accessions.tsv`.

## Reads
`wgsim` (seed 20261117), read 1 only, 150 bp, base error rate 0.005, no indels beyond wgsim defaults
(`-e 0.005 -r 0.001 -R 0.15 -1 150 -2 150 -d 150 -s 0 -S 20261117`; fragment = read length, so each read name gives its exact source coordinates). 500,000 reads per microbial genome, 2,000,000 for the host.

## Processing (identical to the study pipeline)
fastp 0.23.4 (`-q 20 -l 100 --trim_poly_g --cut_tail`); DIAMOND 2.1.10 blastx with the study parameters and database
(`--fast --id 40 --query-cover 80 --evalue 1e-5 --max-target-seqs 1 --max-hsps 1 -b 0.4 -c 4`); a read is assigned
to gene g if its best hit is a g reference with identity >= the per-gene threshold (as `aggregate_markers.py`).

## Original truth definition (v1)
The wgsim read name gives the source contig and coordinates. Each read is labelled with the NCBI RefSeq annotation
(GFF3 `product`) of the CDS it overlaps by >= 75 bp (or "intergenic"). For each marker g, a CDS product is a TRUE
source if it matches the regular expression in `TRUTH` (bench_score.py), defined from protein names before scoring:
  mcrA  : methyl-coenzyme M reductase subunit alpha
  acsB  : acetyl-CoA decarbonylase/synthase complex subunit (alpha/)beta | acetyl-CoA synthase | CO dehydrogenase/acetyl-CoA synthase subunit beta
  dsrA  : dissimilatory-type sulfite reductase subunit alpha | sulfite reductase, dissimilatory-type subunit alpha
  asrA  : anaerobic sulfite reductase subunit A
  nrfA  : ammonia-forming nitrite reductase cytochrome c552 subunit | cytochrome c nitrite reductase | nitrite reductase \(cytochrome; ammonia-forming\) c552 subunit
  narG  : nitrate reductase subunit alpha
  napA  : periplasmic nitrate reductase subunit alpha | nitrate reductase catalytic subunit NapA
  frdA  : fumarate reductase flavoprotein subunit
  rpoB  : DNA-directed RNA polymerase subunit beta(?!')
Because annotation wording varies, every product name behind every assigned read is also tabulated
(`bench_assigned_products.tsv`), so the truth calls can be audited and any regex miss is visible.

## Original metrics (v1), per marker
- precision = assigned reads from TRUE-source CDS / all assigned reads
- false-positive rate = assigned reads from genomes with no TRUE-source CDS, per million simulated reads
- sensitivity = reads overlapping TRUE-source CDS (>= 75 bp) that were assigned to g / all such reads
- host decoy: assigned reads per million host reads, per marker
## Scoring revision (v2)

The post hoc scorer (`bench_score_v2.py`) first checks RefSeq gene symbols and then uses revised product-name patterns. For precision, any overlapping target CDS with at least 75 bp of overlap qualifies. For sensitivity, the denominator is reads with at least 120 bp inside a target CDS. The two scorers also differ in handling overlapping genes, so differences between versions cannot be attributed to annotation wording alone. See `analysis/results/bench/bench_metrics.tsv` (v1) and `bench_metrics_v2.tsv` (v2).
