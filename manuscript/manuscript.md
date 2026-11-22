# Hindgut enrichment of acetogenesis markers in ruminants and faecal variation across vertebrates

**Ananda Prasad**¹ (ORCID [0009-0005-1736-8148](https://orcid.org/0009-0005-1736-8148)) and **Angsuman Das**²\* (ORCID [0000-0001-8146-7697](https://orcid.org/0000-0001-8146-7697))

¹ Department of Zoology, Rajiv Gandhi University, Itanagar, Arunachal Pradesh, India. anandaprasad01@gmail.com
² Department of Zoology, Pachhunga University College (a constituent college of Mizoram University), Aizawl, Mizoram, India. angsuman.das@rgu.ac.in

\*Corresponding author: Angsuman Das, angsuman.das@rgu.ac.in

*Draft v2.3 (2026-11-22). The figures and result tables cited below are available in the repository. Reporting checklists are in `analysis/reporting/`, and reference codes [Lxx] identify records in `literature.csv`.*

## Abstract

**Background.** Microbial consumption of fermentative H₂ can support methanogenesis, reductive acetogenesis and respiration. We asked whether the relative abundance of marker genes for these pathways varies among gut compartments in ruminants or among vertebrates with different fermentation sites.

**Methods.** We quantified marker genes in public shotgun metagenomes using five million forward reads per sequenced run. The primary outcome was the pseudocount-adjusted log ratio of the acetyl-CoA synthase marker *acsB* to the methyl-coenzyme M reductase marker *mcrA*. We used DIAMOND blastx to assign reads to markers. We pooled paired compartment contrasts from 53 ruminants in two studies and fitted Bayesian phylogenetic mixed models to 802 faecal samples from 166 vertebrate species in 179 studies.

**Results.** The ratio was higher in the hindgut in all 53 ruminants. The pooled hindgut-minus-foregut difference was 2.62 natural-log units (95% CI 2.23–3.01), driven principally by higher *acsB* abundance. For the faecal samples, adding fermentation site to the model changed the deviance information criterion by −0.08 (reduced minus full model). The foregut-minus-hindgut contrast was 0.47 log units (95% highest posterior density interval −0.95 to 1.79). Phylogeny accounted for an estimated 39% of the total random-effect and residual variance, rising to 78% with an alternative dated tree. Foregut and hindgut hosts had higher bacterial-*rpoB*-normalised abundances of both markers than hosts classified as simple-gut. Pooled within-study proportionality was close to zero (ρp = 0.03, 95% CI −0.10 to 0.17), with substantial variation among study units.

**Conclusions.** The paired ruminant data confirm distal enrichment of *acsB* relative to *mcrA*. Fermentation site explained little of the faecal log-ratio variation in the fitted cross-host models. The two comparisons represent different anatomical contrasts. DNA marker abundances capture genetic potential with no direct measurement of pathway activity or hydrogen flux.

## Introduction
Fermentation in vertebrate guts releases molecular hydrogen (H₂). If H₂ accumulates, re-oxidation of reduced cofactors becomes thermodynamically unfavourable and fermentation slows, so hydrogenotrophic microbes that keep H₂ partial pressure low are central to gut function [L03]. Several guilds compete for this H₂. Hydrogenotrophic methanogens reduce CO₂ to methane, an energy loss to the host and a greenhouse gas. Reductive acetogens produce acetate, which the host can absorb. Respiratory hydrogenotrophs couple H₂ oxidation to the reduction of sulfate, nitrate, nitrite or fumarate [L03]. In the rumen, H₂ uptake is distributed among these guilds, and alternative uptake pathways are expressed more in low-methane sheep [L03].

Which sink dominates appears to depend on where fermentation happens. Li et al. reported higher caecal than ruminal *acsB* abundance in goats and in a comparative analysis that used the ruminant data of Xie et al. [L02, L04]. Marsupials, in contrast, are heterogeneous: some species resemble high-methane mammals while others are enriched for respiratory H₂-uptake potential [L01]. Across 110 vertebrate species, archaeal diversity tracks host phylogeny more than diet, but that analysis used 16S amplicons and did not measure H₂-sink genes [L07]. Large shotgun catalogues of vertebrate guts now exist [L04, L05, L08], as do curated cross-study metadata [L10]. To our knowledge, no study has tested with phylogenetic and study-level control whether gut fermentation architecture predicts how H₂-disposal potential is partitioned across vertebrates.

We specified four predictions and their planned analyses in the repository before fitting the main models:

- **P1** (within animals; positive control): in paired foregut and hindgut samples from the same animal, ln(*acsB*/*mcrA*) is higher in the hindgut.
- **P2** (across hosts): fermentation site (foregut, hindgut, avian caecal, simple gut) explains variation in faecal ln(*acsB*/*mcrA*) beyond host phylogeny, diet, captivity and technical covariates.
- **P3** (trade-off): within studies, methanogenic and acetogenic potential are negatively proportional.
- **P4** (exploratory): respiratory H₂-uptake potential is higher where methanogenic potential is lower.

P1 reanalyses the two ruminant data sources used by Li et al. [L02, L04] with a common marker pipeline and paired effect estimates. The overlap in data sources makes P1 a positive control rather than an independent replication. P2–P4 extend the question to cross-host and within-study comparisons.

## Methods
### Protocol, registration and amendments
The initial protocol, thresholds and analysis code were committed to the repository in `67eb757` and `eec1a3e` on 15 November 2026. Git tag `prereg-v1` points to the later Gate 1a positive-control commit. Two amendments were committed before P2 marker data were generated, following the positive-control and sampling checks described below.

- **v1.1** (`analysis/PREREG_AMENDMENT_v1.1.md`) designated P1 a positive control, added per-study estimates and a three-level model, introduced four P2 sensitivity or validation analyses, and documented outcome-blind captivity curation.
- **v1.2** (`analysis/PREREG_AMENDMENT_v1.2.md`) changed P2 read sampling after the head-sampling check failed.

The dated commits are the study's protocol record. There was no external registration.

### Data sources and screening
The European Nucleotide Archive (ENA) Portal API was queried on 15 November 2026 for Illumina WGS runs with library source METAGENOMIC and one of seven gut-associated metagenome taxa: gut, feces, rumen, bovine gut, pig gut, chicken gut and fish gut metagenome. Figure S2 and its source tables report screening in a PRISMA 2020 flow adapted for repository data [L31]. A separate checklist adapts STORMS for this secondary analysis [L32]. The screen applied the following criteria.

- **Host:** a structured host taxon ID within Vertebrata, excluding human, mouse and rat.
- **Keywords:** no keywords indicating treatments, inhibitors, in-vitro systems, controls or non-WGS libraries.
- **Read count:** ≥ 5 × 10⁶ reads.
- **Read length:** mean ≥ 100 bp.
- **One run per BioSample:** the deepest.

For P2 we kept faecal samples and capped each host species at 10 runs, drawn round-robin across studies with a fixed seed (20261115). For P1 we selected two studies with paired fermentation-chamber samples from the same animals: Li et al. 2025 (PRJNA1068803; 16 goats, rumen and caecum) [L02] and Xie et al. 2021 (PRJNA657455; seven ruminant species) [L04]. P1 included the rumen, reticulum and omasum as foregut sites and the caecum, colon and rectum as hindgut sites.

### Host traits
**Fermentation site** was assigned from host lineage by fixed rules following Stevens & Hume [L20] (`analysis/scripts/01b_host_traits.py`):

- **Foregut:** Ruminantia, Camelidae, Hippopotamidae, Tayassuidae, Macropodidae and Potoroidae.
- **Hindgut:** Equidae, Rhinocerotidae, Tapiridae, Elephantidae, Sirenia, Lagomorpha, Rodentia, wombats, koala and phalangeriform marsupials, Hyracoidea, Suidae, non-colobine primates and herbivorous Testudines.
- **Avian caecal:** Galliformes, Anseriformes, ratites, Tinamiformes, Gruiformes and Otidiformes.
- **Simple:** all other vertebrates.

**Diet class** (herbivore, omnivore or carnivore) came from EltonTraits 1.0 [L21]; hosts without an EltonTraits entry were represented by a separate missing-diet category in the fitted model (145 samples).

**Captivity** (wild, captive, domestic, laboratory or unknown) was curated from study and BioSample text by keyword rules plus 38 manual study-level decisions, each with an evidence quote (`analysis/config/captivity_manual.tsv`). P2 marker data were unavailable during curation.

**Host phylogeny** was the Open Tree of Life synthetic tree [L34] induced for the matched taxa (`01c_host_tree.py`). Because that tree is topology-only, Grafen branch lengths [L35] were added (ape `compute.brlen`). Hosts identified only to genus without a tree tip were excluded from the phylogenetic models.

### Read sampling
Bandwidth precluded full downloads, so each run was represented by 5 × 10⁶ forward (R1) reads; R1 only was used, as in the Greening-lab protocol [L23]. **P1** used the first 5 × 10⁶ reads ("head" sampling).

A pre-specified gate (Gate 1b) compared head sampling against a uniform random subsample of the full file in two goat samples. For gene groups with at least 100 reads in both subsamples, the rule required the *rpoB*-normalised abundance estimates to agree within 10%. Lower counts were evaluated with a Poisson test. One of ten checks failed. The [NiFe]-hydrogenase group in one rumen sample had 210 versus 189 assigned reads and a 10.3% difference in normalised abundance [`gate1b_head_vs_random.tsv`]. As pre-specified, P2 therefore switched to **multi-offset sampling** (amendment v1.2): 10 blocks of 5 × 10⁵ consecutive spots, with block *k* starting at spot ⌊*k*·*S*/10⌋ + 1, where *k* ranges from 0 to 9 and *S* is the run's spot count.

Reads came from the NCBI SRA copy on AWS. A run was accepted only if every block returned exactly the planned rows. Three fetch implementations were used: remote `vdb-dump`, local extraction from the downloaded `.sra` object, and `fastq-dump`. They were verified to give identical bases and qualities (amendment v1.2 addenda).

### Quality control and marker quantification

- **Quality control:** fastp 0.23.4 [L27] (`-q 20 -l 100 --trim_poly_g --cut_tail`).
- **Search:** DIAMOND 2.1.10 blastx [L26] (`--fast --id 40 --query-cover 80 --evalue 1e-5 --max-target-seqs 1 --max-hsps 1`).
- **Reference:** the Greening-lab compiled metabolic marker database v3 (51 genes, all retained as competitors) [L22]. To this were added 5,635 full-length UniRef90 bacterial RpoB sequences (UniProt 2026_04) [L25] as the single-copy normaliser, and RpoC sequences as decoys (170,171 sequences in total; md5 values in `config.yaml`).
- **Assignment:** a read was assigned to the gene of its best hit if percent identity met that gene's published read-level threshold [L23]. This was 50% for most genes, including *mcrA*, *acsB*, *dsrA*, *nrfA*, *narG*, *napA* and *frdA*, and 60% for [FeFe]- and group-4 [NiFe]-hydrogenases.
- **Abundance:** RPKM = reads / (mean reference length in kb) / (QC reads in millions). Marker abundance relative to bacterial *rpoB* (cpg in the output tables) equals RPKM(gene) divided by RPKM(*rpoB*) (`aggregate_markers.py`). This denominator approximates bacterial genome equivalents and also applies to archaeal markers such as *mcrA*.
- **Sample QC:** a sample passed if it had ≥ 200 *rpoB* reads.

The DIAMOND block size was reduced from 1 to 0.4 during P2 to bound memory. A test run gave identical hit files (amendment v1.2 addendum). Processing used one workstation and four Google Cloud VMs with identical tool versions. A run processed on both gave byte-identical marker tables.

### Outcomes

- **Primary:** ALR_acet = ln[(*acsB* + *c*ₐ)/(*mcrA* + *c*ₘ)]. Each pseudocount is half the smallest positive abundance of its respective marker in the quality-filtered input for that analysis.
- **Secondary:**
  - ALR_resp uses the sum of the six respiratory reductase marker abundances (*dsrA*, *asrA*, *nrfA*, *narG*, *napA*, *frdA*) and a pseudocount calculated for that sum in the numerator;
  - log cpg of each marker;
  - group-1 [NiFe]- and group-B [FeFe]-hydrogenases;
  - *mcrA* detection (≥ 2 reads).

### Statistical analysis

A generative language model assisted with analysis-code development and execution under the authors’ direction (see Declarations).
The statistical analyses ran in R 4.3.3 with metafor 4.8-0 [L29], lme4 1.1-37 and MCMCglmm 2.36 [L28] (`05_meta_analysis.R`; seed 20261115).

**P1.** For each animal we computed the hindgut-minus-foregut difference of chamber means. Differences were averaged within each study × host-species stratum, and strata were pooled by random-effects meta-analysis (REML τ², Hartung–Knapp–Sidik–Jonkman CIs [L30], 95% prediction intervals, leave-one-out). Because P1 has only two source studies, we also report each study separately and a three-level model (strata nested in study). A one-stage linear mixed model was fitted as a cross-check: random animal intercepts and random chamber slopes by stratum. Secondary outcomes were adjusted by Benjamini–Hochberg.

**P2.** We fitted Gaussian MCMCglmm models with random effects for phylogeny (inverse relatedness matrix), host species and study. Fixed effects were diet class, captivity, log₁₀ *rpoB* reads and read length, plus fermentation site in the full model. Chains ran for 130,000 iterations, with 30,000 burn-in and thinning of 100. Priors were parameter-expanded for the random effects (V = 1, ν = 1, α.V = 1000) and inverse-Wishart for the Gaussian residual (V = 1, ν = 0.002). The categorical model fixed the residual variance at 1. Models were compared using the deviance information criterion (ΔDIC = reduced-model DIC minus full-model DIC). We report posterior means and 95% highest posterior density (HPD) intervals for site contrasts, and posterior modes and 95% HPD intervals for the variance shares. We fitted a separate categorical model for *mcrA* detection with the same fixed and random effects.

Pre-specified P2 sensitivity analyses:

- mammals only;
- non-domestic hosts;
- ≤ 3 samples per species;
- ≥ 1,000 *rpoB* reads;
- herbivores only;
- known captivity;
- a sequencing-platform covariate;
- a within-study analysis restricted to studies sampling ≥ 2 fermentation sites, with study as a fixed effect;
- a marsupial hold-out: fit without Metatheria, then predict marsupial samples.

**P3.** Within each quality-filtered P2 study or P1 study × chamber unit with at least 10 samples, we calculated symmetric proportionality (ρp) between *mcrA* and *acsB* after centred log-ratio transformation of an 11-component profile (*rpoB* and ten markers). Before transformation, each component received a unit-specific pseudocount equal to half its smallest positive value, or 10⁻⁶ if entirely zero. We clipped ρp to ±0.999 for Fisher z transformation, estimated variances from 1,000 sample-level bootstrap replicates, and pooled units with positive finite bootstrap variances using REML with Hartung–Knapp intervals. We back-transformed the pooled estimates and intervals, and repeated pooling for units with *mcrA* detected in at least 50% of samples. There was no additional adjustment for host or animal within units [L33].

**Sensitivity of thresholds.** Marker counts were re-derived from the saved DIAMOND hits with all identity thresholds shifted by −10 and +10 percentage points (floor 40%, the search cut-off), and the P1, P3 and primary P2 analyses were repeated (`05c_identity_sensitivity.py`, `05d_sensitivity_sets.R`). With a shift of 0 the re-derivation reproduced all 1,163 original marker tables byte for byte.

***acsB* attribution (post-hoc).** The pre-planned phylogenetic placement of *acsB* reads was not possible because reads were not retained. Instead, each *acsB* read was attributed to the NCBI taxonomy (E-utilities) of its best-hit reference protein. P1, P3 and the primary P2 model were repeated with *acsB* restricted to bacterial or to Bacillota references (`08a`–`08c` scripts).

**Dated-tree sensitivity.** The Open Tree topology was kept (polytomies resolved). Each internal node was dated with the TimeTree 5 API [L36]: the MRCA age of up to 10 analysis species from each child clade, using NCBI taxon IDs. Nodes without an age, or with an age older than their parent, were interpolated between parent and oldest dated descendant. Each node was then constrained to ≤ 98% of its parent's age, giving an ultrametric tree. The primary P2 model was refitted on this tree (`05f_dated_tree.R`).

**In-silico specificity benchmark.** Reads were simulated from 43 RefSeq reference genomes and one host decoy (*Bos taurus*), retrieved with NCBI Datasets [L37]. The genomes were chosen for guild diversity: methanogens, Wood–Ljungdahl carriers, sulfate/sulfite reducers, nitrate/fumarate reducers and common fermenters. Settings: wgsim (seed 20261117); 150-bp reads with fragment = read length, so each read's source coordinates are known; 500,000 reads per microbial genome and 2,000,000 host reads. The reads were processed with the study's fastp, DIAMOND and threshold settings (`analysis/bench/`; protocol written before simulation, `BENCHMARK_PROTOCOL.md`).

Each assigned read was scored against the RefSeq annotation of the coding sequences it overlaps (≥ 75 bp):

- **Precision:** the fraction of assigned reads from a true source gene.
- **Sensitivity:** the fraction of reads with ≥ 120 bp inside a source gene that were assigned to it.
- **False positives per million reads:** counted in genomes lacking the gene.

The pre-written product-name truth patterns (v1) missed RefSeq wording; for example, McrA is annotated "coenzyme-B sulfoethylthiotransferase subunit alpha". Truth was therefore revised post hoc (v2) to use GFF gene symbols, falling back to corrected product patterns. Version 2 also allowed any CDS with at least 75 bp of overlap to count toward precision and required at least 120 bp within a source CDS for the sensitivity denominator. Both versions are reported. The benchmark was planned in amendment v1.1 and run after the main analyses as a methods check.

**P1 sampling check.** Five P1 runs drawn with a fixed seed after Gate 1a but before the main P1 analysis were re-quantified with multi-offset sampling and compared with head sampling (amendment v1.2 §2).

### Data and code availability
The [HydroSink repository](https://github.com/Angsumi/hydrosink) provides code, per-sample results and derived outputs. The analysed data and preceding manuscript draft were archived as [v2.2.0](https://doi.org/10.5281/zenodo.24331026). The [concept DOI](https://doi.org/10.5281/zenodo.24330934) identifies the archive across releases. Run accessions, curated metadata and per-sample marker tables are in `analysis/inventory/` and `analysis/results/per_sample/`; DIAMOND hits are in `analysis/results/hits/`, with SHA-256 manifests (`p1_manifest.sha256`, `p2_manifest.sha256`). Analysis and figure code is in `analysis/scripts/`. Benchmark code is in `analysis/bench/`, and `run_all.sh` reproduces the statistical outputs from the archived inputs.

## Results
### Data set
Screening reduced 136,708 candidate runs to 19,385 eligible runs, of which 8,293 were faecal. Capping at 10 runs per host species gave 930 P2 runs [`prisma_log.tsv`].

Of the 930 runs:

- 926 were quantified. The other 4 runs (one study) consisted entirely of 76-bp reads, below the length filter, and yielded no reads after QC. The screen's read-length estimate had counted both mates of each spot and so missed this.
- 845 passed the *rpoB* QC threshold.
- 802 samples from 166 host species and 179 studies entered the phylogenetic models [`prisma_final.tsv`]. Genus-only hosts without a tree tip were excluded at this step.

The modelled set comprised 207 foregut, 250 hindgut, 67 avian-caecal and 278 simple-gut samples; 655 were mammalian [`table1_p2.tsv`]. Captivity was wild for 274 samples, captive, domestic or laboratory for 263, and unknown for 265.

Read depth was uniform by design: median 4.96 × 10⁶ reads after QC, median 5,707 *rpoB* reads per sample. *mcrA* (≥ 2 reads) was detected in 88% of foregut, 60% of hindgut, 28% of avian-caecal and 12% of simple-gut samples; *acsB* in 97%, 98%, 78% and 67% [`detection_by_site.tsv`; Fig. S1].

Of the 254 P1 runs, 253 passed QC. Both compartments were represented in each of 53 animals: 16 goats from Li 2025 [L02] and 37 animals of seven ruminant species from Xie 2021 [L04] [`table1_p1.tsv`].

### P1 Paired gut compartments in ruminants
ALR_acet was higher in the hindgut than in the foregut in **53 of 53 animals**. The pooled stratum-level difference was **2.62 ln units** (95% CI 2.23–3.01; 95% prediction interval 1.69–3.55; *k* = 8 strata; I² = 70%) [`p1_pooled.tsv`; Fig. 2]. The pooled estimate corresponds to a 13.7-fold difference on the adjusted ratio scale. The median per-animal fold difference was 12.5 [`p1_per_animal.tsv`].

The study-specific estimates were 1.99 (95% CI 1.54–2.44) for Li 2025 and 2.71 (2.33–3.08) for Xie 2021. It was also present in the three-level model with strata nested in study: 2.43 (1.62–3.24) [`p1_per_study.tsv`]. Other checks:

- One-stage mixed model: 2.72 (2.39–3.06) [`p1_mixed_model.tsv`].
- Leave-one-stratum-out estimates: 2.56 to 2.71 [`p1_leave1out.tsv`].
- Rumen versus caecum only: 2.72. Rumen versus rectum: 2.80 [`p1_sensitivity.tsv`].

The contrast was driven mainly by *acsB*, which was higher in all 53 animals. The pooled difference was 2.92 ln units (95% CI 2.15–3.70). The estimate for *mcrA* was smaller and its interval overlapped zero: +0.25 (−0.31 to 0.80; *q* = 0.39). Among secondary markers [`p1_pooled.tsv`]:

- *dsrA* and group-B [FeFe]-hydrogenase were higher in the hindgut.
- *nrfA* (−2.92), *narG* (−1.00), group-1 [NiFe]-hydrogenase (−0.94) and the respiratory sum (−1.33) were lower.
- ALR_resp was therefore lower in the hindgut: −1.58 (−2.53 to −0.62).

The paired data show higher hindgut *acsB* abundance and a lower summed respiratory marker abundance. Several individual respiratory markers differed in direction.

The pre-specified head-versus-multi-offset sampling check was run on 5 P1 runs [`p1_offsets_vs_head.tsv`]. In 22 of 25 run × gene-group comparisons the 95% CI of the ratio included 1. In one rumen sample (SRR12529313), multi-offset sampling gave 32% lower *rpoB*-normalised *acsB* abundance than head sampling (offsets/head = 0.68, 95% CI 0.62–0.74). The difference was −0.39 ln units. This run-level check does not quantify sampling bias in the pooled P1 contrast.

### P2 Faecal marker ratios by fermentation site
In the phylogenetic mixed model, adding fermentation site did not improve fit for ALR_acet: ΔDIC = −0.08 [`p2_model_comparison.tsv`]. No site contrast had a 95% HPD interval excluding zero [`p2_contrasts.tsv`; Fig. 3]:

| Contrast | Posterior mean | 95% HPD |
|---|---|---|
| Foregut − hindgut | 0.47 | −0.95 to 1.79 |
| Hindgut − simple | −0.66 | −2.19 to 1.01 |
| Foregut − simple | −0.19 | −2.03 to 1.76 |
| Avian caecal − simple | 1.13 | −0.88 to 3.17 |

Of the total random-effect and residual variance in ALR_acet, the share attributable to host phylogeny was 0.39 (95% HPD 0.21–0.67), to species 0.15 and to study 0.16 [`p2_variance_partition.tsv`]. The reported intervals for diet, captivity, depth and read length overlapped zero [`p2_fixed_effects.tsv`].

The site contrast intervals included zero in the following pre-specified model refits [`p2_sensitivity.tsv`, `p2_sensitivity_v1_1.tsv`]:

- mammals only (626 samples);
- non-domestic (691);
- ≤ 3 samples per species (355);
- ≥ 1,000 *rpoB* reads (747);
- herbivores only (452);
- known captivity (537);
- with a sequencing-platform covariate.

Restricting to the 12 studies that sampled more than one fermentation site (136 samples, 72 species; study as a fixed effect) also gave no clear site effect: foregut −1.25 (SE 0.89), hindgut −0.96 (SE 0.76) relative to simple [`p2_within_study.tsv`]. The implied foregut-minus-hindgut point difference was −0.29 ln units. These faecal comparisons use different hosts, whereas P1 compares compartments within each animal.

In the marsupial hold-out (19 samples, 4 species), the model predicted no hindgut–foregut difference (+0.22). The observed difference was −3.38, but it rests on a single foregut sample. The species-level Spearman correlation was 0.4 on 4 species [`p2_holdout_marsupial.tsv`]; the small hold-out and single foregut sample limit the assessment.

Foregut and hindgut hosts had higher *rpoB*-normalised abundances of both markers than simple-gut hosts [`p2_contrasts.tsv`]:

| Contrast | ln *acsB* (95% HPD) | ln *mcrA* (95% HPD) |
|---|---|---|
| Hindgut − simple | +3.04 (1.69–4.56) | +3.58 (2.44–4.68) |
| Foregut − simple | +3.91 (2.26–5.36) | +3.97 (2.45–5.32) |

Avian-caecal hosts had more *acsB* than simple-gut hosts (+2.08; 0.36–3.62), but *mcrA* did not differ credibly. In the separate detection model, *mcrA* was more often detected in all three fermenting categories than in simple-gut hosts [`p2_fixed_effects.tsv`].

At the two-read detection threshold, *mcrA* was undetected in 88% of simple-gut samples and neither marker was detected in 115 samples, including 90 from simple-gut hosts [`p2_analysis_set.tsv`]. The pseudocounts have greater influence at low counts and determine the ratio when both counts are zero. The separate abundance and detection results help interpret these samples.

**Post-hoc check (not pre-registered).** The ENA screen excluded house mouse by taxon 10090 only, so 20 *Mus musculus* samples entered P2 under a subspecies or genus taxon, 10 of them from a DSS-colitis model. In addition, one *Elephas maximus* sample lacked a diet label. Removing all *Mus musculus* and harmonising that label gave the same result: 782 samples, ΔDIC = −0.18, foregut − hindgut 0.45 (−0.98 to 1.94) [`posthoc_no_mouse/p2_primary.tsv`].

### P3 Within-study marker proportionality
Across 34 units (30 P2 studies with ≥ 10 samples and 4 P1 study × chamber units), the pooled proportionality between *mcrA* and *acsB* was ρp = 0.03 (95% CI −0.10 to 0.17; 95% prediction interval −0.56 to 0.61; I² = 75%) [`p3_pooled.tsv`; Fig. 4]. In units where *mcrA* was detected in ≥ 50% of samples (*k* = 18), it was −0.01 (−0.20 to 0.18). Unit estimates ranged from −0.80 to 0.72 [`p3_units.tsv`]. The estimates provide little evidence for a consistently negative association on the chosen compositional scale.

### P4 Exploratory results for respiratory markers
Across hosts, ALR_resp was much lower in hindgut (−5.01; −6.50 to −3.77) and foregut (−5.36; −6.92 to −3.69) than in simple-gut hosts [`p2_contrasts.tsv`]. The higher *mcrA* abundance in fermenting hosts contributes to the lower ratio. The ratio model alone cannot quantify a corresponding change in respiratory abundance or test whether the two components covary.

Within animals (P1), the hindgut had lower ALR_resp than the foregut at similar *mcrA* (above). The P4 prediction concerns inverse variation between the components. These ratio contrasts cannot establish that relationship. Cross-assignment of *napA* and *frdA* reads in the benchmark adds uncertainty to the respiratory sum.

### Sensitivity to marker identity thresholds
Shifting every per-gene identity threshold by −10 or +10 percentage points left all conclusions unchanged [`sens_id_m10/`, `sens_id_p10/`]:

| Threshold shift | P1 pooled ALR_acet (95% CI) | P2 ΔDIC | P2 foregut − hindgut (95% HPD) | P3 pooled ρp (95% CI) |
|---|---|---|---|---|
| −10 | 2.61 (2.22–3.00) | −0.17 | 0.48 (−0.86 to 2.02) | 0.01 (−0.12 to 0.15) |
| 0 (primary) | 2.62 (2.23–3.01) | −0.08 | 0.47 (−0.95 to 1.79) | 0.03 (−0.10 to 0.17) |
| +10 | 2.64 (2.26–3.02) | −0.76 | 0.46 (−0.91 to 1.92) | 0.08 (−0.06 to 0.22) |

Under the +10 threshold, 2 P2 samples fell below the *rpoB* QC threshold, leaving 800 samples in the model.

### Taxonomic attribution of *acsB* best hits (post hoc)
The *acsB* reference set includes Wood–Ljungdahl homologues from non-acetogenic lineages: 122 of 371 references are archaeal and 49 come from Thermodesulfobacteriota [`acsb/acsb_ref_taxonomy.tsv`]. Attributing each read to the NCBI taxonomy of its best-hit reference [`acsb/acsb_reads_by_phylum.tsv`]:

- Bacillota received 81–95% of *acsB* reads in every group; 63% of the Bacillota reads hit Lachnospiraceae references [`acsb/acsb_reads_by_family.tsv`].
- Spirochaetota received 4–14%.
- Archaeal references received ≤ 1.5% in each P2 group and 3.6% in P1 foregut samples, mostly Thermoplasmatota.
- Thermodesulfobacteriota received 8.4% in avian-caecal samples and < 1% elsewhere.

We repeated the analyses with *acsB* counted only from bacterial references (228 references) or only from Bacillota references (150). Conclusions were unchanged [`acsb_bacterial/`, `acsb_bacillota/`]:

| *acsB* references | P1 pooled ALR_acet (95% CI) | P2 ΔDIC | P2 foregut − hindgut (95% HPD) | P3 pooled ρp (95% CI) |
|---|---|---|---|---|
| All (primary) | 2.62 (2.23–3.01) | −0.08 | 0.47 (−0.95 to 1.79) | 0.03 (−0.10 to 0.17) |
| Bacterial only | 2.66 (2.28–3.04) | −0.02 | 0.53 (−0.79 to 2.05) | 0.03 (−0.11 to 0.16) |
| Bacillota only | 2.74 (2.36–3.13) | 0.36 | 0.52 (−0.90 to 1.96) | 0.02 (−0.11 to 0.16) |

### Marker specificity benchmark
Across 43 genomes (v2 truth) [`bench/bench_metrics_v2.tsv`]:

| Marker | Carrier genomes | Precision | Sensitivity | False positives per 10⁶ reads in non-carriers |
|---|---|---|---|---|
| *mcrA* | 6 | 0.88 | 0.97 | 0 |
| *acsB* | 11 | 0.84 | 0.96 | 44 |
| *dsrA* | 4 | 0.85 | 0.98 | 0 |
| *asrA* | 10 | 1.00 | 0.77 | 0 |
| *nrfA* | 11 | 1.00 | 0.94 | 0 |
| *narG* | 3 | 0.94 | 0.98 | 4.5 |
| *napA* | 3 | 0.65 | 0.98 | 42 |
| *frdA* | 6 | 0.57 | 0.96 | 16 |
| *rpoB* | 43 | 1.00 | 0.78 | 0 |

The 2 × 10⁶ cattle reads produced no assignment to any marker [`bench/bench_host_decoy_v2.tsv`]. The main error sources were [`bench/bench_assigned_products_v2.tsv`]:

- ***mcrA*:** 407 reads from *mcrG* (the methyl-coenzyme M reductase γ-subunit gene) in methanogen genomes.
- ***acsB*:** 709 of the 712 false reads came from one locus in *Methanomassiliicoccus luminyensis*. It is annotated as a partial pseudogene with a homing-endonuclease domain and sits within the acetyl-CoA decarbonylase/synthase gene cluster, plausibly an intein-interrupted ACS β-subunit gene. Excluding that locus, *acsB* precision was 0.999.
- ***napA* and *narG*:** other molybdopterin oxidoreductases.
- ***frdA*:** succinate dehydrogenase and L-aspartate oxidase. Counting a menaquinol:fumarate reductase flavoprotein annotated *sdhA* as correct raises *frdA* precision to 0.70.

*rpoB* sensitivity was lower because archaeal RNA polymerase genes are not represented in the bacterial normaliser. With the pre-written v1 truth patterns, *mcrA* precision was 0 and *narG*/*napA* precision was 0.23–0.25 [`bench/bench_metrics.tsv`]; the audit table shows these as annotation-wording misses.

Sensitivity was high for both primary markers on this panel, with precision of 0.88 for *mcrA* and 0.84 for *acsB* under the revised truth definition. Cross-assignment for *napA* and *frdA* can inflate the respiratory sum used in ALR_resp. Error rates in the benchmark apply to this genome panel and its annotation rules.

### Dated-tree sensitivity
TimeTree dated 152 of 165 nodes; 13 were interpolated and 22 were adjusted for conflicts with the topology (root 435 Myr) [`dated_tree/node_ages.tsv`]. With the dated tree, fermentation site still did not improve fit: ΔDIC = 0.99; foregut − hindgut 0.70 (95% HPD −0.96 to 2.29); hindgut − simple −1.14 (−3.11 to 0.48) [`dated_tree/p2_primary.tsv`]. The phylogenetic variance share was higher, 0.78 (0.61–0.87), and the species share lower [`dated_tree/p2_variance_partition.tsv`].

### Reproducibility
The 17 P1 runs quantified before DIAMOND hits were archived were processed again from scratch. The second pass used a different read source (NCBI SRA on AWS instead of ENA files) and a different DIAMOND block size. All 17 marker tables were byte-identical to the originals [`p1_rerun_reproducibility.txt`]. For the archived analysis, the complete statistical analysis was re-run from the per-sample outputs with a single script (`run_all.sh`; 57 min on 16 cores). All 103 result tables were byte-identical to those reported here [`rerun_comparison.tsv`].

## Discussion

### Paired compartments in ruminants

The two public ruminant datasets showed a higher *acsB*:*mcrA* marker ratio in every hindgut comparison. Li et al. had previously reported caecal enrichment of *acsB* in these data sources [L02, L04]. The paired and pooled estimates here quantify the consistency and magnitude of that pattern under a shared pipeline. The hindgut difference was driven mainly by *acsB*. The *mcrA* interval included both lower and higher values, limiting inference about methanogenic potential on its own.

Most hindgut *acsB* reads matched Bacillota references, and the paired difference persisted when counts were restricted to Bacillota or bacterial references. These assignments are consistent with enrichment of bacterial acetyl-CoA synthase genes. A genome-resolved rumen survey used *acsB* alongside other Wood–Ljungdahl genes to identify putative acetogens [L17]. Best-hit taxonomy and gene counts cannot establish reductive acetogenesis in the source organisms or distinguish all uses of the Wood–Ljungdahl pathway. Several respiratory markers also differed between chambers. Their variable specificity and the lack of activity measurements limit conclusions about which pathways consumed H₂.

### Fermentation site and faecal markers across hosts

The fitted P2 models gave little evidence that fermentation site accounted for faecal log-ratio variation after adjustment for phylogeny, study and measured covariates. On the two main tree specifications, the intervals for the foregut-minus-hindgut contrast overlapped zero and lay above −2.62 ln units, the direction and magnitude that would match the paired P1 hindgut advantage. Smaller cross-host differences remain compatible with the estimates. This comparison is descriptive because P1 contrasts chambers within ruminants and P2 compares distal-gut material across hosts.

Faecal samples mainly represent the distal gut. Ruminant rectal samples in P1 resembled the caecal samples in their enrichment relative to rumen samples. Consequently, P2 can only capture anatomical differences that are reflected in faeces. Foregut and hindgut hosts nevertheless had higher *rpoB*-normalised *mcrA* and *acsB* abundances than simple-gut hosts. At low counts, the ratio is sensitive to pseudocount choice, particularly in the simple-gut group, which had few *mcrA* detections.

Estimated phylogenetic variance was substantial under both tree specifications, and the difference between 39% and 78% shows its sensitivity to branch lengths and variance allocation. The archived model summaries provide no multi-chain convergence assessment. Site categories were assigned by lineage, so the model has limited leverage to distinguish anatomy from other lineage-associated host traits. The within-study fit reduces differences between contributing studies but covers only 12 studies and does not remove differences between hosts within those studies.

### Within-study association and respiratory markers

The mean proportionality estimate was near zero, with large differences among study units. Centred log-ratio coordinates depend on the marker panel and its geometric mean. Sample-level resampling does not account for repeated observations from the same animal. The results provide little evidence for a consistent marker trade-off on this relative scale. They cannot establish whether methanogens and acetogens compete for H₂ or how the pathways would respond to an intervention.

Higher *mcrA* abundance in fermenting hosts contributes to the lower respiratory-to-methanogenic ratio. The ratio does not estimate inverse covariation, and the benchmark found cross-assignment of homologues to *napA* and *frdA*. Direct measurements of expression or H₂ flux and compartment-resolved sampling across vertebrate lineages would address questions that remain open.

### Limitations

P1 draws on two studies of ruminants. P2 is a capped, observational sample of public faecal metagenomes. Diet was grouped into broad classes and missing in 145 samples, and captivity was unknown in 265 samples. Most inventory studies sampled a single fermentation category. The initial screen admitted 20 mouse samples, including samples from a colitis model, and four libraries whose individual reads were shorter than the threshold. One analyst managed the automated screening and metadata curation. The no-mouse sensitivity analysis was performed after the primary analysis.

All abundance estimates depend on five million forward reads per run before quality control, the marker database, sequence-assignment thresholds and the bacterial *rpoB* denominator. P1 used head sampling, and the limited comparison with multi-offset sampling found a difference for *acsB* in one run. The planned phylogenetic placement of *acsB* could not be completed because the original reads were not retained. The post hoc best-hit attribution and specificity benchmark have limited taxonomic resolution and cannot establish error rates in mixed gut communities. Marker abundances describe DNA recovered from the sampled communities. Neither per-study risk of bias nor certainty of evidence was formally graded.

### Conclusion

The paired ruminant samples consistently showed hindgut enrichment of *acsB* relative to *mcrA*. Across the sampled vertebrates, fermentation site added little explanatory value for the faecal marker ratio under the fitted models. These results apply to the sampled compartments and faecal communities and describe marker-gene abundance.

## Declarations
**Funding.** The authors funded the study.

**Competing interests.** The authors declare no competing interests.

**Author contributions (CRediT).** Ananda Prasad: Conceptualization, Methodology (overall research design), Writing – original draft. Angsuman Das: Formal analysis, Investigation. Both authors reviewed and approved the preceding manuscript draft.

**Acknowledgements.** We thank the authors of the original studies for depositing their sequencing data and metadata in public archives, the developers of the Greening-lab marker gene databases, and the maintainers of ENA/SRA, Open Tree of Life and TimeTree.

**Use of generative language models.** Under the authors’ direction, a generative language model helped write and run analysis code, execute analyses, and draft and edit manuscript text, including the Discussion. A generative language model also assisted with the editorial and methodological review of this revision. The authors defined the research question and study design. They reviewed the code, results and text in the archived version and take responsibility for the work.

## Figure legends
**Figure 1. Study design.** (a) Predicted differences in marker abundance between foregut and hindgut chambers and relative marker associations for P3 and P4. (b) Design: P1, a within-animal foregut-versus-hindgut contrast in ruminants used as a positive control; P2, a cross-vertebrate phylogenetic mixed model on faecal samples; P3, within-study proportionality, with the robustness checks performed. (c) Processing pipeline.

**Figure 2. Within-animal foregut-versus-hindgut contrast (P1).** (a) Pooled hindgut − foregut differences in ln units for the primary outcome (ALR_acet, the pseudocount-adjusted *acsB*:*mcrA* log ratio) and secondary outcomes; points are REML estimates with 95% Hartung–Knapp confidence intervals (black) and 95% prediction intervals (grey). (b) ALR_acet difference for each study × host-species stratum (mean and approximate 95% CI; *n* = animals); the grey band is the pooled 95% CI. Source: `p1_pooled.tsv`, `p1_stratum_effects.tsv`.

**Figure 3. Faecal ALR_acet across vertebrates (P2).** (a) Species-mean ALR_acet by fermentation site; point area increases with the number of samples (maximum 10) and horizontal bars are medians. (b) Posterior differences between fermentation sites (mean and 95% HPD interval) from the phylogenetic mixed model with phylogeny, species and study random effects and diet, captivity, depth and read length as fixed effects. Source: `p2_species_means.tsv`, `p2_contrasts.tsv`.

**Figure 4. Proportionality between *mcrA* and *acsB* within studies (P3).** Proportionality ρp with 95% intervals derived from the bootstrap variance on the Fisher z scale for each P2 study with ≥ 10 samples (grey) and each P1 study × chamber unit (blue); point size reflects sample size. The vertical line and grey band show the pooled estimate and its 95% CI. Negative ρp indicates negative association between the centred log-ratio coordinates. Source: `p3_units.tsv`, `p3_pooled.tsv`.

**Figure S1. Marker detection by fermentation site.** Share of P2 samples with ≥ 2 reads for each marker gene. Source: `detection_by_site.tsv`.

**Figure S2. PRISMA flow diagram.** Screening of ENA runs to the P2 analysis set, with exclusions at each step; the P1 paired-compartment set is summarised below. Source: `prisma_log.tsv`, `prisma_final.tsv`.

## References

- [L01] Bowerman KL, et al. (2026) *mSystems* 11:e01608-25. doi:10.1128/msystems.01608-25
- [L02] Li Q, et al. (2025) *Microbiome* 13:28. doi:10.1186/s40168-024-02018-1
- [L03] Greening C, et al. (2019) *ISME J* 13:2617–2632. doi:10.1038/s41396-019-0464-2
- [L04] Xie F, et al. (2021) *Microbiome* 9:137. doi:10.1186/s40168-021-01078-x
- [L05] Youngblut ND, et al. (2020) *mSystems* 5:e01045-20. doi:10.1128/mSystems.01045-20
- [L07] Youngblut ND, et al. (2021) *Nat Microbiol* 6:1443–1454. doi:10.1038/s41564-021-00980-2
- [L08] Levin D, et al. (2021) *Science* 372:eabb5352. doi:10.1126/science.abb5352
- [L10] Kuhn M, et al. (2026) *Nucleic Acids Res* 54:D826–D834. doi:10.1093/nar/gkaf1118
- [L17] Li Q, et al. (2025) *ISME J* 19:wraf183. doi:10.1093/ismejo/wraf183
- [L20] Stevens CE, Hume ID (1998) *Physiol Rev* 78:393–427. doi:10.1152/physrev.1998.78.2.393
- [L21] Wilman H, et al. (2014) *Ecology* 95:2027. doi:10.1890/13-1917.1
- [L22] Leung PM, Greening C (2021) Greening lab metabolic marker gene databases. Monash Bridges. doi:10.26180/c.5230745
- [L23] Ricci F, et al. (2025) *Nat Commun* 16:11634. doi:10.1038/s41467-025-66734-4
- [L25] UniProt Consortium (2025) *Nucleic Acids Res* 53:D609–D617. doi:10.1093/nar/gkae1010
- [L26] Buchfink B, et al. (2021) *Nat Methods* 18:366–368. doi:10.1038/s41592-021-01101-x
- [L27] Chen S, et al. (2018) *Bioinformatics* 34:i884–i890. doi:10.1093/bioinformatics/bty560
- [L28] Hadfield JD (2010) *J Stat Softw* 33(2). doi:10.18637/jss.v033.i02
- [L29] Viechtbauer W (2010) *J Stat Softw* 36(3). doi:10.18637/jss.v036.i03
- [L30] IntHout J, et al. (2014) *BMC Med Res Methodol* 14:25. doi:10.1186/1471-2288-14-25
- [L31] Page MJ, et al. (2021) *J Clin Epidemiol* 134:178–189. doi:10.1016/j.jclinepi.2021.03.001
- [L32] Mirzayi C, et al. (2021) *Nat Med* 27:1885–1892. doi:10.1038/s41591-021-01552-x
- [L33] Erb I, Notredame C (2016) *Theory Biosci* 135:21–36. doi:10.1007/s12064-015-0220-8
- [L34] Hinchliff CE, et al. (2015) *Proc Natl Acad Sci USA* 112:12764–12769. doi:10.1073/pnas.1423041112
- [L35] Grafen A (1989) *Philos Trans R Soc B* 326:119–157. doi:10.1098/rstb.1989.0106
- [L36] Kumar S, et al. (2022) *Mol Biol Evol* 39:msac174. doi:10.1093/molbev/msac174
- [L37] O'Leary NA, et al. (2024) *Sci Data* 11:732. doi:10.1038/s41597-024-03571-y
