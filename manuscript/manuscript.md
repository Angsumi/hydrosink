# [TITLE: set after P2 results] Gut fermentation architecture and hydrogen-sink potential across vertebrates: a pre-registered, phylogenetically controlled synthesis of public metagenomes

*Draft v1 (2026-11-17). Every number in this draft comes from a script-generated file in `analysis/results/` (file named in brackets). Items marked **[NOT DONE]** were pre-specified but not executed; they are listed under Limitations. Reference codes [Lxx] map to `literature.csv`; only entries marked VERIFIED there are cited.*

---

## Abstract
[RESULTS-DEPENDENT; written last]

## Introduction
Fermentation in vertebrate guts releases molecular hydrogen (H₂). If H₂ accumulates, re-oxidation of reduced cofactors becomes thermodynamically unfavourable and fermentation slows, so hydrogenotrophic microbes that keep H₂ partial pressure low are central to gut function [L03]. Several guilds compete for this H₂. Hydrogenotrophic methanogens reduce CO₂ to methane, an energy loss to the host and a greenhouse gas. Reductive acetogens produce acetate, which the host can absorb. Respiratory hydrogenotrophs couple H₂ oxidation to the reduction of sulfate, nitrate, nitrite or fumarate [L03]. In the rumen, H₂ uptake is distributed among these guilds, and alternative uptake pathways are expressed more in low-methane sheep [L03].

Which sink dominates appears to depend on where fermentation happens. In 16 goats, the reductive-acetogenesis marker *acsB* was about 12-fold more abundant in the caecum than in the rumen, and the authors reported the same direction in several ruminant species [L02]. Marsupials, in contrast, are heterogeneous: some species resemble high-methane mammals while others are enriched for respiratory H₂-uptake potential [L01]. Across 110 vertebrate species, archaeal diversity tracks host phylogeny more than diet, but that analysis used 16S amplicons and did not measure H₂-sink genes [L07]. Large shotgun catalogues of vertebrate guts now exist [L04, L05, L08], as do curated cross-study metadata [L10]. To our knowledge, no study has tested with phylogenetic and study-level control whether gut fermentation architecture predicts how H₂-disposal potential is partitioned across vertebrates.

We tested four predictions that were fixed, with their analysis code, before the marker data were inspected (git tags `prereg-v1`, `prereg-v1.1`, `prereg-v1.2`):
- **P1** (within animals; positive control): in paired foregut and hindgut samples from the same animal, ln(*acsB*/*mcrA*) is higher in the hindgut.
- **P2** (across hosts): fermentation site (foregut, hindgut, avian caecal, simple gut) explains variation in faecal ln(*acsB*/*mcrA*) beyond host phylogeny, diet, captivity and technical covariates.
- **P3** (trade-off): within studies, methanogenic and acetogenic potential are negatively proportional.
- **P4** (exploratory): respiratory H₂-uptake potential is higher where methanogenic potential is lower.

P1 replicates the published compartment result [L02] with the same and additional ruminant data [L04]. It serves to calibrate the marker pipeline and is not presented as a new finding. The novelty claim rests on P2–P4.

## Methods
### Protocol, registration and amendments
Hypotheses, inclusion criteria, marker thresholds, outcomes and models were fixed in `analysis/config/config.yaml`, `research_plan.md` and `analysis/scripts/05_meta_analysis.R`, and committed (git tag `prereg-v1`) before any marker data from the main sets were inspected. Two amendments were committed before any P2 marker data existed:
- **v1.1** (`analysis/PREREG_AMENDMENT_v1.1.md`): P1 re-labelled a positive control; per-study P1 estimates and a three-level model added; within-study, known-captivity, sequencing-platform and marsupial hold-out analyses added; outcome-blind curation of captivity metadata.
- **v1.2** (`analysis/PREREG_AMENDMENT_v1.2.md`): change of P2 read sampling after a pre-specified gate failed (below).

The study was not registered on an external registry; the time-stamped git history is the record.

### Data sources and screening
The European Nucleotide Archive (ENA) Portal API was queried on 15 November 2026 for Illumina WGS runs with library source METAGENOMIC and one of seven gut-associated metagenome taxa: gut, feces, rumen, bovine gut, pig gut, chicken gut and fish gut metagenome. Screening is shown as a PRISMA 2020 flow [L31] [`analysis/prisma_log.tsv`, `prisma_final.tsv`]. Criteria:
- **Host:** a structured host taxon ID within Vertebrata, excluding human, mouse and rat.
- **Keywords:** no keywords indicating treatments, inhibitors, in-vitro systems, controls or non-WGS libraries.
- **Read count:** ≥ 5 × 10⁶ reads.
- **Read length:** mean ≥ 100 bp.
- **One run per BioSample:** the deepest.

For P2 we kept faecal samples and capped each host species at 10 runs, drawn round-robin across studies with a fixed seed (20261115). For P1 we used the two studies with paired fermentation-chamber samples from the same animals: Li et al. 2025 (PRJNA1068803; 16 goats, rumen and caecum) [L02] and Xie et al. 2021 (PRJNA657455; seven ruminant species) [L04]. P1 was restricted to fermentation chambers: rumen, reticulum and omasum (foregut) versus caecum, colon and rectum (hindgut).

### Host traits (assigned before marker data were inspected)
**Fermentation site** was assigned from host lineage by fixed rules following Stevens & Hume [L20] (`analysis/scripts/01b_host_traits.py`):
- **Foregut:** Ruminantia, Camelidae, Hippopotamidae, Tayassuidae, Macropodidae and Potoroidae.
- **Hindgut:** Equidae, Rhinocerotidae, Tapiridae, Elephantidae, Sirenia, Lagomorpha, Rodentia, wombats, koala and phalangeriform marsupials, Hyracoidea, Suidae, non-colobine primates and herbivorous Testudines.
- **Avian caecal:** Galliformes, Anseriformes, ratites, Tinamiformes, Gruiformes and Otidiformes.
- **Simple:** all other vertebrates.

**Diet class** (herbivore, omnivore or carnivore) came from EltonTraits 1.0 [L21]; hosts without an EltonTraits entry were coded "unknown".

**Captivity** (wild, captive, domestic, laboratory or unknown) was curated from study and BioSample text by keyword rules plus 38 manual study-level decisions, each with an evidence quote (`analysis/config/captivity_manual.tsv`). This was done blind to marker data.

**Host phylogeny** was the Open Tree of Life synthetic tree induced for the matched taxa (`01c_host_tree.py`). Because that tree is topology-only, Grafen branch lengths were added (ape `compute.brlen`). Hosts identified only to genus without a tree tip were excluded from the phylogenetic models.

### Read sampling
Bandwidth precluded full downloads, so each run was represented by 5 × 10⁶ forward (R1) reads; R1 only was used, as in the Greening-lab protocol [L23]. **P1** used the first 5 × 10⁶ reads ("head" sampling).

A pre-specified gate (Gate 1b) compared head sampling against a uniform random subsample of the full file in two goat samples. The rule required every gene group to agree within 10% when ≥ 100 reads were available. One of ten checks failed: [NiFe]-hydrogenase in one rumen sample, 210 versus 189 reads, +10.3% [`gate1b_head_vs_random.tsv`]. As pre-specified, P2 therefore switched to **multi-offset sampling** (amendment v1.2): 10 blocks of 5 × 10⁵ consecutive spots, with block *k* starting at spot ⌊*k*·*S*/10⌋ + 1, where *S* is the run's spot count.

Reads came from the NCBI SRA copy on AWS. A run was accepted only if every block returned exactly the planned rows. Three fetch implementations were used: remote `vdb-dump`, local extraction from the downloaded `.sra` object, and `fastq-dump`. They were verified to give identical bases and qualities (amendment v1.2 addenda).

### Quality control and marker quantification
- **Quality control:** fastp 0.23.4 [L27] (`-q 20 -l 100 --trim_poly_g --cut_tail`).
- **Search:** DIAMOND 2.1.10 blastx [L26] (`--fast --id 40 --query-cover 80 --evalue 1e-5 --max-target-seqs 1 --max-hsps 1`).
- **Reference:** the Greening-lab compiled metabolic marker database v3 (51 genes, all retained as competitors) [L22]. To this were added 5,635 full-length UniRef90 bacterial RpoB sequences (UniProt 2026_04) [L25] as the single-copy normaliser, and RpoC sequences as decoys (170,171 sequences in total; md5 values in `config.yaml`).
- **Assignment:** a read was assigned to the gene of its best hit if percent identity met that gene's published read-level threshold [L23]. This was 50% for most genes, including *mcrA*, *acsB*, *dsrA*, *nrfA*, *narG*, *napA* and *frdA*, and 60% for [FeFe]- and group-4 [NiFe]-hydrogenases.
- **Abundance:** RPKM = reads / (mean reference length in kb) / (QC reads in millions). Copies per genome (cpg) = RPKM(gene) / RPKM(*rpoB*) (`aggregate_markers.py`).
- **Sample QC:** a sample passed if it had ≥ 200 *rpoB* reads.

The DIAMOND block size was reduced from 1 to 0.4 during P2 to bound memory. A test run gave identical hit files (amendment v1.2 addendum). Processing used one workstation and four Google Cloud VMs with identical tool versions. A run processed on both gave byte-identical marker tables.

### Outcomes
- **Primary:** ALR_acet = ln[(*acsB* + *c*)/(*mcrA* + *c*)], where *c* is half the minimum non-zero cpg of each marker in the analysed set.
- **Secondary:**
  - ALR_resp, the same ratio with the summed respiratory reductases (*dsrA*, *asrA*, *nrfA*, *narG*, *napA*, *frdA*) in the numerator;
  - log cpg of each marker;
  - group-1 [NiFe]- and group-B [FeFe]-hydrogenases;
  - *mcrA* detection (≥ 2 reads).

### Statistical analysis
All analyses ran in R 4.3.3 with metafor 4.8-0 [L29], lme4 1.1-37 and MCMCglmm 2.36 [L28] (`05_meta_analysis.R`; seed 20261115).

**P1.** For each animal we computed the hindgut-minus-foregut difference of chamber means. Differences were averaged within each study × host-species stratum, and strata were pooled by random-effects meta-analysis (REML τ², Hartung–Knapp–Sidik–Jonkman CIs [L30], 95% prediction intervals, leave-one-out). Because P1 has only two source studies, we also report each study separately and a three-level model (strata nested in study). A one-stage linear mixed model was fitted as a cross-check: random animal intercepts and random chamber slopes by stratum. Secondary outcomes were adjusted by Benjamini–Hochberg.

**P2.** We fitted Gaussian MCMCglmm models with random effects for phylogeny (inverse relatedness matrix), host species and study. Fixed effects were diet class, captivity, log₁₀ *rpoB* reads and read length, plus fermentation site in the full model. Chains ran for 130,000 iterations, with 30,000 burn-in and thinning of 100. Priors were parameter-expanded for the random effects (V = 1, ν = 1, α.V = 1000) and inverse-Wishart for the residual (V = 1, ν = 0.002). Models with and without fermentation site were compared by DIC. We report posterior contrasts between sites (mean, 95% HPD interval) and the phylogenetic share of random-effect plus residual variance. *mcrA* detection was modelled with the same structure (categorical family; hurdle component).

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

**P3.** For each P2 study with ≥ 10 samples, and each P1 study × chamber, we computed proportionality ρp between *mcrA* and *acsB* on CLR-transformed marker profiles [rpoB and 10 markers]. Values were Fisher-z transformed with bootstrap variances (1,000 resamples) and pooled by REML with Hartung–Knapp. We also report a subset restricted to units with *mcrA* detected in ≥ 50% of samples.

**Sensitivity of thresholds.** Marker counts were re-derived from the saved DIAMOND hits with all identity thresholds shifted by −10 and +10 percentage points (floor 40%, the search cut-off), and the P1, P3 and primary P2 analyses were repeated (`05c_identity_sensitivity.py`, `05d_sensitivity_sets.R`). With a shift of 0 the re-derivation reproduced all 1,163 original marker tables byte for byte.

**P1 sampling check.** Five P1 runs drawn with a fixed seed before the data were analysed were re-quantified with multi-offset sampling and compared with head sampling (amendment v1.2 §2).

### Data and code availability
Run accessions, curated metadata and per-sample marker tables are in `analysis/inventory/` and `analysis/results/per_sample/`; DIAMOND hits are in `analysis/results/hits/`, with SHA-256 manifests (`p1_manifest.sha256`, `p2_manifest.sha256`). All code is in `analysis/scripts/`. [Archive with DOI at submission.]

## Results
[RESULTS-DEPENDENT]

## Discussion
[RESULTS-DEPENDENT]

## References
[assembled from literature.csv, VERIFIED entries only]
