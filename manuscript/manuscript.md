# Hindgut enrichment of acetogenic potential within ruminants does not generalise across vertebrates: a pre-registered, phylogenetically controlled synthesis of gut hydrogen-sink genes in public metagenomes

*Draft v1 (2026-11-17). Every number in this draft comes from a script-generated file in `analysis/results/` (file named in brackets). Items marked **[NOT DONE]** were pre-specified but not executed; they are listed under Limitations. Reference codes [Lxx] map to `literature.csv`; only entries marked VERIFIED there are cited.*

---

## Abstract
**Background.** Gut microbes dispose of fermentative H₂ mainly through methanogenesis, reductive acetogenesis and respiratory pathways. Within ruminants, the hindgut is enriched for acetogens relative to the rumen, but whether gut fermentation architecture governs H₂-sink partitioning across vertebrates is untested.
**Methods.** We fixed four predictions, outcomes and analysis code before inspecting data. Public shotgun gut metagenomes were screened from ENA, and 5 × 10⁶ forward reads per run were searched against a curated marker database. Abundances were expressed as copies per genome (*rpoB*-normalised). The primary outcome was ln(*acsB*/*mcrA*). Within-animal foregut-versus-hindgut contrasts (53 ruminants, two studies) were pooled by random-effects meta-analysis. Cross-host faecal data (802 samples, 166 species, 179 studies) were analysed with Bayesian phylogenetic mixed models.
**Results.** In all 53 animals the hindgut had a higher *acsB*:*mcrA* ratio than the foregut: pooled 2.62 ln units (95% CI 2.23–3.01), about 12-fold, driven by *acsB*. Across vertebrates, fermentation site did not predict the faecal ratio (ΔDIC −0.08; foregut − hindgut 0.47, 95% HPD −0.95 to 1.79), in every pre-specified sensitivity analysis. Host phylogeny accounted for 39% of variance (95% HPD 21–67%). Foregut and hindgut hosts carried 21–53-fold more of both *mcrA* and *acsB* than simple guts. Methanogenic and acetogenic potential were not negatively proportional within studies (pooled ρp 0.03, 95% CI −0.10 to 0.17).
**Conclusions.** The hindgut acetogenic shift is a robust within-ruminant regularity, but it does not translate into a cross-vertebrate faecal pattern. Gut architecture predicts the amount of hydrogenotrophic potential, while its partitioning follows host phylogeny. These results constrain generalisations from ruminant compartments to vertebrate methane biology.

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

***acsB* attribution (post-hoc).** The pre-planned phylogenetic placement of *acsB* reads was not possible because reads were not retained. Instead, each *acsB* read was attributed to the NCBI taxonomy (E-utilities) of its best-hit reference protein. P1, P3 and the primary P2 model were repeated with *acsB* restricted to bacterial or to Bacillota references (`08a`–`08c` scripts).

**P1 sampling check.** Five P1 runs drawn with a fixed seed before the data were analysed were re-quantified with multi-offset sampling and compared with head sampling (amendment v1.2 §2).

### Data and code availability
Run accessions, curated metadata and per-sample marker tables are in `analysis/inventory/` and `analysis/results/per_sample/`; DIAMOND hits are in `analysis/results/hits/`, with SHA-256 manifests (`p1_manifest.sha256`, `p2_manifest.sha256`). All code is in `analysis/scripts/`. [Archive with DOI at submission.]

## Results
### Data set
Screening reduced 136,708 candidate runs to 19,385 eligible runs, of which 8,293 were faecal. Capping at 10 runs per host species gave 930 P2 runs [`prisma_log.tsv`].

Of the 930 runs:
- 926 were quantified. The other 4 runs (one study) consisted entirely of 76-bp reads, below the length filter, and yielded no reads after QC. The screen's read-length estimate had counted both mates of each spot and so missed this.
- 845 passed the *rpoB* QC threshold.
- 802 samples from 166 host species and 179 studies entered the phylogenetic models [`prisma_final.tsv`]. Genus-only hosts without a tree tip were excluded at this step.

The modelled set comprised 207 foregut, 250 hindgut, 67 avian-caecal and 278 simple-gut samples; 655 were mammalian [`table1_p2.tsv`]. Captivity was wild for 274 samples, captive, domestic or laboratory for 263, and unknown for 265.

Read depth was uniform by design: median 4.96 × 10⁶ reads after QC, median 5,707 *rpoB* reads per sample. *mcrA* (≥ 2 reads) was detected in 88% of foregut, 60% of hindgut, 28% of avian-caecal and 12% of simple-gut samples; *acsB* in 97%, 98%, 78% and 67% [`detection_by_site.tsv`; Fig. S1].

P1 comprised 253 of 254 runs passing QC, from 53 animals with both chambers sampled: 16 goats from Li 2025 [L02] and 37 animals of seven ruminant species from Xie 2021 [L04] [`table1_p1.tsv`].

### P1: the hindgut is enriched for acetogenic relative to methanogenic potential in every animal
ALR_acet was higher in the hindgut than in the foregut in **53 of 53 animals**. The pooled stratum-level difference was **2.62 ln units** (95% CI 2.23–3.01; 95% prediction interval 1.69–3.55; *k* = 8 strata; I² = 70%) [`p1_pooled.tsv`; Fig. 2]. This corresponds to a median per-animal 12.5-fold higher *acsB*:*mcrA* ratio [`p1_per_animal.tsv`].

The effect was present in each study separately: Li 2025, 1.99 (1.54–2.44); Xie 2021, 2.71 (2.33–3.08). It was also present in the three-level model with strata nested in study: 2.43 (1.62–3.24) [`p1_per_study.tsv`]. Other checks:
- One-stage mixed model: 2.72 (2.39–3.06) [`p1_mixed_model.tsv`].
- Leave-one-stratum-out estimates: 2.56 to 2.71 [`p1_leave1out.tsv`].
- Rumen versus caecum only: 2.72. Rumen versus rectum: 2.80 [`p1_sensitivity.tsv`].

The contrast was driven by *acsB*: +2.92 ln units (2.15–3.70), higher in 53/53 animals. *mcrA* did not differ: +0.25 (−0.31 to 0.80; *q* = 0.39). Among secondary markers [`p1_pooled.tsv`]:
- *dsrA* and group-B [FeFe]-hydrogenase were higher in the hindgut.
- *nrfA* (−2.92), *narG* (−1.00), group-1 [NiFe]-hydrogenase (−0.94) and the respiratory sum (−1.33) were lower.
- ALR_resp was therefore lower in the hindgut: −1.58 (−2.53 to −0.62).

Thus, within ruminants, the hindgut has more acetogenic but less respiratory H₂-uptake potential than the foregut, with similar methanogenic potential.

The pre-specified head-versus-multi-offset sampling check was run on 5 P1 runs [`p1_offsets_vs_head.tsv`]. In 22 of 25 run × gene-group comparisons the 95% CI of the ratio included 1. In one rumen sample (SRR12529313), head sampling gave 32% more *acsB* than multi-offset sampling (ratio 0.68, 95% CI 0.62–0.74). That is 0.39 ln units, small relative to the P1 contrast but larger than counting noise.

### P2: fermentation site does not predict faecal ALR_acet across vertebrates
In the phylogenetic mixed model, adding fermentation site did not improve fit for ALR_acet: ΔDIC = −0.08 [`p2_model_comparison.tsv`]. No site contrast had a 95% HPD interval excluding zero [`p2_contrasts.tsv`; Fig. 3]:

| Contrast | Posterior mean | 95% HPD |
|---|---|---|
| Foregut − hindgut | 0.47 | −0.95 to 1.79 |
| Hindgut − simple | −0.66 | −2.19 to 1.01 |
| Foregut − simple | −0.19 | −2.03 to 1.76 |
| Avian caecal − simple | 1.13 | −0.88 to 3.17 |

Of the variance in ALR_acet, the share attributable to host phylogeny was 0.39 (95% HPD 0.21–0.67), to species 0.15 and to study 0.16 [`p2_variance_partition.tsv`]. Diet, captivity, depth and read length had no credible effects [`p2_fixed_effects.tsv`].

The null result for the site contrasts was consistent across every pre-specified sensitivity analysis; all contrast intervals included zero [`p2_sensitivity.tsv`, `p2_sensitivity_v1_1.tsv`]:
- mammals only (626 samples);
- non-domestic (691);
- ≤ 3 samples per species (355);
- ≥ 1,000 *rpoB* reads (747);
- herbivores only (452);
- known captivity (537);
- with a sequencing-platform covariate.

Restricting to the 12 studies that sampled more than one fermentation site (136 samples, 72 species; study as a fixed effect) also gave no clear site effect: foregut −1.25 (SE 0.89), hindgut −0.96 (SE 0.76) relative to simple [`p2_within_study.tsv`]. These within-study estimates point in the direction opposite to P1.

In the marsupial hold-out (19 samples, 4 species), the model predicted no hindgut–foregut difference (+0.22). The observed difference was −3.38, but it rests on a single foregut sample. The species-level Spearman correlation was 0.4 on 4 species [`p2_holdout_marsupial.tsv`]; this validation is uninformative.

Fermentation site did predict the absolute abundance of both markers, in the same direction [`p2_contrasts.tsv`]:

| Contrast | ln *acsB* (95% HPD) | ln *mcrA* (95% HPD) |
|---|---|---|
| Hindgut − simple | +3.04 (1.69–4.56) | +3.58 (2.44–4.68) |
| Foregut − simple | +3.91 (2.26–5.36) | +3.97 (2.45–5.32) |

Avian-caecal hosts had more *acsB* than simple-gut hosts (+2.08; 0.36–3.62), but *mcrA* did not differ credibly. In the hurdle model, *mcrA* detection was strongly associated with all three fermenting categories relative to simple guts [`p2_fixed_effects.tsv`]. Fermenting guts therefore carry more of both sinks, and ALR_acet does not separate them.

Because *mcrA* was undetected in 88% of simple-gut samples, ALR_acet in simple guts is often bounded by the pseudocount. 115 samples had neither marker detected, 90 of them from simple guts [`p2_analysis_set.tsv`]. The log-abundance and hurdle models are therefore the more informative description of that category.

**Post-hoc check (not pre-registered).** The ENA screen excluded house mouse by taxon 10090 only, so 20 *Mus musculus* samples entered P2 under a subspecies or genus taxon, 10 of them from a DSS-colitis model. In addition, one *Elephas maximus* sample lacked a diet label. Removing all *Mus musculus* and harmonising that label gave the same result: 782 samples, ΔDIC = −0.18, foregut − hindgut 0.45 (−0.98 to 1.94) [`posthoc_no_mouse/p2_primary.tsv`].

### P3: no methanogen–acetogen trade-off
Across 34 units (30 P2 studies with ≥ 10 samples and 4 P1 study × chamber units), the pooled proportionality between *mcrA* and *acsB* was ρp = 0.03 (95% CI −0.10 to 0.17; 95% prediction interval −0.56 to 0.61; I² = 75%) [`p3_pooled.tsv`; Fig. 4]. In units where *mcrA* was detected in ≥ 50% of samples (*k* = 18), it was −0.01 (−0.20 to 0.18). Unit estimates ranged from −0.80 to 0.72 [`p3_units.tsv`]. P3 is not supported: the data are consistent with no consistent proportional relationship and do not show a trade-off.

### P4 (exploratory): respiratory potential relative to methanogenic potential
Across hosts, ALR_resp was much lower in hindgut (−5.01; −6.50 to −3.77) and foregut (−5.36; −6.92 to −3.69) than in simple-gut hosts [`p2_contrasts.tsv`]. Comparing with the *mcrA* contrasts above, most of this shift reflects higher *mcrA* in fermenting guts, not lower respiratory potential.

Within animals (P1), the hindgut had lower ALR_resp than the foregut at similar *mcrA* (above). Thus the cross-host pattern is directionally consistent with P4, but within animals respiratory and methanogenic potential did not vary inversely.

### Sensitivity to marker identity thresholds
Shifting every per-gene identity threshold by −10 or +10 percentage points left all conclusions unchanged [`sens_id_m10/`, `sens_id_p10/`]:

| Threshold shift | P1 pooled ALR_acet (95% CI) | P2 ΔDIC | P2 foregut − hindgut (95% HPD) | P3 pooled ρp (95% CI) |
|---|---|---|---|---|
| −10 | 2.61 (2.22–3.00) | −0.17 | 0.48 (−0.86 to 2.02) | 0.01 (−0.12 to 0.15) |
| 0 (primary) | 2.62 (2.23–3.01) | −0.08 | 0.47 (−0.95 to 1.79) | 0.03 (−0.10 to 0.17) |
| +10 | 2.64 (2.26–3.02) | −0.76 | 0.46 (−0.91 to 1.92) | 0.08 (−0.06 to 0.22) |

Under the +10 threshold, 2 P2 samples fell below the *rpoB* QC threshold, leaving 800 samples in the model.

### Taxonomic origin of *acsB* reads (post-hoc)
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

### Reproducibility
The 17 P1 runs quantified before DIAMOND hits were archived were processed again from scratch. The second pass used a different read source (NCBI SRA on AWS instead of ENA files) and a different DIAMOND block size. All 17 marker tables were byte-identical to the originals [`p1_rerun_reproducibility.txt`].

## Discussion
**Summary.** The pre-registered positive control reproduced the published within-ruminant result [L02]. In every one of 53 animals across eight ruminant strata, the hindgut carried more acetogenic relative to methanogenic potential than the foregut, by about 12-fold. The test of generality failed. Across 166 vertebrate species, fermentation site did not predict the faecal *acsB*:*mcrA* ratio once host phylogeny, species and study were modelled, and this held in every pre-specified sensitivity analysis. Instead, fermenting guts carried more of both markers. Host phylogeny accounted for the largest share of variance in the ratio, as previously found for archaeal diversity with 16S data [L07]. Methanogenic and acetogenic potential were not negatively proportional within studies.

**Interpretation.** The research framework anticipated this outcome: a compartment effect inside animals with phylogenetic structure dominating between them (`research_plan.md` §1.5). Three non-exclusive explanations fit the data.
1. **Sample type.** P2 used faeces, which integrate the whole tract. In a ruminant, the faecal community reflects the hindgut and not the rumen, so "foregut" hosts are not observed at their foregut in P2. The P1 rumen-versus-rectum contrast (2.80 ln units) shows that this alone does not abolish within-animal differences. The cross-host faecal comparison, however, contrasts hosts, not chambers.
2. **Anatomy acts on load more than on partitioning.** Fermenting guts carry more *mcrA* and more *acsB*, as expected if longer retention and more fermentable substrate support larger hydrogenotrophic populations of both guilds.
3. **Phylogeny and its correlates.** Lineage-specific microbial associations [L07] could fix the methanogen:acetogen balance more strongly than gut architecture does.

With the present design we cannot distinguish these explanations. The within-study analysis, the only one free of between-study confounding, also showed no site effect, but it rests on 12 studies.

**Relation to prior work.** Marsupials show heterogeneous H₂-management potential across species [L01], and our cross-host results generalise that heterogeneity: categories defined by gut anatomy do not delimit H₂-sink partitioning. In rumen, H₂ uptake is shared among several guilds [L03]. The absence of a trade-off between *mcrA* and *acsB* is consistent with this, but gene abundance does not measure flux.

**Limitations.**
- **Marker abundance is genetic potential, not activity or flux.** *acsB* also occurs outside reductive acetogens.
- ***acsB* attribution.** The pre-planned phylogenetic placement of *acsB* reads could not be done because reads were deleted after quantification. It was replaced post hoc by best-hit taxonomic attribution and taxon-restricted re-analyses, which left all conclusions unchanged. Best-hit attribution cannot show whether a Bacillota *acsB* carrier is itself a reductive acetogen.
- **Normaliser.** The normaliser is bacterial *rpoB*, so *mcrA* copies per genome are relative to bacterial genomes.
- **Read sampling.** Each sample is 5 × 10⁶ forward reads. P1 used head sampling, which differed from multi-offset sampling for *acsB* in one of five tested runs.
- **Pre-specified gate.** Gate 1b failed its pre-specified rule on one of ten checks, and the protocol was changed accordingly (amendment v1.2).
- **P1 coverage.** P1 has two source studies and only ruminants.
- **Study imbalance in P2.** 176 of 189 inventory studies sampled a single fermentation category.
- **Metadata.** Captivity is unknown for a third of samples. Diet is coarse and missing for hosts without EltonTraits entries. Fermentation categories are assigned by lineage rules and do not capture within-category variation, e.g. among simple-gut herbivores.
- **Screening errors.** Mouse samples entered P2 through a screening gap (handled post hoc above). The screen's read-length estimate also missed 76-bp paired libraries.
- **Phylogeny.** The tree is a topology with Grafen branch lengths. The dated-tree sensitivity was **[NOT DONE]**, as was the in-silico marker specificity benchmark **[NOT DONE]**.

**Conclusion.** Within ruminants, the hindgut is consistently enriched for acetogenic relative to methanogenic potential. Across vertebrates, faecal H₂-sink partitioning potential follows host phylogeny more than gut fermentation architecture, and methanogenic and acetogenic potential do not trade off. Fermentation architecture predicts how much hydrogenotrophic potential a gut carries, not how it is divided.

## References
- [L01] Bowerman KL, et al. (2025) *mSystems* 11:e01608-25. doi:10.1128/msystems.01608-25
- [L02] Li Q, et al. (2025) *Microbiome* 13:28. doi:10.1186/s40168-024-02018-1
- [L03] Greening C, et al. (2019) *ISME J* 13:2617–2632. doi:10.1038/s41396-019-0464-2
- [L04] Xie F, et al. (2021) *Microbiome* 9:137. doi:10.1186/s40168-021-01078-x
- [L05] Youngblut ND, et al. (2020) *mSystems* 5:e01045-20. doi:10.1128/mSystems.01045-20
- [L07] Youngblut ND, et al. (2021) *Nat Microbiol* 6:1443–1454. doi:10.1038/s41564-021-00980-2
- [L08] Levin D, et al. (2021) *Science* 372:eabb5352. doi:10.1126/science.abb5352
- [L10] Kuhn M, et al. (2026) *Nucleic Acids Res* 54:D826–D834. doi:10.1093/nar/gkaf1118
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
- [TO ADD AND VERIFY] Proportionality ρp (Erb & Notredame 2016); Open Tree of Life synthetic tree and TNRS; Grafen branch lengths.
