# Research plan — HydroSink: a cross-vertebrate meta-analysis of gut hydrogen-sink partitioning

## Protocol record

The original planning and execution notes below record decisions available at the time. References to superseded scripts and working logs point to material in the [archived source snapshot](https://github.com/Angsumi/hydrosink/tree/v2.2.0) or the earlier protocol commits. The executed methods and results are in [the manuscript](manuscript/manuscript.md), [the reproducibility instructions](README.md) and the two protocol amendments.

The initial protocol and analysis commits are `67eb757` and `eec1a3e`. The `prereg-v1` tag follows those commits and includes Gate 1a results. P1 became a positive control for a previously published pattern in the same source datasets. The final P2 models used faeces, a phylogenetic covariance matrix and DIC. P3 pooled study-unit proportionality without removing host effects within units. The planned *acsB* placement was replaced with post hoc best-hit attribution.

### Status legend for the historical notes
**[VERIFIED]** checked against a primary source during this session ·
**[HYPOTHESIS]** proposed, untested ·
**[PLANNED]** analysis designed, code written, **not executed on real data** ·
**[TO VERIFY]** recalled or partially sourced; must be checked before use ·
**[NOT DONE]** nothing produced yet.

> **Historical status update (session 2, 2026-11-15, evening).** This machine has network access to ENA/NCBI and R 4.6. The environment, the marker database, the inventory and the PRISMA screening are now **executed** (see §6). Quantification has started with Gate 1a (the positive control). **No biological result is reported anywhere until it is computed.** The note below describes session 1.
>
> **Historical status (session 1).** Phases 1–2 (literature and dataset feasibility) were carried out with web search. Phase 3 code is written and has been checked for syntax. Only the Python random-effects module was test-run, and only on **synthetic data** to confirm the code works. **No public metagenome has been downloaded or profiled, and no biological result exists yet.** The execution environment had no network route to ENA/SRA (the HTTPS proxy refused the connection with a 403), and no R installation.

---

## Phase 1 — Literature and gap discovery

### 1.1 Search strategy (executed in part)
Sources: PubMed/PMC, publisher pages, bioRxiv/Research Square, DOAJ, the EMBL SPIRE/Metalog resources (web search; Web of Science and Scopus are not reachable from this environment).
Core queries used: hydrogen sink partitioning across animal guts; reductive acetogenesis in the hindgut; hydrogenases across mammal metagenomes; captive vs wild functional metagenomics; vertebrate gut metagenome catalogues.
**[TO DO before submission]** a formal, dated, reproducible Boolean search in PubMed + Web of Science + Scopus + bioRxiv, with a PRISMA 2020 flow log (`analysis/prisma_log.tsv`). The draft search string is in §2.2.

### 1.2 What the evidence establishes [VERIFIED; see literature.csv]
1. **Within ruminants, the hindgut favours acetogenesis.** In 16 goats, the reductive-acetogenesis marker *acsB* was ~12-fold more abundant in the cecum than in the rumen, while the rumen was enriched for methanogens and hydrogenotrophic respiratory bacteria. The authors report the same hindgut enrichment of acetogens in several ruminant species using the Xie et al. 2021 data (Li et al. 2025, *Microbiome*).
2. **Marsupials are not homogeneous in H₂ management.** Of 33 faecal samples from 14 mostly captive species, some grouped with high-methane mammals and others were enriched for bacterial H₂-uptake hydrogenases and alternative electron acceptors. The comparison set was cattle, capybara, guinea pig, horse and rabbit. n = 1–5 per species (Bowerman et al. 2025, *mSystems*).
3. **In rumen, H₂ uptake is shared among many guilds.** Methanogens, acetogens and fumarate, nitrate and sulfate reducers all take part. Alternative uptake transcripts are upregulated in low-methane sheep (Greening et al. 2019, *ISME J*).
4. **Host phylogeny explains archaeal diversity more strongly than diet across 110 vertebrate species.** This was shown with archaea-specific 16S amplicons only, not functional genes (Youngblut et al. 2021).
5. **Large shotgun catalogues of vertebrate guts exist** (Youngblut et al. 2020: 180 mostly wild species; Levin et al. 2021: >180 wild species; Xie et al. 2021: 370 samples, 10 GIT regions, 7 ruminants). None of them tests H₂-sink partitioning across hosts with phylogenetic control.
6. **A curated, harmonised metadata layer now exists.** Metalog has ~11,400 animal metagenome samples annotated with host species and captivity status, linked to SPIRE (Kuhn et al. 2026, *NAR*). This is what makes a cross-study design practical.

### 1.3 Gap
Evidence on which H₂ sink dominates comes from (a) single host lineages (ruminants, marsupials, humans), (b) within-host compartment contrasts in one species, or (c) cross-host surveys that used 16S, not function. **No study has tested, across vertebrates and with phylogenetic and study-level control, whether gut fermentation architecture predicts how H₂ disposal is split between methanogenesis and the alternative sinks, or whether methanogenic and acetogenic capacity trade off consistently.**

### 1.4 Candidate hypotheses (scored 1–5)

| # | Hypothesis | Novelty | Biological significance | Data availability | Feasibility | Publication potential | Total |
|---|---|---|---|---|---|---|---|
| H1 | **Fermentation architecture governs H₂-sink partitioning.** Foregut chambers favour methanogenesis; hindgut and caecal compartments favour reductive acetogenesis and respiratory hydrogenotrophy. Methanogenic and acetogenic capacity trade off across vertebrates, independent of host phylogeny and diet. | 4 (generalises single-lineage findings; no cross-vertebrate functional test exists) | 5 (methane mitigation, energy harvest, the theory of syntrophic H₂ economics) | 4 (Xie 2021 within-host compartments; ruminant, hindgut-fermenter, marsupial and broad vertebrate sets) | 4 (a targeted marker-gene approach is cheap; reads can be subsampled) | 4 (*ISME J* / *Microbiome*; *Nat Microbiol* if the trade-off holds with phylogenetic control) | **21** |
| H2 | Captivity causes convergent loss of fibre-degradation (CAZyme) capacity across species. | 3 (Alberdi 2021 covered the taxonomic side with 16S; no shotgun functional synthesis) | 4 (conservation and rewilding) | 2 (few shotgun wild/captive pairs; captivity is confounded with diet and study) | 3 | 3 | 15 |
| H3 | ARG burden in wild animals tracks proximity to humans and livestock, mediated by mobile elements. | 2 (crowded field) | 4 | 3 | 3 | 3 | 15 |

### 1.5 Selected question: **H1** [HYPOTHESIS]
**Rationale.** H1 rests on a thermodynamic mechanism (H₂ partial pressure, residence time and substrate supply differ systematically between foregut and hindgut). It has a natural within-host experiment: paired compartments from the same animal remove host genotype, diet and study as confounders. Its predictions can be falsified with a small, well-annotated set of marker genes. It also asks a mechanism question rather than producing a descriptive survey.

**Formal predictions**
- **P1 (within-host, primary):** In paired foregut vs hindgut samples from the same animal, log(acsB / mcrA) is higher in the hindgut, and the pooled within-study effect is greater than 0.
- **P2 (across hosts):** Fermentation site (foregut / caecal-colonic hindgut / simple gut) explains variance in the sink log-ratio beyond host phylogeny (Pagel's λ random effect), diet class, captivity and technical covariates.
- **P3 (trade-off):** Within studies, methanogenic and acetogenic capacities covary negatively. Measured on a compositional scale (proportionality ρ on CLR values), after removing study and host effects.
- **P4 (respiratory sinks, exploratory):** Respiratory H₂ uptake (group 1 [NiFe]-hydrogenase together with *nrfA* / *dsrAB* reductive type / *frdA*) is higher where methanogenic capacity is lower.

**What would falsify H1:** P1 pooled effect ≤ 0, or heterogeneous in sign across studies (95% prediction interval mostly ≤ 0); fermentation site adds nothing once phylogeny is in the model (ΔELPD within 2 SE); ρ(mcrA, acsB) ≥ 0.

---

## Phase 2 — Systematic dataset mining

### 2.1 Sources
NCBI SRA/BioProject, ENA, MGnify, **Metalog** (curated animal metadata including captivity status), **SPIRE** (pre-computed genes, MAGs and mOTUs profiles). Seed accessions come from the primary papers; see `datasets.csv`, where each row carries a verification status.

### 2.2 Draft search string (to run on ENA/SRA and literature databases)
`("shotgun metagenom*" OR "whole metagenome") AND (gut OR rumen OR caecum OR cecum OR colon OR faec* OR fec*) AND (mammal* OR bird* OR reptil* OR ruminant* OR marsupial* OR primate* OR equid* OR rodent* OR lagomorph*) NOT (human OR "Homo sapiens" OR mouse OR "Mus musculus")`
ENA filter: `library_strategy="WGS" AND library_source="METAGENOMIC" AND host_tax_id != 9606`.

### 2.3 Inclusion and exclusion criteria (PRISMA 2020-aligned)
**Include:** shotgun metagenomes (WGS/METAGENOMIC) from the GI contents or faeces of non-human vertebrates; Illumina reads ≥ 100 bp; ≥ 1 M read pairs after QC; host species identifiable to species level; sample site recorded (compartment, or faeces); publicly downloadable raw reads.
**Exclude:** amplicon, metatranscriptome or virome-enriched libraries; in-vitro fermenters and enrichment cultures; animals given methanogenesis inhibitors, ionophores or nitrate/sulfate/fumarate supplements (these perturb the outcome directly; retained only in a sensitivity set); antibiotic-treated animals; germ-free or gnotobiotic animals; duplicate uploads of the same biosample; pooled samples spanning several hosts; samples with > 50% host reads after decontamination.
**Overlap control:** deduplicate by BioSample and by MD5 of the first 10⁵ reads (catches re-deposits). When Youngblut 2020 re-used public projects, keep the original project record.
**Risk-of-bias items per study:** sample type (content vs faeces), storage and preservation, DNA extraction kit, library kit, sequencing depth, captivity, diet documentation, sampling design (paired vs unpaired).

### 2.4 Feasibility assessment
- **Paired within-host compartment data (P1):** Xie 2021 (PRJNA657455; 7 ruminant species × 10 regions; 370 samples) and Li 2025 (PRJNA1068803; 16 goats, rumen and cecum). **Sufficient for P1 on its own.** Limitation: ruminants only.
- **Cross-host (P2/P3):** hindgut fermenters (horse PRJNA590977; rabbit, guinea pig and capybara through the projects cited by Bowerman 2025), marsupials (PRJNA1167890), broad vertebrate taxa (Youngblut 2020 PRJEB38078). The Metalog animal subset (~11.4k samples) is the main source for expansion. Expected after filtering: **several hundred to ~2,000 samples, ≥ 40 host species** [estimate, to confirm by running `01_build_inventory.py`].
- **Compute:** marker-only DIAMOND blastx on reads subsampled to 5 M pairs per sample costs roughly 5–15 CPU-min per sample against a small (< 50k sequence) marker database. That is about 300–500 CPU-h for 2,000 samples, which is modest.
- **Key risks:** (i) faeces vs gut contents differ (handled as a covariate plus a contents-only sensitivity set); (ii) *acsB* homologues outside acetogens (handled by a phylogenetic placement filter); (iii) archaea at < 1% in many hosts give zero-inflated mcrA (handled with a hurdle / zero-inflated model); (iv) diet and fermentation site are collinear across hosts, which is why the within-host P1 design is primary.
- **Verdict:** **Feasible**, provided the analysis runs on infrastructure with SRA/ENA access.

---

## Phase 3 — Reproducible analysis [PLANNED; code in `analysis/`]

| Step | Tool | Script |
|---|---|---|
| Inventory and metadata harmonisation | ENA Portal API, Metalog download | `scripts/01_build_inventory.py` |
| Download and subsample | `fasterq-dump`/ENA FTP, `seqkit sample` | `scripts/02_fetch_subsample.sh` |
| QC | `fastp` | `02_fetch_subsample.sh` |
| Marker database | HydDB classes (Søndergaard et al. 2016), curated KOfam HMMs for mcrA (K00399), acsB (K14138), dsrA/B (K11180/K11181), nrfA (K03385), frdA (K00244), narG (K00370), napA (K02567) | `scripts/03_build_marker_db.md` |
| Quantification | DIAMOND blastx (≥ 50% id, ≥ 80% query cover, e < 1e-10) plus group-specific thresholds | `scripts/04_quantify_markers.sh` |
| Normalisation | genome equivalents from universal single-copy genes (SingleM / MicrobeCensus); abundance as copies per genome equivalent | `04_quantify_markers.sh` |
| Taxonomy | mOTUs3 (use SPIRE pre-computed profiles where available) | — |
| Host phylogeny | TimeTree 5 / VertLife consensus; Pagel's λ | `scripts/05_meta_analysis.R` |
| Models | brms phylogenetic GLMM (primary); metafor REML + Hartung–Knapp two-stage random effects (P1); proportionality ρ (P3) | `05_meta_analysis.R`, `scripts/random_effects.py` |
| Multiplicity | Benjamini–Hochberg across marker families (q < 0.05) | both |
| Sensitivity | leave-one-study-out; contents-only; wild-only; depth-matched (2 M pairs); alternative identity thresholds (40/60%); excluding birds; excluding any inhibitor/supplement studies | `05_meta_analysis.R` |
| Validation | train on all but marsupials, then predict the marsupial sink ratio (out-of-sample hold-out; Bowerman 2025 as an independent test set) | `05_meta_analysis.R` |

Executed in this session: `random_effects.py` self-test on **simulated** data (it exercises the DerSimonian–Laird / REML and Hartung–Knapp code paths only). Log: `analysis/selftest_log.txt`.

## Phase 4 — Interpretation [NOT DONE]
Planned framing: association vs causation (gene potential ≠ flux, so we will compare against the published in-vitro H₂/CH₄ measurements in Li 2025); the thermodynamic explanation (methanogens' lower H₂ threshold vs acetogens' tolerance of higher pH₂ in hindguts); and generalisability (ruminant-heavy sampling).

## Phase 5 — Manuscript [skeleton written; results sections deliberately empty]
`manuscript/manuscript.md`: Title, Abstract (structured, results placeholders), Introduction, Methods (complete), Results (planned analyses only), Discussion (framework), References (verified only).

## Decision log
| Date | Decision | Reason |
|---|---|---|
| 2026-11-15 | Selected H1 over H2 and H3 | Highest score; a within-host design breaks diet–phylogeny confounding |
| 2026-11-15 | Marker-gene approach instead of full assembly | Compute cost; the outcome is defined by a few genes |
| 2026-11-15 | Within-host compartment contrast (P1) is primary | Strongest design against confounding |
| 2026-11-15 | Halted before data processing | No network route to SRA/ENA and no R in the execution environment |

## Original next steps before execution
1. `conda env create -f analysis/env/environment.yml`
2. `python analysis/scripts/01_build_inventory.py --seeds datasets.csv --metalog <metalog_animal.tsv> --out analysis/inventory.tsv`, then complete the PRISMA log
3. Build the marker database (`03_build_marker_db.md`); validate it on Li 2025 to reproduce their ~12× acsB cecum enrichment (positive control)
4. Run steps 02 → 05; fill in the manuscript Results from the outputs only


---

## 6. Execution log — session 2 (2026-11-15)

### 6.1 Infrastructure (verified)
- 16 cores, 26 GB RAM, 199 GB disk. Env `hydrosink` (micromamba; `analysis/env/environment.yml`): DIAMOND 2.1.10, fastp 0.23.4, R 4.3 with metafor 4.8.0, MCMCglmm 2.36, ape, phytools, glmmTMB.
- **Bandwidth ≈ 4 MB/s total** (ENA HTTPS ≈ 0.6–1 MB/s per stream; NCBI S3 similar in aggregate). Full runs (7–18 GB) are infeasible at scale.

### 6.2 Decisions (added to the decision log)
| Date | Decision | Reason |
|---|---|---|
| 2026-11-15 | Continue H1; no restart of Phase 1 | 2026 novelty re-check (web + Europe PMC) found no cross-vertebrate, phylogenetically controlled functional test of H₂-sink partitioning |
| 2026-11-15 | Stream the first 5 M **forward** reads per run instead of a full download | Bandwidth; the Greening-lab protocol counts forward reads only. Head-of-file bias is tested in Gate 1b |
| 2026-11-15 | Marker reference = Greening lab compiled DB v3 (51 genes, all kept as competitors) instead of a from-scratch KOfam build | Curated, published, CC-BY; per-gene read thresholds published (Ricci et al. 2025, L23) |
| 2026-11-15 | Normaliser = full-length UniRef90 bacterial RpoB (5,635 seqs) with RpoC decoys; rejected UniRef50 because its representatives were mostly fragments (median 135 aa) | Accurate gene length for RPKM; copies per genome = RPKM/RPKM_rpoB. Archaeal RpoB is not included (a slight overestimate of copies per genome where archaea are abundant; noted as a limitation) |
| 2026-11-15 | SingleM and brms dropped; MCMCglmm used for the phylogenetic GLMM | Install reliability; equivalent model class |
| 2026-11-15 | P2 restricted to **faecal** samples | The only sample type shared by all hosts; gut-content samples go into a sensitivity analysis |
| 2026-11-15 | P1 restricted to fermentation chambers (rumen, reticulum, omasum vs cecum, colon, rectum) | Abomasum and small intestine are not fermentation chambers and have low biomass / high host DNA |
| 2026-11-15 | Runs without a structured `host_tax_id` (48,106) excluded | Free-text hosts are mostly domestic pig, cattle, chicken and dog, which are already at the species cap. Limitation: rare wild hosts labelled only in free text are missed |
| 2026-11-15 | Fermentation site assigned by lineage rules (Stevens & Hume 1998, L20) **before** any marker data were inspected | Blinded predictor assignment; `analysis/scripts/01b_host_traits.py` |

### 6.3 Executed outputs so far
- `analysis/db/` — hydrosink.faa / .dmnd / labels.tsv (170,171 sequences); sources and md5 values in config.yaml.
- `analysis/inventory/raw/` — ENA candidate run tables (search date in `search_date.txt`), host lineages, EltonTraits.
- `analysis/inventory/host_traits.tsv` — 488 vertebrate host taxa with fermentation site and diet.
- `analysis/prisma_log.tsv` — screening counts: 136,708 runs identified → 19,385 eligible → P2 930 faecal runs (183 species, 189 studies); P1 254 runs (53 animals; Li 2025 + Xie 2021).
- Pipeline smoke test (SRR28204437, goat cecum): 5 M reads, 4,704 rpoB reads; acsB 0.30 and mcrA 0.004 copies per genome. 14 min wall time (12 min streaming, 2 min DIAMOND on 14 threads).
