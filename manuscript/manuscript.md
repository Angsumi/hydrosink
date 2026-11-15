# Gut fermentation architecture governs hydrogen-sink partitioning across vertebrates: a phylogenetically controlled meta-analysis of public shotgun metagenomes

*Pre-analysis manuscript (registered-report style). The Methods are complete. **The Results and the quantitative parts of the Abstract and Discussion are intentionally empty: no data have been analysed yet.** Text in [brackets] marks content that must come from executed analyses.*

Target journals: *ISME Journal* / *Microbiome* (*Nature Microbiology* if P2 and P3 hold with phylogenetic control). Reporting: PRISMA 2020 (adapted for data synthesis), STORMS checklist for microbiome studies.

---

## Abstract
**Background.** Molecular hydrogen (H₂) produced by fermentation must be consumed for fermentation to stay thermodynamically favourable. Methanogens, reductive acetogens and respiratory hydrogenotrophs compete for it, and who wins shapes methane emissions and host energy harvest. Single-lineage studies point to acetogenesis dominating ruminant hindguts and to varied H₂ management in marsupials, but no cross-vertebrate test exists.
**Hypothesis.** Gut fermentation architecture, rather than host phylogeny or diet alone, determines how H₂ disposal is split: methanogenesis in foregut chambers, acetogenesis and respiratory sinks in hindguts. Methanogenic and acetogenic capacity trade off.
**Methods.** Systematic, PRISMA-aligned mining of public shotgun gut metagenomes from non-human vertebrates (ENA/SRA, Metalog/SPIRE). Curated marker-gene quantification (mcrA, acsB, dsrAB, nrfA, narG/napA, frdA, HydDB hydrogenase classes) normalised to genome equivalents. A within-animal paired compartment contrast pooled by random-effects meta-analysis. Cross-host Bayesian phylogenetic mixed models. Compositional proportionality for trade-offs. Out-of-sample validation on marsupials.
**Results.** [n studies, n samples, n host species; pooled within-host effect with 95% CI and prediction interval; ΔELPD for fermentation site; phylogenetic signal; pooled ρ(mcrA, acsB); validation performance.]
**Conclusions.** [To be written from the results.]

## Introduction
Fermentation in vertebrate guts releases H₂. If H₂ accumulates, NADH re-oxidation is inhibited and fermentation stalls, so hydrogenotrophs keep partial pressures low [L03]. Several guilds compete for this H₂. Hydrogenotrophic methanogens convert it to methane, which is an energy loss for the host and a greenhouse gas. Reductive acetogens make acetate, which the host can absorb. Fumarate, nitrate and sulfate reducers couple H₂ oxidation to respiration [L03]. In the rumen, alternative uptake pathways are upregulated in low-methane sheep [L03]. Within goats, acsB is about 12-fold enriched in the cecum relative to the rumen, and the authors report the same direction across several ruminant species [L02]. Marsupials show diverse H₂ management between and within species [L01].

Two gaps remain. First, every functional test so far is confined to one host lineage. Cross-vertebrate work on methanogens relies on 16S amplicons and finds host phylogeny to be a strong determinant [L07], which leaves open whether functional sink partitioning follows phylogeny or gut architecture. Second, diet, phylogeny and gut anatomy are collinear across species, so cross-sectional comparisons cannot separate them. Paired compartments from the same animal provide a natural control. Large vertebrate shotgun catalogues [L04, L05, L08] and curated cross-study metadata [L10] now make a phylogenetically controlled synthesis possible.

We test four pre-specified predictions (P1–P4; research_plan.md §1.5).

## Methods
### Protocol and registration
Hypotheses, inclusion criteria, thresholds and models were fixed in `analysis/config/config.yaml` and `research_plan.md` before any data were downloaded. [Register on OSF before data access; add DOI.]

### Data sources and search
ENA/SRA were queried with the strategy in research_plan.md §2.2. Seed projects came from the primary literature (datasets.csv). The Metalog animal subset [L10] supplied curated host and captivity metadata. Search date: [date].

### Eligibility
As in research_plan.md §2.3. Duplicates were removed by BioSample and by a read-level MD5 prefix. Each study was scored for risk of bias on sample type, preservation, extraction kit, depth, captivity, diet documentation and design.

### Processing
Reads were QC'd with fastp (Q20, length ≥ 60, deduplication) and subsampled to 5 × 10⁶ pairs (seed 20261115; 2 × 10⁶ for sensitivity). Marker reads were identified with DIAMOND blastx against a curated, guild-labelled reference (`03_build_marker_db.md`). Ambiguous families (acsB, dsrAB, frdA/sdhA) were resolved by phylogenetic placement, and guild thresholds were chosen to give simulated-read precision ≥ 0.95. Abundance was expressed as reads per kilobase per genome equivalent, with genome equivalents from SingleM single-copy marker coverage. Taxonomy came from mOTUs3 (SPIRE profiles where available). Host traits: diet class (EltonTraits), fermentation site (foregut / hindgut / simple / avian caecal; assigned from the literature and blinded to outcome), and captivity (Metalog). Host phylogeny: TimeTree/VertLife consensus.

### Outcomes
Primary: ALR_acet = log[(acsB + c)/(mcrA + c)], with c set to half the minimum non-zero value per marker. Secondary: ALR_resp (summed respiratory reductases vs mcrA), mcrA presence (hurdle), and hydrogenase class profile.

### Statistical analysis
**P1:** Within-animal differences (hindgut minus foregut) were pooled by study. This used random-effects meta-analysis (REML τ², Hartung–Knapp–Sidik–Jonkman CI, 95% prediction interval, I², leave-one-out), and in parallel a Bayesian mixed model with random slopes by study and random intercepts for animal and host. Python and R implementations were cross-checked. **P2:** brms Gaussian models with a phylogenetic random effect (Brownian correlation matrix), host-species, and study random effects. Fixed effects: diet, captivity, sample type, log depth and read length. Models with and without fermentation site were compared by PSIS-LOO; phylogenetic signal is reported as the variance ratio (λ analogue). **P3:** Within-study proportionality ρ(mcrA, acsB) on CLR-transformed abundances, Fisher-z pooled by random effects. **Multiplicity:** Benjamini–Hochberg across marker families (q < 0.05). **Sensitivity:** contents-only, wild-only, excluding birds, depth-matched, identity thresholds 40/60%, excluding any perturbation studies, leave-one-study-out. **Validation:** model trained without marsupials, then used to predict marsupial ALR_acet (Spearman ρ, RMSE). **Positive control:** reproduce the acsB cecum enrichment in Li et al. [L02] before the main analyses.

### Causality statement
Marker abundance measures genetic potential, not flux. We interpret associations as consistent with, not proof of, differential H₂ disposal, and compare them with published direct measurements [L02].

### Data and code availability
All accessions are in datasets.csv and `analysis/inventory.tsv`. The code and environment are in `analysis/` [archive on Zenodo with a DOI at submission].

## Results
*Not yet available. Planned order:*
1. Search and screening (PRISMA flow, Fig. 2) and dataset characteristics (Table 1).
2. Positive control and marker benchmark (Supplementary Fig. S1).
3. P1 within-host contrast: forest plot (Fig. 3).
4. P2 cross-host model, phylogenetic signal, and hurdle results (Fig. 4).
5. P3 trade-off (Fig. 5) and P4 respiratory sinks.
6. Sensitivity and marsupial validation (Supplementary Figs. S2–S4).

## Discussion (framework, to be completed from results)
- If P1 and P2 hold: H₂ economics follow anatomy convergently across lineages. This would suggest hindgut acetogenesis as a general alternative electron sink and has implications for methane mitigation (promoting acetogens).
- If P1 holds but P2 does not: compartment effects exist within hosts, but phylogenetic inertia dominates between hosts. This would be consistent with the 16S-based archaeal results [L07].
- Limitations to address: ruminant-heavy paired data; faeces as a proxy for the hindgut; small n for many wild species; gene potential ≠ activity; database bias against non-model hosts [L05, L08]; captivity confounding [L09].

## References
Only sources verified during planning are cited; full list in `literature.csv`. Items marked TO VERIFY there must be checked before submission.
- [L01] Bowerman KL, et al. (2025) *mSystems* 11:e01608-25. doi:10.1128/msystems.01608-25
- [L02] Li Q, et al. (2025) *Microbiome* 13:28. doi:10.1186/s40168-024-02018-1
- [L03] Greening C, et al. (2019) *ISME J* 13:2617–2632. doi:10.1038/s41396-019-0464-2
- [L04] Xie F, et al. (2021) *Microbiome* 9:137. doi:10.1186/s40168-021-01078-x
- [L05] Youngblut ND, et al. (2020) *mSystems* 5:e01045-20. doi:10.1128/mSystems.01045-20
- [L07] Youngblut ND, et al. (2021) Vertebrate host phylogeny influences gut archaeal diversity. *Nat Microbiol* [DOI to verify]
- [L08] Levin D, et al. (2021) *Science* 372:eabb5352. doi:10.1126/science.abb5352
- [L09] Alberdi A, et al. (2021) *Sci Rep* 11:22660. doi:10.1038/s41598-021-02015-6
- [L10] Kuhn M, et al. (2026) *Nucleic Acids Res* 54:D826–D834. doi:10.1093/nar/gkaf1118
