# PRISMA 2020 checklist — HydroSink (adapted to a synthesis of public metagenome data)

Item numbers and topics follow the official PRISMA 2020 checklist (Page et al. 2021 [L31]; file
`prisma_src/PRISMA_2020_checklist.docx`, prisma-statement.org). The unit of synthesis is a sequencing run (sample) and a
data-generating study (BioProject), not a published effect estimate; items written for effect-size reviews are
answered for that design. "MS" = `manuscript/manuscript.md`.

| # | Topic | Where / response |
|---|---|---|
| 1 | Title | MS title: identified as a pre-registered synthesis of public metagenomes (not labelled "systematic review"; the search is systematic over data repositories, not literature). |
| 2 | Abstract | MS Abstract (structured). PRISMA-Abstracts items not separately checked. |
| 3 | Rationale | MS Introduction ¶1–2. |
| 4 | Objectives | MS Introduction (P1–P4). |
| 5 | Eligibility criteria | MS Methods "Data sources and screening"; `analysis/config/config.yaml`; `research_plan.md` §2.3. Grouping for syntheses: P1 (paired compartments), P2 (faecal, cross-host), P3 (studies ≥ 10 samples). |
| 6 | Information sources | ENA Portal API, searched 15 Nov 2026 (`analysis/inventory/raw/search_date.txt`); seed studies from literature (`datasets.csv`). No literature databases searched for data. |
| 7 | Search strategy | ENA query parameters in `config.yaml` (`inclusion`) and `analysis/scripts/01_build_inventory.py`. |
| 8 | Selection process | Automated rules (`01_build_inventory.py`); one analyst plus scripted rules; no independent dual screening (limitation). Captivity curation: rules + 38 manual study-level decisions with evidence quotes (`config/captivity_manual.tsv`). |
| 9 | Data collection process | Reads fetched programmatically from NCBI SRA (AWS); metadata from ENA/BioSample APIs; no contact with study authors. |
| 10a | Data items — outcomes | MS Methods "Outcomes" (marker copies per genome; ALR_acet primary). All runs passing QC were used; none selected by result. |
| 10b | Data items — other variables | Fermentation site, diet (EltonTraits), captivity, host phylogeny, depth, read length, platform (MS Methods). Missing diet coded "unknown". |
| 11 | Study risk of bias assessment | **Not done.** Risk-of-bias items were listed in `research_plan.md` §2.3 but not scored per study. Technical study effects are modelled (study random effect; platform, depth, read-length covariates). |
| 12 | Effect measures | Difference in ln ratio (ALR) between chambers/sites; proportionality ρp (MS Methods). |
| 13a | Synthesis eligibility | MS Methods (P1 strata; P2 analysis set; P3 units ≥ 10 samples). |
| 13b | Data preparation | Pseudocount = half minimum non-zero; QC rpoB ≥ 200; genus-only hosts without tree tip excluded (MS Methods). |
| 13c | Tabulation / display | Figs 2–4, S1–S2; Table 1 (`results/table1_*.tsv`). |
| 13d | Synthesis methods | REML + Hartung–Knapp (P1), MCMCglmm phylogenetic model (P2), Fisher-z pooled ρp (P3) (MS Methods). |
| 13e | Heterogeneity | I², prediction intervals, three-level model (P1); variance partition (P2); I² (P3). |
| 13f | Sensitivity analyses | Pre-specified list + identity thresholds, dated tree, acsB attribution, post-hoc no-mouse set (MS Methods/Results). |
| 14 | Reporting bias assessment | Not applicable in the usual sense (raw data, not published estimates); not assessed. Selection of studies by data availability in ENA is a limitation. |
| 15 | Certainty assessment | **Not done** (no GRADE-type assessment). |
| 16a | Study selection | Fig. S2; `analysis/prisma_log.tsv`, `results/prisma_final.tsv`. |
| 16b | Excluded studies | Exclusions are counted per rule (Fig. S2); excluded studies are not individually cited (tens of thousands of runs). |
| 17 | Study characteristics | `results/table1_p2.tsv`, `table1_p1.tsv`; per-run metadata `results/p2_analysis_set.tsv`. |
| 18 | Risk of bias in studies | **Not done** (see 11). |
| 19 | Results of individual studies | P1 per stratum (`p1_stratum_effects.tsv`, Fig. 2b); P3 per unit (`p3_units.tsv`, Fig. 4); P2 per sample in `p2_analysis_set.tsv`. |
| 20a | Characteristics of syntheses | MS Results "Data set". |
| 20b | Results of syntheses | MS Results P1–P4 (estimates, CIs/HPDs, prediction intervals). |
| 20c | Heterogeneity results | MS Results (I², variance partition). |
| 20d | Sensitivity results | MS Results (sensitivity tables). |
| 21 | Reporting biases | Not assessed (see 14). |
| 22 | Certainty of evidence | Not assessed (see 15). |
| 23a | Interpretation | MS Discussion. |
| 23b | Limitations of evidence | MS Discussion "Limitations". |
| 23c | Limitations of review process | MS Discussion "Limitations" (single analyst, automated screening, screening errors found). |
| 23d | Implications | MS Discussion "Conclusion". |
| 24a | Registration | **Not registered** in an external registry; time-stamped git tags `prereg-v1`, `-v1.1`, `-v1.2`. |
| 24b | Protocol | `research_plan.md`, `analysis/config/config.yaml`, `analysis/scripts/05_meta_analysis.R` at tag `prereg-v1`. |
| 24c | Amendments | `analysis/PREREG_AMENDMENT_v1.1.md`, `v1.2.md`; post-hoc analyses labelled in MS. |
| 25 | Support | Self-funded; no external funding (MS Declarations). |
| 26 | Competing interests | None declared (MS Declarations). |
| 27 | Availability of data, code, materials | MS "Data and code availability"; all code in `analysis/scripts/`, per-sample outputs and hits with SHA-256 manifests. |
