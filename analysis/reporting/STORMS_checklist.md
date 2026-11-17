# STORMS checklist (v1.03) — HydroSink

Items follow STORMS v1.03 (Mirzayi et al. 2021 [L32]; `storms_src/STORMS_Excel_1.03.xlsx`, Zenodo 10.5281/zenodo.5714305).
STORMS is written for primary human microbiome studies; this is a secondary analysis of public animal metagenomes, so
laboratory and participant items refer to the source studies and are answered "not available / source study".
"MS" = `manuscript/manuscript.md`.

| # | Item | Response |
|---|---|---|
| 1 | Abstract | MS Abstract (structured). |
| 1.1 | Study design | Abstract: pre-registered synthesis; within-animal meta-analysis + cross-host phylogenetic models. |
| 1.2 | Sequencing methods | Abstract: read-based marker-gene quantification (DIAMOND blastx) of shotgun reads. |
| 1.3 | Specimens | Abstract: paired foregut/hindgut contents (P1), faeces (P2). |
| 2 | Background and rationale | MS Introduction. |
| 2.1 | Hypotheses | MS Introduction (P1–P4, pre-specified). |
| 3 | Study design | MS Methods "Protocol"; secondary analysis of public data. |
| 3.1 | Participants | Host animals of public studies; population = non-human vertebrates with public gut shotgun metagenomes (MS Methods). |
| 3.2 | Geographic location | Per-run country in `results/p2_matrix.tsv` (`country`, where provided by submitters). |
| 3.3 | Relevant dates | Search date 15 Nov 2026; source sampling dates not harmonised (not available consistently). |
| 3.4 | Eligibility criteria | MS Methods "Data sources and screening". |
| 3.5 | Antibiotics usage | Treatment studies excluded by keywords; individual antibiotic exposure not available. |
| 3.6 | Analytic sample size | Fig. S2; `results/prisma_final.tsv`. |
| 3.7 | Longitudinal studies | Not applicable (one run per BioSample). |
| 3.8 | Matching | P1 pairs chambers within animal (MS Methods). |
| 3.9 | Ethics | Not applicable (public data); source-study ethics apply. |
| 4 | Laboratory methods | Source studies; not re-done. |
| 4.1 | Specimen collection | Body site from ENA metadata (faeces; P1 chambers from sample aliases). |
| 4.2 | Shipping | Not available (source studies). |
| 4.3 | Storage | Not available (source studies). |
| 4.4 | DNA extraction | Not harmonised (source studies); study modelled as random effect. |
| 4.5 | Host DNA depletion / microbial enrichment | Not known per study; normalisation to rpoB reduces dependence on host-read fraction. |
| 4.6 | Primer selection | Not applicable (shotgun). |
| 4.7 | Positive controls | Pipeline positive control: P1 reproduction of [L02]; in-silico benchmark (simulated reads from 43 genomes). |
| 4.8 | Negative controls | Not available from source studies; host-genome decoy in the in-silico benchmark. |
| 4.9 | Contaminant mitigation | Not performed beyond QC; limitation. |
| 4.10 | Replication | Reproducibility: 17 runs re-processed from scratch, byte-identical outputs. |
| 4.11 | Sequencing strategy | Illumina shotgun (WGS, METAGENOMIC). |
| 4.12 | Sequencing methods (quantification) | Relative: copies per genome (marker RPKM / rpoB RPKM). |
| 4.13 | Batch effects | Study random effect; platform covariate sensitivity; within-study analysis. |
| 4.14 | Metatranscriptomics | Not applicable. |
| 4.15 | Metaproteomics | Not applicable. |
| 4.16 | Metabolomics | Not applicable. |
| 5 | Data sources / measurement | Host traits: fermentation site (lineage rules [L20]), diet (EltonTraits [L21]), captivity (curated), phylogeny (Open Tree [L34]; TimeTree [L36] for sensitivity). |
| 6 | Causal inference | MS Discussion: associations of genetic potential only; confounding of anatomy with phylogeny, diet and study discussed; within-study analysis. |
| 6.1 | Selection bias | MS Limitations (repository availability, species cap, screening errors). |
| 7 | Bioinformatic and statistical methods | MS Methods (ALR with pseudocount; CLR for ρp). |
| 7.1 | Quality control | fastp settings; rpoB ≥ 200 reads. |
| 7.2 | Sequence analysis | DIAMOND blastx against Greening-lab DB v3 + RpoB [L22, L25]; per-gene thresholds [L23]. |
| 7.3 | Statistical methods | MS Methods "Statistical analysis". |
| 7.4 | Longitudinal analysis | Not applicable. |
| 7.5 | Subgroup analysis | P2 sensitivity subsets (mammals, herbivores, non-domestic, known captivity). |
| 7.6 | Missing data | Diet missing → "unknown" category; captivity unknown → category; genus-only hosts without tree tip excluded. |
| 7.7 | Sensitivity analyses | MS Methods/Results (pre-specified and post-hoc, labelled). |
| 7.8 | Findings | All pre-specified outcomes reported regardless of result. |
| 7.9 | Software | MS Methods (versions); `analysis/results/sessionInfo.txt`; `analysis/env/environment.yml`. |
| 8 | Reproducible research | Full pipeline scripts and git history; `README.md`, `run_all.sh`. |
| 8.1 | Raw data access | ENA/SRA run accessions in `analysis/inventory/p1_runs.tsv`, `p2_runs.tsv`. |
| 8.2 | Processed data access | `analysis/results/per_sample/`, `analysis/results/hits/` (SHA-256 manifests). |
| 8.3 | Participant data access | Host metadata per run: `results/p2_analysis_set.tsv`, `p1_matrix.tsv`. |
| 8.4 | Source code access | `analysis/scripts/` [archive DOI at submission]. |
| 8.5 | Full results | `analysis/results/*.tsv`. |
| 9 | Descriptive data | Table 1 (`table1_p2.tsv`, `table1_p1.tsv`). |
| 10 | Microbiome data | MS Results; Fig. S1. |
| 10.1 | Taxonomy | Not a taxonomic study; acsB read attribution uses NCBI taxonomy (`results/acsb/`). |
| 10.2 | Differential abundance | Marker contrasts with direction and intervals (`p1_pooled.tsv`, `p2_contrasts.tsv`). |
| 10.3 | Other data types | Functional potential (marker genes) only. |
| 10.4 | Other statistical analysis | P3 proportionality; benchmark; reproducibility audit. |
| 11 | Key results | MS Discussion "Summary". |
| 12 | Interpretation | MS Discussion "Interpretation". |
| 13 | Limitations | MS Discussion "Limitations". |
| 13.1 | Bias | MS Limitations. |
| 13.2 | Generalizability | MS Limitations (ruminant-only P1; faeces only P2). |
| 14 | Ongoing/future work | MS Discussion [brief; to expand]. |
| 15 | Funding | [AUTHORS TO COMPLETE] |
| 15.1 | Acknowledgements | [AUTHORS TO COMPLETE] |
| 15.2 | Conflicts of interest | [AUTHORS TO COMPLETE] |
| 16 | Supplements | Figs S1–S2, supplementary tables = `analysis/results/` files named in MS. |
| 17 | Supplementary data | Per-sample marker tables (all genes) in `results/per_sample/`. |
