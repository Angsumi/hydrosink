# Pre-registration amendment v1.1 (2026-11-16)

**Status at time of amendment:** 0 / 930 P2 samples quantified, no `p1_matrix.tsv` / `p2_matrix.tsv` assembled;
P1: 56 / 254 per-sample marker tables exist but have not been assembled or inspected. Only the Gate 1a positive
control (Li 2025 goats; acsB, mcrA per pair) has been examined. Amended code was tested only on synthetic random data
(`scratchpad`, not part of the project). Base version: git tag `prereg-v1`.

Motivation: an external critique of the design (six points, summarised below). Every change is listed with the file
that implements it. Nothing in the primary outcome (ALR_acet), primary contrast, marker thresholds or search
parameters changed.

## 1. Claims and framing (manuscript, no code change)
| Critique | Response |
|---|---|
| Li/Xie paired result already published (Li et al. 2025 Microbiome) | P1 is re-labelled **positive control / calibration**, not a novel finding. Novelty claim rests on P2 only. |
| Bowerman et al. 2025 (mSystems) already compared H2-management across marsupials and other mammals | Cited as the closest prior work; our claim is a **larger, phylogenetically controlled test of generality**, framed as a tendency (prediction intervals and within-category variance reported), not a rule. Never "first cross-host comparison". |
| Faeces do not represent foregut metabolism | P2 conclusions restricted to **faecal genetic potential**. Whole-gut H2 disposal claims only from P1 compartment data. |
| Single markers do not show flux/pathway completeness | Wording restricted to "genetic potential / marker abundance". |
| Gate 1a validates acsB only; mcrA contrast is null (8/16 pairs, p = 0.60) | Reported as such. mcrA/rpoB mixes archaeal marker with a bacterial normaliser (DB rpoB is bacterial); this is stated as a limitation. |

## 2. Analysis-code changes (`analysis/scripts/05_meta_analysis.R`)
1. **P1 = two studies.** New output `p1_per_study.tsv`: per-study estimates (Li 2025; Xie 2021) and a three-level model
   (`rma.mv`, strata nested in study). The stratum-level pooling (`p1_pooled.tsv`) is kept but labelled descriptive;
   study x species strata are not presented as independent replication.
2. **Within-study design check** (`p2_within_study.tsv`): only studies sampling >= 2 fermentation categories
   (12 studies in the inventory); study as a fixed effect, species random; estimates the category contrast purely
   within studies. Addresses study/anatomy confounding (176 / 189 studies contain one category).
3. **Marsupial holdout** (`p2_holdout_marsupial*.tsv`): model fitted without Metatheria; marsupial samples predicted
   from fixed effects with full predictive variance; reports species-level Spearman, MAE, 95 % coverage and
   observed vs predicted hindgut-foregut difference. **Limitation:** the inventory contains only 21 marsupial samples
   from 5 species (species cap), so this is a weak validation and will be described as such.
4. **Known-captivity subset** and **sequencing-platform covariate** (`p2_sensitivity_v1_1.tsv`).
5. Results directory and MCMC length can be overridden by environment variables (testing only; defaults unchanged).

## 3. Metadata curation (outcome-blind; `analysis/scripts/01d_curate_metadata.py` -> `inventory/p2_curation.tsv`)
- Sources: ENA study title/description (187 / 189 studies) and BioSample attributes (530 / 930 samples retrieved).
- Keyword rules (wild / captive / domestic / laboratory), sample text first, then study text; "wild-type" excluded;
  study-level "wild" never overrides a domestic species. Manual per-study (or study x species) decisions with an
  evidence quote: `analysis/config/captivity_manual.tsv` (38 entries). Rule calls were reviewed per study and
  misclassifications corrected (e.g. "wild-type goats", "Public Health Laboratory", mixed wild/captive studies -> unknown).
- Captivity: unknown **751 -> 296** of 930 (wild 327, domestic 126, captive 123, laboratory 58).
- Host labels: 5 genus-only labels resolved where unambiguous (`analysis/config/species_manual.tsv`: lab Mus/Rattus,
  monotypic genera). Genus-only hosts without a tree tip are **excluded from P2 models** by construction (47 samples
  total excluded at the phylogeny step), so they reduce sample size but cannot bias the phylogenetic model.
- Integrated in `05a_assemble.py` (curated values replace v1 values; v1 values retained as `captivity_v1`, `species_v1`).

## 4. Pipeline robustness (no effect on results)
- RAM guard (max 3 concurrent DIAMOND jobs, `-b 1 -c 4`); resumable byte-range downloads (same first-N reads).

## 5. Planned but NOT yet implemented (will be reported as not done if not completed)
- In-silico marker specificity benchmark (reads simulated from genomes with/without each marker).
- Phylogenetic placement of acsB hits (acetogen vs non-acetogen carriers).
- Dated-tree sensitivity (VertLife posterior trees) for the P2 phylogeny.
- Gate 1b (head-of-file vs random subsample) — running; decision rule committed before data (`gate1b_compare.py`).
- Rewrite of `manuscript/manuscript.md` Methods from the executed pipeline (current text is the v0 draft).
