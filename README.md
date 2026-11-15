# HydroSink — reproducibility package

**Status: plan + pipeline. No real data have been processed and no biological results exist yet.**

| Path | Contents | State |
|---|---|---|
| `research_plan.md` | Gap analysis, 3 scored hypotheses, selected question (H1), PRISMA criteria, feasibility, analysis plan, decision log | Complete |
| `literature.csv` | 19 references with DOI/URL and a verification status | 12 verified, 7 partial / to verify |
| `datasets.csv` | 15 seed dataset records with accessions, metadata fields and source links | Accessions verified from papers/SPIRE; run-level metadata to extract with step 01 |
| `analysis/env/environment.yml` | Conda environment | Written, not solved |
| `analysis/config/config.yaml` | Pre-registered thresholds | Complete |
| `analysis/scripts/01_build_inventory.py` | ENA query + inclusion filters + PRISMA counts | Syntax-checked; needs network |
| `analysis/scripts/02_fetch_subsample.sh` | Download, fastp QC, paired subsampling | Syntax-checked |
| `analysis/scripts/03_build_marker_db.md` | Marker reference curation and benchmark protocol | Protocol |
| `analysis/scripts/04_quantify_markers.sh`, `aggregate_markers.py` | DIAMOND + SingleM genome-equivalent normalisation | Syntax-checked |
| `analysis/scripts/random_effects.py` | Two-stage REML/HKSJ meta-analysis (P1) | **Self-test passed on simulated data** (`analysis/selftest_log.txt`) |
| `analysis/scripts/05_meta_analysis.R` | brms phylogenetic GLMM, metafor, proportionality, sensitivity, validation | Written; not run (no R in authoring env) |
| `figures/Fig1_concept_design.*` | Conceptual design schematic (no data) | Done |
| `manuscript/manuscript.md` | Full manuscript skeleton; Methods complete; Results empty by design | Draft |

## Run order
```bash
conda env create -f analysis/env/environment.yml && conda activate hydrosink
python analysis/scripts/01_build_inventory.py --seeds datasets.csv --config analysis/config/config.yaml \
       --metalog metalog_animal.tsv --out analysis/inventory.tsv --prisma analysis/prisma_log.tsv
# build marker DB per 03_build_marker_db.md, then run the positive control on PRJNA1068803 first
for r in $(cut -f1 analysis/inventory.tsv | tail -n +2); do
  bash analysis/scripts/02_fetch_subsample.sh $r work/$r
  bash analysis/scripts/04_quantify_markers.sh $r work/$r/sub db/hydrosink_markers work/$r
done
# concatenate work/*/*.markers.tsv -> analysis/markers_merged.tsv
Rscript analysis/scripts/05_meta_analysis.R
python analysis/scripts/random_effects.py --input analysis/p1_pairs.tsv --out analysis/results/p1_python.json  # cross-check
```
