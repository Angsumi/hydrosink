#!/usr/bin/env bash
# HydroSink — regenerate every table and figure from the per-sample outputs (analysis/results/per_sample, hits).
# Re-downloading/quantifying reads is NOT done here (see README "Full pipeline"); this reproduces all statistics.
# Runtime: ~2-3 h on 16 cores (MCMCglmm chains dominate). Requires the hydrosink environment (analysis/env).
set -euo pipefail
cd "$(dirname "$0")"
E=${HS_ENV:-/home/angsuman/extra_spac/software_libs/conda_envs/hydrosink/bin}
PY="$E/python3 -I"; RS="$E/Rscript"
sha256sum --quiet -c analysis/results/p1_manifest.sha256                 # integrity of P1 inputs
(cd analysis/results && sha256sum --quiet -c p2_manifest.sha256)          # integrity of P2 inputs
$PY analysis/scripts/05a_assemble.py                                   # p1_matrix.tsv, p2_matrix.tsv
$RS analysis/scripts/05_meta_analysis.R all                            # pre-registered P1, P2, P3
$PY analysis/scripts/05e_p1_offsets_compare.py                         # P1 head vs multi-offset (5 runs)
for s in -10 10; do $PY analysis/scripts/05c_identity_sensitivity.py $s; done
$PY analysis/scripts/05g_posthoc_no_mouse.py
for d in sens_id_m10 sens_id_p10 posthoc_no_mouse; do HS_RESULTS=analysis/results/$d $RS analysis/scripts/05d_sensitivity_sets.R; done
$PY analysis/scripts/08a_acsb_ref_counts.py && $PY analysis/scripts/08b_acsb_restricted_sets.py && $PY analysis/scripts/08c_acsb_composition.py
for d in acsb_bacterial acsb_bacillota; do HS_RESULTS=analysis/results/$d $RS analysis/scripts/05d_sensitivity_sets.R; done
$RS analysis/scripts/05f_dated_tree.R                                  # needs network (TimeTree API)
$PY analysis/scripts/06_descriptives.py
$PY analysis/scripts/07_figures.py
echo "done: analysis/results/*.tsv, figures/*"
# In-silico benchmark (separate; needs genome downloads, ~30 min on 56 vCPU): analysis/bench/{bench_fetch.py,bench_all.sh}
