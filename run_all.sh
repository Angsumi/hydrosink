#!/usr/bin/env bash
# HydroSink — regenerate every table and figure from the per-sample outputs (analysis/results/per_sample, hits).
# Re-downloading/quantifying reads is NOT done here (see README "Full pipeline"); this reproduces all statistics.
# Independent steps run in parallel (each MCMCglmm chain is single-threaded); ~1 h on 16 cores, critical path = P2,
# plus ~1.5 h for the post-hoc four-chain P2 convergence check (05h; set HS_SKIP_CONVERGENCE=1 to skip it).
# P1, P2 and P3 run as separate R sessions, as in the original execution (P3 bootstrap draws depend on RNG state).
set -euo pipefail
cd "$(dirname "$0")"
E=${HS_ENV:-/home/angsuman/extra_spac/software_libs/conda_envs/hydrosink/bin}
PY="$E/python3 -I"; RS="$E/Rscript"; L=analysis/logs/stats/run_all; mkdir -p "$L"
sha256sum --quiet -c analysis/results/p1_manifest.sha256                 # integrity of P1 inputs
(cd analysis/results && sha256sum --quiet -c p2_manifest.sha256)          # integrity of P2 inputs
$PY analysis/scripts/05a_assemble.py                                      # p1_matrix.tsv, p2_matrix.tsv
pids=()
run() { local name=$1; shift; ( "$@" ) > "$L/$name.log" 2>&1 & pids+=($!); echo "started $name"; }
run p1    $RS analysis/scripts/05_meta_analysis.R p1
run p2    $RS analysis/scripts/05_meta_analysis.R p2
run p3    $RS analysis/scripts/05_meta_analysis.R p3
run offs  $PY analysis/scripts/05e_p1_offsets_compare.py
run idm   bash -c "$PY analysis/scripts/05c_identity_sensitivity.py -10 && HS_RESULTS=analysis/results/sens_id_m10 $RS analysis/scripts/05d_sensitivity_sets.R"
run idp   bash -c "$PY analysis/scripts/05c_identity_sensitivity.py 10 && HS_RESULTS=analysis/results/sens_id_p10 $RS analysis/scripts/05d_sensitivity_sets.R"
run mouse bash -c "$PY analysis/scripts/05g_posthoc_no_mouse.py && HS_RESULTS=analysis/results/posthoc_no_mouse $RS analysis/scripts/05d_sensitivity_sets.R"
run acsb  bash -c "$PY analysis/scripts/08a_acsb_ref_counts.py && $PY analysis/scripts/08b_acsb_restricted_sets.py && \
  { HS_RESULTS=analysis/results/acsb_bacterial $RS analysis/scripts/05d_sensitivity_sets.R > $L/acsb_bacterial.log 2>&1 & \
    HS_RESULTS=analysis/results/acsb_bacillota $RS analysis/scripts/05d_sensitivity_sets.R > $L/acsb_bacillota.log 2>&1 & wait; }"
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
[ $fail -eq 0 ] || { echo "a parallel step failed: see $L/"; exit 1; }
$RS analysis/scripts/05f_dated_tree.R > "$L/dated_tree.log" 2>&1        # needs p2_analysis_set.tsv (P2); TimeTree answers cached
$PY analysis/scripts/08c_acsb_composition.py                            # needs p2_analysis_set.tsv
cpids=()
if [ "${HS_SKIP_CONVERGENCE:-0}" != 1 ]; then                           # post hoc: 4 chains per P2 model, 802 and 782 sets
  for s in "analysis/results n802" "analysis/results/posthoc_no_mouse n782"; do set -- $s
    HS_RESULTS=$1 HS_LABEL=$2 $RS analysis/scripts/05h_convergence.R > "$L/convergence_$2.log" 2>&1
  done
  for s in "analysis/results n802" "analysis/results/posthoc_no_mouse n782"; do set -- $s   # detection model, 10x longer chains
    HS_RESULTS=$1 HS_LABEL=$2_detect_long HS_JOBS='mcrA_present|full' HS_NITT=1300000 HS_CORES=8 \
      $RS analysis/scripts/05h_convergence.R > "$L/convergence_$2_detect_long.log" 2>&1 & cpids+=($!)
  done
  for p in "${cpids[@]}"; do wait "$p" || { echo "convergence step failed: see $L/"; exit 1; }; done
fi
$PY analysis/scripts/06_descriptives.py
$PY analysis/scripts/07_figures.py
echo "done: analysis/results/*.tsv, figures/*"
# In-silico benchmark (separate; needs genome downloads, ~30 min on 56 vCPU): analysis/bench/{bench_fetch.py,bench_all.sh}
