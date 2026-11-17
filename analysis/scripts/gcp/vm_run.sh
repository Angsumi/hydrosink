#!/usr/bin/env bash
# Runs on a Google Cloud worker VM: P2 protocol (amendment v1.2) on the run list given as $1
# (a slice of analysis/inventory/p2_remote_claims.txt), 3 passes (later passes retry failures; done runs skipped).
# Sizing from the VM's vCPUs: one DIAMOND search per 6 vCPUs (6 threads, -b 0.4, ~2.8 GB each), 2 extra runs
# in flight so downloads overlap with searches. Search parameters are those of 02_fetch_quantify.sh, unchanged.
LIST=${1:?run list}
cd "/home/angsuman/extra_spac/Metagenomic Paper" || exit 1
NC=$(nproc); SL=$(( NC / 6 )); [ "$SL" -lt 1 ] && SL=1
export HS_REMOTE=1 HS_SAMPLING=offsets HS_SOURCE=aws HS_FETCH=sra HS_DMND_SLOTS=$SL HS_DMND_THREADS=6 HS_DMND_BLOCK=0.4
for pass in 0 1 2; do
  echo "pass $pass $(date -Is) list=$LIST slots=$SL"
  bash analysis/scripts/run_queue.sh "$LIST" $(( SL + 4 )) 5000000 3
done
date -Is > analysis/logs/ALL_DONE
