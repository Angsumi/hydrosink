#!/usr/bin/env bash
# Wait for the P1 queue to exit, then re-run it (resumable: done runs are skipped) up to 2 more passes.
# P1 only: P2 must not run with head sampling (amendment v1.2); P2 gets its own queue with HS_SAMPLING=offsets.
cd "$(dirname "$0")/../.."
LIST=${1:-analysis/inventory/p1_queue.txt}
while pgrep -f "run_queue.sh $LIST" >/dev/null; do sleep 300; done
for pass in 1 2; do
  echo "retry pass $pass $(date -Is)"
  HS_SAMPLING=head bash analysis/scripts/run_queue.sh "$LIST" 8 5000000 3
done
