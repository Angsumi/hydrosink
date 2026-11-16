#!/usr/bin/env bash
# Wait for the main queue to exit, then re-run it (resumable: done runs are skipped) up to 2 more passes.
cd "$(dirname "$0")/../.."
while pgrep -f "run_queue.sh analysis/inventory/main_queue.txt" >/dev/null; do sleep 300; done
for pass in 1 2; do
  echo "retry pass $pass $(date -Is)"
  bash analysis/scripts/run_queue.sh analysis/inventory/main_queue.txt 8 5000000 2
done
