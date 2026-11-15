#!/usr/bin/env bash
# Run 02_fetch_quantify.sh over a list of run accessions, P samples in parallel.
# Resumable: finished runs (non-empty per_sample/<RUN>.markers.tsv) are skipped.
# Usage: run_queue.sh RUNLIST.txt [P=4] [N_READS=5000000] [THREADS_PER_JOB=4]
set -uo pipefail
LIST=$1; P=${2:-4}; N=${3:-5000000}; export THREADS=${4:-4}
HERE="$(cd "$(dirname "$0")" && pwd)"
echo "[$(date -Is)] queue start: $(wc -l < "$LIST") runs, P=$P, N=$N, threads/job=$THREADS"
grep -v '^\s*$' "$LIST" | xargs -P "$P" -I{} bash -c '"$0" {} "$1" > /dev/null 2>&1 || echo "fail {}"' "$HERE/02_fetch_quantify.sh" "$N"
echo "[$(date -Is)] queue end"
