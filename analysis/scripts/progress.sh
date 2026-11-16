#!/usr/bin/env bash
# Quick progress report for the quantification queue.
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"; I=$ROOT/analysis/inventory; P=$ROOT/analysis/results/per_sample
cnt() { local n=0; while read -r r; do [ -s "$P/$r.markers.tsv" ] && n=$((n+1)); done < <(tail -n +2 "$1" | cut -f1); echo "$n/$(($(wc -l < "$1")-1))"; }
echo "time     : $(date '+%F %T')"
echo "Li gate1a: $(cnt "$I/li2025_runs.tsv")"
echo "P1       : $(cnt "$I/p1_runs.tsv")"
echo "P2       : $(cnt "$I/p2_runs.tsv")"
echo "failures : $(cut -f1 "$ROOT/analysis/logs/failures.tsv" 2>/dev/null | sort -u | wc -l)  partial-flags: $(wc -l < "$ROOT/analysis/logs/flags.tsv" 2>/dev/null || echo 0)"
echo "running  : $(pgrep -fc 'diamond blastx') diamond, $(ls "$ROOT/analysis/work" 2>/dev/null | wc -l) work dirs"
since=$(date -d '-1 hour' -Is); echo "last hour: $(awk -F'\t' -v s="$since" '$2>=s' "$ROOT/analysis/logs/done.tsv" | wc -l) samples"
