#!/usr/bin/env bash
# Quick progress report for the quantification queue.
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"; I=$ROOT/analysis/inventory; P=$ROOT/analysis/results/per_sample
cnt() { local n=0; while read -r r; do [ -s "$P/$r.markers.tsv" ] && n=$((n+1)); done < <(tail -n +2 "$1" | cut -f1); echo "$n/$(($(wc -l < "$1")-1))"; }
echo "time     : $(date '+%F %T')"
echo "Li gate1a: $(cnt "$I/li2025_runs.tsv")"
echo "P1       : $(cnt "$I/p1_runs.tsv")"
echo "P2       : $(cnt "$I/p2_runs.tsv")"
echo "failures : $(cut -f1 "$ROOT/analysis/logs/failures.tsv" 2>/dev/null | sort -u | wc -l)  partial-flags: $(cat "$ROOT/analysis/logs/flags.tsv" 2>/dev/null | wc -l)"
echo "running  : $(pgrep -fc 'diamond blastx') diamond, $(ls "$ROOT/analysis/work" 2>/dev/null | wc -l) work dirs"
since=$(date -d '-1 hour' -Is); echo "last hour: $(awk -F'\t' -v s="$since" '$2>=s' "$ROOT/analysis/logs/done.tsv" | wc -l) samples"
echo "services : $(for s in hs-queue hs-retry hs-gate1b; do printf '%s=%s ' $s "$(systemctl --user is-active $s 2>/dev/null)"; done)"
since10=$(date -d '-10 min' -Is); n10=$(awk -F'\t' -v s="$since10" '$2>=s' "$ROOT/analysis/logs/done.tsv" | wc -l)
left=$(( $(grep -c . "$ROOT/analysis/inventory/main_queue.txt") - $(ls "$P"/*.markers.tsv 2>/dev/null | wc -l) ))
rate=$(awk -F'\t' -v s="$since" '$2>=s' "$ROOT/analysis/logs/done.tsv" | wc -l)
echo "last 10m : $n10 samples   remaining: $left   ETA at last-hour rate: $( [ "$rate" -gt 0 ] && echo "$((left/rate)) h" || echo n/a)"
G=$ROOT/analysis/work/gate1b; A=$(pgrep -x aria2c | head -1)
if [ -n "$A" ]; then echo "gate1b   : downloading $(for f in "$G"/*.aria2; do b=${f##*/}; printf "%s " "${b%%.*}"; done) (written this session $(( $(awk '/^write_bytes/{print $2}' /proc/$A/io)/1000000 )) MB of ~8800)";
else echo "gate1b   : $(tail -n1 "$ROOT/analysis/logs/gate1b.log" 2>/dev/null)"; fi
