#!/usr/bin/env bash
# Second runner (same protocol as run.sh) for the half of the not-yet-started runs; separate output tree root2.
D="$(cd "$(dirname "$0")" && pwd)"
export HS_ROOT="$D/root2" HS_SAMPLING=head HS_SOURCE=aws HS_FETCH=vdb HS_DMND_SLOTS=2 HS_DMND_THREADS=4 HS_DMND_BLOCK=0.4 HS_LOCK_DIR="$D/root2/analysis/locks"
mkdir -p "$HS_LOCK_DIR"
cat "$D/runs2.txt" | xargs -P 3 -I{} bash "$D/root2/analysis/scripts/02_fetch_quantify.sh" {} 5000000
