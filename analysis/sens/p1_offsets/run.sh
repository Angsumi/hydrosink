#!/usr/bin/env bash
# Pre-specified (amendment v1.2 §2): re-quantify 5 seeded P1 runs with multi-offset sampling, separate output tree.
D="$(cd "$(dirname "$0")" && pwd)"
export HS_ROOT="$D/root" HS_SAMPLING=offsets HS_SOURCE=aws HS_FETCH=vdb HS_DMND_SLOTS=2 HS_DMND_THREADS=6 HS_DMND_BLOCK=0.4 HS_LOCK_DIR="$D/root/analysis/locks"
mkdir -p "$HS_LOCK_DIR"
cat "$D/runs.txt" | xargs -P 2 -I{} bash "$D/root/analysis/scripts/02_fetch_quantify.sh" {} 5000000
