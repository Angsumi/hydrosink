#!/usr/bin/env bash
# Re-run the 17 Gate-1a P1 runs (processed before DIAMOND hits were kept) with the identical P1 protocol (head, first
# 5 M R1 reads) in a separate output tree: regenerates hits for the identity-threshold sensitivity and serves as a
# reproducibility audit (new marker table must equal the original).
D="$(cd "$(dirname "$0")" && pwd)"
export HS_ROOT="$D/root" HS_SAMPLING=head HS_SOURCE=aws HS_FETCH=vdb HS_DMND_SLOTS=2 HS_DMND_THREADS=4 HS_DMND_BLOCK=0.4 HS_LOCK_DIR="$D/root/analysis/locks"
mkdir -p "$HS_LOCK_DIR"
cat "$D/runs.txt" | xargs -P 2 -I{} bash "$D/root/analysis/scripts/02_fetch_quantify.sh" {} 5000000
