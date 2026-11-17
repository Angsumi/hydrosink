#!/usr/bin/env bash
# Fetch forward reads (R1) of one run from the NCBI SRA copy on AWS, in 10 blocks fetched in parallel.
# fastq-dump reads the run remotely and transfers only the requested spots (no full-file download).
# Usage: fetch_aws_blocks.sh RUN N MODE OUT.fq TOTAL_SPOTS
#   MODE=head    : spots 1..N, split into 10 consecutive blocks and concatenated in order
#                  -> identical to the first N R1 reads of the ENA FASTQ (P1 protocol, prereg v1).
#   MODE=offsets : 10 blocks of N/10 spots starting at floor(k*S/10)+1, k=0..9 (P2 protocol, amendment v1.2).
# Each block is retried up to 3 times; exits non-zero unless every block has exactly the expected reads.
set -uo pipefail
RUN=$1; N=$2; MODE=$3; OUT=$4; S=$5
SRA=${HS_SRA_ENV:-/home/angsuman/extra_spac/software_libs/conda_envs/sratools/bin}
K=10; B=$((N / K)); PAR=${HS_BLOCK_PAR:-5}
FQD=${HS_FQD_SLOTS:-8}; LOCKS=${HS_LOCK_DIR:-$(cd "$(dirname "$0")/.." && pwd)/locks}; mkdir -p "$LOCKS"
[ "$S" -gt 0 ] 2>/dev/null || { echo "need TOTAL_SPOTS"; exit 2; }
# offsets with S <= N would make blocks overlap (duplicate reads): take all S spots instead
[ "$MODE" = offsets ] && [ "$S" -le "$N" ] && MODE=head
TMP="$OUT.blocks"; rm -rf "$TMP"; mkdir -p "$TMP"
: > "$TMP/plan"
for k in $(seq 0 $((K - 1))); do
  if [ "$MODE" = head ]; then a=$((k * B + 1)); elif [ "$MODE" = offsets ]; then a=$((k * S / K + 1)); else echo "bad MODE"; exit 2; fi
  b=$((a + B - 1)); [ "$b" -gt "$S" ] && b=$S
  [ "$a" -le "$S" ] && echo "$k $a $b" >> "$TMP/plan"
done
fetch() {  # k a b
  local k=$1 a=$2 b=$3 want=$(( 4 * ($3 - $2 + 1) )) t f
  for t in 1 2 3; do
    rm -rf "$TMP/b$k"
    # machine-wide cap on concurrent fastq-dump (each caches up to ~0.7 GB RAM when reading remotely)
    while :; do for s in $(seq 1 "$FQD"); do exec {fd}>"$LOCKS/fqd.$s"; flock -n "$fd" && break 2; exec {fd}>&-; done; sleep 5; done
    timeout 1800 "$SRA/fastq-dump" -N "$a" -X "$b" --split-files -O "$TMP/b$k" "$RUN" > "$TMP/b$k.log" 2>&1
    exec {fd}>&-
    f="$TMP/b$k/${RUN}_1.fastq"; [ -f "$f" ] || f="$TMP/b$k/${RUN}.fastq"
    if [ -f "$f" ] && [ "$(wc -l < "$f")" -eq "$want" ]; then mv "$f" "$TMP/r1.$k.fq"; rm -rf "$TMP/b$k"; return 0; fi
    sleep $((10 * t))
  done
  return 1
}
while read -r k a b; do
  fetch "$k" "$a" "$b" &
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do wait -n; done
done < "$TMP/plan"
wait
EXP=0; : > "$OUT"
while read -r k a b; do
  [ -f "$TMP/r1.$k.fq" ] || { echo "$RUN block $k ($a-$b) failed"; exit 1; }
  cat "$TMP/r1.$k.fq" >> "$OUT"; EXP=$((EXP + b - a + 1))
done < "$TMP/plan"
L=$(wc -l < "$OUT")
[ "$L" -eq $((4 * EXP)) ] || { echo "$RUN got $((L / 4)) reads, expected $EXP"; exit 1; }
rm -rf "$TMP"
echo "$RUN $MODE $EXP reads"
