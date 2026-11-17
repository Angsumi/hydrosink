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
# Fast paths (no NAME column; fastq-dump spends minutes per process fetching original spot names remotely;
# bases/qualities are identical). The 10 blocks are split into HS_VDB_PAR groups extracted in parallel
# (vdb-dump | vdb_tab_to_r1fq.py per group, each checking its own rows), then concatenated in block order.
#   HS_FETCH=vdb : vdb-dump reads the run remotely from the SRA copy on AWS (default; used on this PC).
#   HS_FETCH=sra : download the whole .sra file from AWS once (size checked), extract locally, delete it
#                  (cloud VMs: ~60 MB/s download; same SRA object, so the same reads).
EXP=0; while read -r k a b; do EXP=$((EXP + b - a + 1)); done < "$TMP/plan"
FETCH=${HS_FETCH:-vdb}
if [ "$FETCH" = vdb ] || [ "$FETCH" = sra ]; then
  PY=${HS_PY:-/home/angsuman/extra_spac/software_libs/conda_envs/hydrosink/bin/python3.11}
  CONV="$(cd "$(dirname "$0")" && pwd)/vdb_tab_to_r1fq.py"
  VP=${HS_VDB_PAR:-$([ "$FETCH" = sra ] && echo 5 || echo 2)}  # remote vdb-dump ~0.5 GB RAM each: 2 on this PC
  for t in 1 2 3; do
    SRC="$RUN"; ok=1
    if [ "$FETCH" = sra ]; then
      SRC="$TMP/$RUN.sra"; U="https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN"
      CL=$(curl -sI --max-time 60 "$U" | awk -F': ' 'tolower($1)=="content-length"{print $2+0}')
      if [ -z "$CL" ] || [ "$CL" -gt $(( ${HS_SRA_MAX_GB:-30} * 1073741824 )) ]; then
        SRC="$RUN"; echo "no size or > ${HS_SRA_MAX_GB:-30} GB: reading remotely instead"   # avoid huge downloads
      else
        # disk reservation: wait until free space minus other runs' reservations covers this file + 5 GB
        NEED=$(( CL / 1048576 + 5120 ))
        while :; do
          exec {rl}>"$LOCKS/disk.lock"; flock "$rl"
          FREE=$(df -Pm "$TMP" | awk 'NR==2{print $4}')
          RES=$(cat "$LOCKS"/res.* 2>/dev/null | awk '{s+=$1}END{print s+0}')
          if [ $(( FREE - RES )) -ge "$NEED" ]; then echo "$NEED" > "$LOCKS/res.$RUN"; exec {rl}>&-; break; fi
          exec {rl}>&-; sleep 30
        done
        curl -sf --max-time 3600 -o "$SRC" "$U" || ok=0
        [ "$ok" = 1 ] && [ "$(stat -c %s "$SRC" 2>/dev/null || echo 0)" -eq "$CL" ] || ok=0
      fi
    fi
    if [ "$ok" = 1 ]; then
      awk -v n="$VP" -v d="$TMP" '{print > (d "/plan.g" ((NR-1) % n))}' "$TMP/plan"
      pids=()
      for g in $(seq 0 $((VP - 1))); do
        [ -s "$TMP/plan.g$g" ] || continue
        R=$(awk '{printf "%s%s-%s", (NR>1?",":""), $2, $3}' "$TMP/plan.g$g")
        ( timeout 3600 "$SRA/vdb-dump" -I -R "$R" -C READ,QUALITY,READ_LEN -f tab "$SRC" 2> "$TMP/vdb.g$g.err" \
            | "$PY" -I "$CONV" "$RUN" "$TMP/plan.g$g" > "$TMP/g$g.fq" 2> "$TMP/conv.g$g.err"
          st=("${PIPESTATUS[@]}"); [ "${st[0]}" -eq 0 ] && [ "${st[1]}" -eq 0 ] ) &
        pids+=($!)
      done
      for p in "${pids[@]}"; do wait "$p" || ok=0; done
    fi
    if [ "$ok" = 1 ]; then
      # reassemble in block order: block k is line k of plan, assigned to group (k mod VP), in order within group
      python3 -I - "$TMP" "$VP" "$OUT" <<'EOF' || ok=0
import sys, itertools
tmp, vp, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
plan = [tuple(map(int, l.split())) for l in open(f"{tmp}/plan")]
fh = {g: open(f"{tmp}/g{g}.fq") for g in range(min(vp, len(plan)))}
with open(out, "w") as o:
    for i, (k, a, b) in enumerate(plan):
        f = fh[i % vp]
        for _ in range(4 * (b - a + 1)):
            l = f.readline()
            if not l:
                sys.exit("short group file")
            o.write(l)
for f in fh.values():
    if f.readline():
        sys.exit("extra lines in group file")
EOF
    fi
    rm -f "$TMP/$RUN.sra" "$TMP"/g*.fq "$LOCKS/res.$RUN"
    if [ "$ok" = 1 ] && [ "$(wc -l < "$OUT")" -eq $((4 * EXP)) ]; then
      rm -rf "$TMP"; echo "$RUN $MODE $EXP reads ($FETCH, $VP parallel groups)"; exit 0
    fi
    echo "$FETCH attempt $t failed: $(cat "$TMP"/conv.g*.err 2>/dev/null | tail -c 200)"; sleep $((20 * t))
  done
  echo "$FETCH failed 3x; falling back to per-block fastq-dump"
fi
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
