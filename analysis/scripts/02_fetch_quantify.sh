#!/usr/bin/env bash
# Stream the first N forward reads of one ENA run, QC, quantify H2-sink markers, delete reads.
# Usage: 02_fetch_quantify.sh RUN_ACCESSION [N_READS]
# Env: HS_ROOT (project root), THREADS (default 4), HS_DB (DIAMOND db prefix)
# Rationale: bandwidth (~4 MB/s) forbids full downloads; Greening-lab protocol uses forward reads
# only, so only R1 is streamed. Head-of-file sampling is validated against random subsampling
# in Gate 1b (analysis/results/gate1b_head_vs_random.tsv).
set -uo pipefail
RUN=$1; N=${2:-5000000}; MINR=${HS_MIN_READS:-2000000}
ROOT=${HS_ROOT:-"$(cd "$(dirname "$0")/../.." && pwd)"}
ENV=${HS_ENV:-/home/angsuman/extra_spac/software_libs/conda_envs/hydrosink/bin}
DB=${HS_DB:-$ROOT/analysis/db/hydrosink}
T=${THREADS:-4}
WORK=$ROOT/analysis/work/$RUN; OUT=$ROOT/analysis/results/per_sample; LOG=$ROOT/analysis/logs
mkdir -p "$OUT" "$LOG" "$OUT/../hits"
[ -s "$OUT/$RUN.markers.tsv" ] && { echo "$RUN already done"; exit 0; }
mkdir -p "$WORK"
# Resumable download (v2, 2026-11-16): the compressed head of R1 is saved on the SSD (raw.gz) and
# extended with HTTP byte ranges, so an interrupted run continues where it stopped instead of
# restarting. Reads used are still exactly the first N R1 reads, i.e. identical to v1.
if [ ! -f "$WORK/r1.ok" ]; then
  rm -f "$WORK/try.fq" "$WORK/r1.fq"
  REP=$(curl -s -m 60 --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp,fastq_bytes,read_count&format=tsv" | tail -n1)
  read -r URL BYTES RC < <(echo "$REP" | awk -F'\t' '{n=split($2,u,";"); split($3,b,";"); for(i=1;i<=n;i++) if(u[i]~/(_1\.fastq\.gz|[^_][0-9]\.fastq\.gz)$/){print u[i], b[i], $4; exit}}')
  [ -z "${URL:-}" ] && { echo -e "$RUN\tNO_URL\t$(date -Is)" >> "$LOG/failures.tsv"; exit 1; }
  # v3 (2026-11-16): default source is the NCBI SRA copy on AWS, fetched as 10 parallel blocks
  # (fetch_aws_blocks.sh). HS_SAMPLING=head -> first N spots (P1; same reads as the ENA head, verified
  # identical marker counts on SRR12529199); HS_SAMPLING=offsets -> 10 evenly spaced blocks (P2, v1.2).
  SRC=${HS_SOURCE:-aws}; MODE=${HS_SAMPLING:-head}
  if [ "$SRC" = aws ]; then
    rm -f "$WORK/raw.gz" "$WORK/part"
    if bash "$ROOT/analysis/scripts/fetch_aws_blocks.sh" "$RUN" "$N" "$MODE" "$WORK/r1.fq" "$RC" > "$WORK/aws.log" 2>&1; then
      touch "$WORK/r1.ok"; echo -e "$RUN\tAWS_SRA\t$MODE\t$(date -Is)" >> "$LOG/source.tsv"
    elif [ "$MODE" != head ]; then
      echo -e "$RUN\tAWS_FAIL_$MODE\t$(date -Is)" >> "$LOG/failures.tsv"; exit 1
    else
      echo -e "$RUN\tAWS_FAIL_ENA_FALLBACK\t$(date -Is)" >> "$LOG/flags.tsv"; rm -f "$WORK/r1.fq"
    fi
  fi
  [ "$MODE" != head ] && [ ! -f "$WORK/r1.ok" ] && exit 1
  if [ ! -f "$WORK/r1.ok" ]; then  # ENA head download (fallback, or HS_SOURCE=ena)
  if [ "$RC" -gt "$N" ]; then TGT=$(( BYTES / RC * N * 103 / 100 + 2000000 )); else TGT=$BYTES; fi
  [ "$TGT" -gt "$BYTES" ] && TGT=$BYTES
  STALL=0
  while :; do
    HAVE=$(stat -c %s "$WORK/raw.gz" 2>/dev/null || echo 0)
    if [ "$HAVE" -lt "$TGT" ]; then
      # One attempt per call, written to a part file and appended only on HTTP 206 (no curl --retry:
      # it restarts the range after a drop and duplicated bytes, corrupting raw.gz). Retries = this loop.
      CODE=$(curl -s --speed-limit 2000 --speed-time 60 -r "$HAVE-$((TGT-1))" -o "$WORK/part" -w '%{http_code}' "https://$URL")
      [ "$CODE" = "206" ] && [ -s "$WORK/part" ] && cat "$WORK/part" >> "$WORK/raw.gz"
      rm -f "$WORK/part"
      NOW=$(stat -c %s "$WORK/raw.gz" 2>/dev/null || echo 0)
      if [ "$NOW" -lt "$TGT" ]; then [ "$NOW" -eq "$HAVE" ] && STALL=$((STALL+1)); [ "$STALL" -ge 20 ] && break; sleep 30; continue; fi
    fi
    # corrupt (not merely truncated) gzip -> discard and start over
    if "$ENV/pigz" -dc "$WORK/raw.gz" 2>&1 >/dev/null | grep -q "invalid"; then
      echo -e "$RUN\tCORRUPT_RESTART\t$(date -Is)" >> "$LOG/flags.tsv"; rm -f "$WORK/raw.gz"; continue; fi
    LN=$("$ENV/pigz" -dc "$WORK/raw.gz" 2>/dev/null | head -n $((4*N)) | wc -l)
    [ "$LN" -ge $((4*N)) ] && break
    [ "$TGT" -ge "$BYTES" ] && break
    TGT=$(( TGT * 115 / 100 )); [ "$TGT" -gt "$BYTES" ] && TGT=$BYTES
  done
  "$ENV/pigz" -dc "$WORK/raw.gz" 2>/dev/null | head -n $((4*N)) > "$WORK/r1.fq"
  L=$(wc -l < "$WORK/r1.fq"); L=$((L - L % 4)); head -n "$L" "$WORK/r1.fq" > "$WORK/r1.tmp" && mv "$WORK/r1.tmp" "$WORK/r1.fq"
  if [ "$L" -eq $((4*N)) ] || [ "$TGT" -ge "$BYTES" ]; then touch "$WORK/r1.ok"; rm -f "$WORK/raw.gz"; fi
  fi
fi
L=$(wc -l < "$WORK/r1.fq"); L=$((L - L % 4))
[ "$L" -lt $((4*N)) ] && [ "$L" -ge $((4*MINR)) ] && echo -e "$RUN\tPARTIAL_$((L/4))\t$(date -Is)" >> "$LOG/flags.tsv"
if [ "$L" -lt $((4*MINR)) ]; then
  echo -e "$RUN\tSHORT_STREAM_$L\t$(date -Is)" >> "$LOG/failures.tsv"; exit 1; fi

"$ENV/fastp" -i "$WORK/r1.fq" -o "$WORK/r1.qc.fq" -q 20 -l 100 --trim_poly_g --cut_tail \
   -w "$T" -j "$OUT/$RUN.fastp.json" -h /dev/null 2> "$WORK/fastp.log"
# RAM guard: at most HS_DMND_SLOTS (default 3) DIAMOND jobs at once; smaller blocks (-b 1 -c 4).
# Search parameters (sensitivity, thresholds) are unchanged, so results are identical.
# DIAMOND threads: HS_DMND_THREADS (default 5; 3 slots x 5 = 15 of 16 cores). Thread count does not change hits.
SLOTS=${HS_DMND_SLOTS:-3}; mkdir -p "$ROOT/analysis/locks"
while :; do for s in $(seq 1 "$SLOTS"); do exec 9>"$ROOT/analysis/locks/dmnd.$s"; flock -n 9 && break 2; done; sleep 15; done
"$ENV/diamond" blastx -q "$WORK/r1.qc.fq" -d "$DB" -o "$WORK/hits.m8" --fast \
   --id 40 --query-cover 80 --evalue 1e-5 --max-target-seqs 1 --max-hsps 1 -p "${HS_DMND_THREADS:-5}" -b 1 -c 4 \
   --tmpdir "$WORK" --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore --quiet
exec 9>&-
"$ENV/python3.11" -I "$ROOT/analysis/scripts/aggregate_markers.py" --m8 "$WORK/hits.m8" \
   --labels "$DB.labels.tsv" --fastp "$OUT/$RUN.fastp.json" --run "$RUN" --n_streamed $((L/4)) \
   --out "$OUT/$RUN.markers.tsv"
rc=$?
gzip -c "$WORK/hits.m8" > "$OUT/../hits/$RUN.m8.gz"
[ $rc -eq 0 ] && echo -e "$RUN\t$(date -Is)\t$((L/4))" >> "$LOG/done.tsv"
rm -rf "$WORK"
exit $rc
