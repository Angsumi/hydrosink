#!/usr/bin/env bash
# Stream the first N forward reads of one ENA run, QC, quantify H2-sink markers, delete reads.
# Usage: 02_fetch_quantify.sh RUN_ACCESSION [N_READS]
# Env: HS_ROOT (project root), THREADS (default 4), HS_DB (DIAMOND db prefix)
# Rationale: bandwidth (~4 MB/s) forbids full downloads; Greening-lab protocol uses forward reads
# only, so only R1 is streamed. Head-of-file sampling is validated against random subsampling
# in Gate 1b (analysis/results/gate1b_head_vs_random.tsv).
set -uo pipefail
RUN=$1; N=${2:-5000000}
ROOT=${HS_ROOT:-"$(cd "$(dirname "$0")/../.." && pwd)"}
ENV=/home/angsuman/extra_spac/software_libs/conda_envs/hydrosink/bin
DB=${HS_DB:-$ROOT/analysis/db/hydrosink}
T=${THREADS:-4}
WORK=$ROOT/analysis/work/$RUN; OUT=$ROOT/analysis/results/per_sample; LOG=$ROOT/analysis/logs
mkdir -p "$OUT" "$LOG" "$OUT/../hits"
[ -s "$OUT/$RUN.markers.tsv" ] && { echo "$RUN already done"; exit 0; }
rm -rf "$WORK"; mkdir -p "$WORK"

URL=$(curl -s -m 60 --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp&format=tsv" | tail -n1 | cut -f2 | tr ';' '\n' | grep -E '(_1\.fastq\.gz|[^_][0-9]\.fastq\.gz)$' | head -n1)
[ -z "$URL" ] && { echo -e "$RUN\tNO_URL\t$(date -Is)" >> "$LOG/failures.tsv"; exit 1; }

# stream: curl gets SIGPIPE when head closes -> expected; check line count instead
# complete = exactly N reads; after 4 failed attempts the longest stream is kept if >= 2 M reads (flagged)
BEST=0
for attempt in 1 2 3 4; do
  curl -s --retry 5 --speed-limit 2000 --speed-time 60 "https://$URL" 2>/dev/null | "$ENV/pigz" -dc 2>/dev/null | head -n $((4*N)) > "$WORK/try.fq"
  L=$(wc -l < "$WORK/try.fq"); L=$((L - L % 4))
  if [ "$L" -gt "$BEST" ]; then head -n "$L" "$WORK/try.fq" > "$WORK/r1.fq"; BEST=$L; fi
  rm -f "$WORK/try.fq"
  [ "$BEST" -eq $((4*N)) ] && break
  sleep 30
done
L=$BEST
[ "$L" -lt $((4*N)) ] && [ "$L" -ge 8000000 ] && echo -e "$RUN\tPARTIAL_$((L/4))\t$(date -Is)" >> "$LOG/flags.tsv"
if [ "$L" -lt 8000000 ]; then
  echo -e "$RUN\tSHORT_STREAM_$L\t$(date -Is)" >> "$LOG/failures.tsv"; rm -rf "$WORK"; exit 1; fi

"$ENV/fastp" -i "$WORK/r1.fq" -o "$WORK/r1.qc.fq" -q 20 -l 100 --trim_poly_g --cut_tail \
   -w "$T" -j "$OUT/$RUN.fastp.json" -h /dev/null 2> "$WORK/fastp.log"
"$ENV/diamond" blastx -q "$WORK/r1.qc.fq" -d "$DB" -o "$WORK/hits.m8" --fast \
   --id 40 --query-cover 80 --evalue 1e-5 --max-target-seqs 1 --max-hsps 1 -p "$T" \
   --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore --quiet
"$ENV/python3.11" -I "$ROOT/analysis/scripts/aggregate_markers.py" --m8 "$WORK/hits.m8" \
   --labels "$DB.labels.tsv" --fastp "$OUT/$RUN.fastp.json" --run "$RUN" --n_streamed $((L/4)) \
   --out "$OUT/$RUN.markers.tsv"
rc=$?
gzip -c "$WORK/hits.m8" > "$OUT/../hits/$RUN.m8.gz"
[ $rc -eq 0 ] && echo -e "$RUN\t$(date -Is)\t$((L/4))" >> "$LOG/done.tsv"
rm -rf "$WORK"
exit $rc
