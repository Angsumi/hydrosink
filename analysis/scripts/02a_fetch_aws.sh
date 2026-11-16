#!/usr/bin/env bash
# AWS source for one run (used on Google Colab VMs). Downloads the run's .sra object from the NCBI
# SRA Open Data bucket (s3://sra-pub-run-odp, public, no credentials), extracts the FIRST N spots'
# forward reads (fastq-dump -X N --split-files -> _1), then hands over to 02_fetch_quantify.sh,
# which skips its own ENA download when r1.ok exists and runs the identical fastp/DIAMOND/aggregate
# steps. Reads used are the first N R1 reads in SRA spot order, the same order ENA FASTQs follow.
# Usage: 02a_fetch_aws.sh RUN [N]    Env: HS_ENV (tool bin dir), THREADS
set -uo pipefail
RUN=$1; N=${2:-5000000}
ROOT=${HS_ROOT:-"$(cd "$(dirname "$0")/../.." && pwd)"}
ENV=${HS_ENV:?set HS_ENV}
OUT=$ROOT/analysis/results/per_sample; LOG=$ROOT/analysis/logs; WORK=$ROOT/analysis/work/$RUN
mkdir -p "$OUT" "$LOG" "$WORK"
[ -s "$OUT/$RUN.markers.tsv" ] && { echo "$RUN already done"; exit 0; }
if [ ! -f "$WORK/r1.ok" ]; then
  "$ENV/aws" s3 cp --no-sign-request --only-show-errors "s3://sra-pub-run-odp/sra/$RUN/$RUN" "$WORK/$RUN.sra" \
    || { echo -e "$RUN\tNO_AWS\t$(date -Is)" >> "$LOG/failures.tsv"; exit 1; }
  "$ENV/fastq-dump" -X "$N" --split-files -O "$WORK" "$WORK/$RUN.sra" > "$WORK/dump.log" 2>&1
  rm -f "$WORK/$RUN.sra" "$WORK/${RUN}_2.fastq"
  [ -s "$WORK/${RUN}_1.fastq" ] || { echo -e "$RUN\tAWS_DUMP_FAIL\t$(date -Is)" >> "$LOG/failures.tsv"; exit 1; }
  mv "$WORK/${RUN}_1.fastq" "$WORK/r1.fq"; touch "$WORK/r1.ok"
  echo -e "$RUN\tAWS_SRA_ODP\t$(date -Is)" >> "$LOG/source.tsv"
fi
exec bash "$ROOT/analysis/scripts/02_fetch_quantify.sh" "$RUN" "$N"
