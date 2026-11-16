#!/usr/bin/env bash
# Gate 1b: is head-of-file sampling (first N R1 reads, as used in 02_fetch_quantify.sh) unbiased
# relative to a uniform random subsample of the whole R1 file?
# Runs: goat 14, cecum (SRR28204439, Ce14) and rumen (SRR28204452, Rum14) from PRJNA1068803.
# For each run: download full R1 (resumable, size-verified against ENA fastq_bytes), draw a random
# 5 M-read subsample (seqkit sample, seed 20261116), run the identical fastp/DIAMOND/aggregate steps,
# then compare against the existing head-of-file result (results/per_sample/RUN.markers.tsv).
# Pre-specified decision rule (written before any random-subsample result was seen) is in
# gate1b_compare.py. Full files are deleted after use.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ENV=/home/angsuman/extra_spac/software_libs/conda_envs/hydrosink/bin
DB=$ROOT/analysis/db/hydrosink; T=4; N=5000000; SEED=20261116
G=$ROOT/analysis/results/gate1b; W=$ROOT/analysis/work/gate1b; mkdir -p "$G" "$W"

for RUN in SRR28204439 SRR28204452; do
  [ -s "$G/$RUN.random.markers.tsv" ] && continue
  REP=$(curl -s -m 60 --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp,fastq_bytes,read_count&format=tsv" | tail -n1)
  URL=$(echo "$REP" | cut -f2 | cut -d';' -f1); SIZE=$(echo "$REP" | cut -f3 | cut -d';' -f1)
  F=$W/$RUN.R1.fastq.gz
  until [ -f "$F" ] && [ ! -f "$F.aria2" ] && [ "$(stat -c %s "$F")" -eq "$SIZE" ]; do
    aria2c -q -c --file-allocation=none -x 16 -s 16 -k 20M --max-tries=0 --retry-wait=20 --lowest-speed-limit=2K -d "$W" -o "$(basename "$F")" "https://$URL"; sleep 20
  done
  echo "$RUN downloaded $(date -Is)"
  TOTAL=$("$ENV/pigz" -dc "$F" | awk 'END{print NR/4}')
  P=$(awk -v n=$N -v t="$TOTAL" 'BEGIN{printf "%.6f", 1.05*n/t}')
  "$ENV/pigz" -dc "$F" | "$ENV/seqkit" sample -p "$P" -s $SEED -j 2 2>/dev/null | head -n $((4*N)) > "$W/$RUN.r1.fq"
  L=$(wc -l < "$W/$RUN.r1.fq")
  "$ENV/fastp" -i "$W/$RUN.r1.fq" -o "$W/$RUN.qc.fq" -q 20 -l 100 --trim_poly_g --cut_tail \
     -w $T -j "$G/$RUN.random.fastp.json" -h /dev/null 2> /dev/null
  "$ENV/diamond" blastx -q "$W/$RUN.qc.fq" -d "$DB" -o "$W/$RUN.m8" --fast \
     --id 40 --query-cover 80 --evalue 1e-5 --max-target-seqs 1 --max-hsps 1 -p $T \
     --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore --quiet
  "$ENV/python3.11" -I "$ROOT/analysis/scripts/aggregate_markers.py" --m8 "$W/$RUN.m8" \
     --labels "$DB.labels.tsv" --fastp "$G/$RUN.random.fastp.json" --run "$RUN" --n_streamed $((L/4)) \
     --out "$G/$RUN.random.markers.tsv"
  gzip -c "$W/$RUN.m8" > "$G/$RUN.random.m8.gz"
  echo -e "$RUN\ttotal_reads=$TOTAL\tsampled=$((L/4))\tseed=$SEED" >> "$G/sampling_log.tsv"
  rm -f "$F" "$W/$RUN.r1.fq" "$W/$RUN.qc.fq" "$W/$RUN.m8"
done
"$ENV/python3.11" -I "$ROOT/analysis/scripts/gate1b_compare.py"
