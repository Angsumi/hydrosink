#!/usr/bin/env bash
# Download, QC and subsample one run. Usage: 02_fetch_subsample.sh RUN_ACCESSION OUTDIR [PAIRS]
# Requires sra-tools, fastp, seqkit and network access to SRA. Not executed in the authoring environment.
set -euo pipefail
RUN=$1; OUT=$2; PAIRS=${3:-5000000}; SEED=20261115; THREADS=${THREADS:-8}
mkdir -p "$OUT/raw" "$OUT/qc" "$OUT/sub"
prefetch --max-size 100G -O "$OUT/raw" "$RUN"
fasterq-dump --split-files -e "$THREADS" -O "$OUT/raw" "$OUT/raw/$RUN"
fastp -i "$OUT/raw/${RUN}_1.fastq" -I "$OUT/raw/${RUN}_2.fastq" \
      -o "$OUT/qc/${RUN}_1.fq.gz" -O "$OUT/qc/${RUN}_2.fq.gz" \
      --detect_adapter_for_pe -q 20 -l 60 --dedup -w "$THREADS" \
      -j "$OUT/qc/${RUN}.fastp.json" -h "$OUT/qc/${RUN}.fastp.html"
# Paired subsampling with an identical seed keeps mates synchronised
seqkit sample -s "$SEED" -n "$PAIRS" -2 "$OUT/qc/${RUN}_1.fq.gz" -o "$OUT/sub/${RUN}_1.fq.gz"
seqkit sample -s "$SEED" -n "$PAIRS" -2 "$OUT/qc/${RUN}_2.fq.gz" -o "$OUT/sub/${RUN}_2.fq.gz"
rm -f "$OUT/raw/${RUN}"_*.fastq
echo -e "$RUN\t$(date -Is)\tok" >> "$OUT/fetch_log.tsv"
