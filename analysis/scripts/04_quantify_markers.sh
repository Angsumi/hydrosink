#!/usr/bin/env bash
# Quantify H2-sink markers per sample and normalise to genome equivalents.
# Usage: 04_quantify_markers.sh RUN SUBDIR DBPREFIX OUTDIR
# Not executed in the authoring environment.
set -euo pipefail
RUN=$1; SUB=$2; DB=$3; OUT=$4; THREADS=${THREADS:-8}
mkdir -p "$OUT/diamond" "$OUT/singlem"
# Forward read only for marker counting (avoids double counting pairs)
diamond blastx -q "$SUB/${RUN}_1.fq.gz" -d "$DB" -o "$OUT/diamond/${RUN}.m8" \
  --fast --id 40 --query-cover 80 --evalue 1e-10 --max-target-seqs 1 -p "$THREADS" \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore
# Genome equivalents: SingleM 'pipe' gives per-sample coverage of single-copy marker genes
singlem pipe -1 "$SUB/${RUN}_1.fq.gz" -2 "$SUB/${RUN}_2.fq.gz" \
  --otu-table "$OUT/singlem/${RUN}.otu.tsv" --threads "$THREADS"
python "$(dirname "$0")/aggregate_markers.py" --m8 "$OUT/diamond/${RUN}.m8" \
  --labels "${DB}.labels.tsv" --singlem "$OUT/singlem/${RUN}.otu.tsv" \
  --run "$RUN" --out "$OUT/${RUN}.markers.tsv"
