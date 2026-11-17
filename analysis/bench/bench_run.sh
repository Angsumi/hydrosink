#!/usr/bin/env bash
# Benchmark step 2 (one genome): simulate reads with wgsim, then fastp + DIAMOND exactly as the study pipeline.
# Usage: bench_run.sh ACCESSION GROUP   (run from analysis/bench; see BENCHMARK_PROTOCOL.md)
set -euo pipefail
ACC=$1; GRP=$2
SW=/home/angsuman/extra_spac/software_libs; ENV=$SW/conda_envs/hydrosink/bin; BE=$SW/conda_envs/bench/bin
B="$(cd "$(dirname "$0")" && pwd)"; DB="$B/../db/hydrosink"; W="$B/work/$ACC"; O="$B/../results/bench/runs"; mkdir -p "$W" "$O"
[ -s "$O/$ACC.m8.gz" ] && exit 0
N=500000; [ "$GRP" = host_decoy ] && N=2000000
"$BE/wgsim" -e 0.005 -r 0.001 -R 0.15 -1 150 -2 150 -d 150 -s 0 -N "$N" -S 20261117 "$B/genomes/$ACC/genome.fna" "$W/r1.fq" /dev/null > /dev/null 2> "$W/wgsim.log"
"$ENV/fastp" -i "$W/r1.fq" -o "$W/r1.qc.fq" -q 20 -l 100 --trim_poly_g --cut_tail -w 4 -j "$O/$ACC.fastp.json" -h /dev/null 2> "$W/fastp.log"
"$ENV/diamond" blastx -q "$W/r1.qc.fq" -d "$DB" -o "$W/hits.m8" --fast --id 40 --query-cover 80 --evalue 1e-5 \
  --max-target-seqs 1 --max-hsps 1 -p 6 -b 0.4 -c 4 --tmpdir "$W" --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore --quiet
gzip -c "$W/hits.m8" > "$O/$ACC.m8.gz.tmp" && mv "$O/$ACC.m8.gz.tmp" "$O/$ACC.m8.gz"
awk 'NR%4==1{print substr($1,2)}' "$W/r1.fq" | gzip -c > "$O/$ACC.names.gz"
rm -rf "$W"
