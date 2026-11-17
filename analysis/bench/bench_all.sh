#!/usr/bin/env bash
# Benchmark driver (VM): wait for bench_fetch.py, run bench_run.sh for every genome (9 parallel x 6 threads), score.
B="$(cd "$(dirname "$0")" && pwd)"; cd "$B"
while pgrep -f "[b]ench_fetch.py" >/dev/null; do sleep 20; done
tail -n +2 ../results/bench/genome_accessions.tsv | awk -F'\t' '$3!=""{print $3" "$2}' | xargs -P 9 -n 2 bash bench_run.sh
/home/angsuman/extra_spac/software_libs/conda_envs/bench/bin/python3 -u bench_score.py > score.log 2>&1
echo BENCH_DONE >> score.log
