#!/usr/bin/env python3
"""Convert `vdb-dump -I -C READ,QUALITY,READ_LEN -f tab` output (stdin) to FASTQ of read 1 (stdout).

Read 1 = first read segment of the spot (offset 0, length READ_LEN[0]) - the same segment fastq-dump
--split-files writes to <RUN>_1.fastq. Bases and Phred+33 qualities are identical to fastq-dump; the
defline uses the SRA row number instead of the original spot name (the NAME column is not fetched,
because reading it remotely is what made fastq-dump slow).
Usage: vdb_tab_to_r1fq.py RUN PLANFILE  (PLANFILE lines: k a b) - checks every planned row appears once, in order.
"""
import sys

run, plan = sys.argv[1], sys.argv[2]
rows = []
for line in open(plan):
    k, a, b = map(int, line.split())
    rows.append((a, b))
Q = {str(i): chr(i + 33) for i in range(94)}
out = sys.stdout
want_iter = (r for a, b in rows for r in range(a, b + 1))
n = 0
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    if len(f) != 4:
        sys.exit(f"bad line {n + 1}: {len(f)} fields")
    rid = int(f[0])
    exp = next(want_iter, None)
    if rid != exp:
        sys.exit(f"row {rid} where {exp} expected")
    l = int(f[3].split(", ")[0])
    if l <= 0:
        sys.exit(f"row {rid}: empty read 1")
    seq = f[1][:l]
    qual = "".join(map(Q.__getitem__, f[2].split(", ")[:l]))
    if len(seq) != l or len(qual) != l:
        sys.exit(f"row {rid}: length mismatch")
    hdr = f"{run}.{rid} {rid} length={l}"
    out.write(f"@{hdr}\n{seq}\n+{hdr}\n{qual}\n")
    n += 1
if next(want_iter, None) is not None:
    sys.exit(f"only {n} rows received")
sys.stderr.write(f"{run}: {n} reads\n")
