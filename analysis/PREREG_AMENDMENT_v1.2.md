# Pre-registration amendment v1.2 (2026-11-16)

**Status at time of amendment:** 0 / 930 P2 samples quantified (no P2 marker table exists); P1: 107 / 254 per-sample
marker tables exist, not assembled or analysed. Results inspected so far: Gate 1a (prereg v1), Gate 1b (this document).
Base version: git tag `prereg-v1.1`. Nothing in the primary outcome (ALR_acet), contrasts, models, marker thresholds or
search parameters changes.

## 1. Gate 1b outcome: FAIL (pre-specified rule)
Gate 1b compared, per gene group, estimated copies per genome from the first 5 M R1 reads ("head", the v1 protocol)
with a uniform random 5 M-read subsample of the full R1 file (seed 20261116), in goat 14 caecum (SRR28204439) and
rumen (SRR28204452). Rule (`gate1b_compare.py`, committed before data): every gene group must pass;
>= 100 reads in both subsamples -> |head/random - 1| <= 10 %; otherwise a Poisson test with p >= 0.05.

| Sample | mcrA | acsB | resp. sum | [NiFe] | [FeFe] |
|---|---|---|---|---|---|
| Ce14 | -3.8 % pass | +3.2 % pass | -0.3 % pass | -8.2 % pass | -7.3 % pass |
| Rum14 | -8.5 % pass | -0.8 % pass | +5.7 % pass | **+10.3 % FAIL** (210 vs 189 reads) | -5.8 % pass |

Source: `analysis/results/gate1b_head_vs_random.tsv` (commit 58da916). **Result: FAIL (9 / 10 checks).**

Interpretation, stated without changing the result: with ~200 reads per subsample the Poisson standard error of a
ratio is about 10 % (sqrt(1/210 + 1/189) = 0.10), and the exact conditional test for 210 vs 189 (scaled by rpoB reads,
4395 vs 4362) gives p = 0.34. The 10 % criterion for genes with 100-1000 reads was therefore too strict to distinguish
bias from counting noise; this was a design error in the rule. The FAIL is **not** re-evaluated post hoc and is
reported in the manuscript as such. The pre-specified consequence (switch to multi-offset sampling) is applied below.

Technical note: the first Rum14 full-file download was size-correct but failed the gzip CRC; its partial output is
quarantined in `analysis/results/gate1b/invalid_corrupt_download/`, and the file was re-downloaded and verified
against the ENA MD5 before the result above was produced.

## 2. Sampling protocol from now on
### P2 (all 930 runs): multi-offset sampling (new)
- 10 blocks of 500,000 consecutive spots; block k (k = 0..9) starts at spot floor(k x S / 10) + 1, where S is the
  run's total spot count (ENA `read_count`). Forward reads (R1) only, as before. Deterministic (no random seed).
- If S <= 5 M, all spots are used (no P2 run is affected: minimum S in the P2 inventory is 5,002,285).
- Implementation: `analysis/scripts/fetch_aws_blocks.sh` (MODE=offsets), called by `02_fetch_quantify.sh` with
  `HS_SAMPLING=offsets`. A run is accepted only if every block returns exactly the expected number of reads.
- Runs that cannot be fetched this way are **not** replaced by head sampling; they are recorded in
  `logs/failures.tsv` (`AWS_FAIL_offsets`) and reported as "not retrievable" in the PRISMA flow.
- Verified before any P2 quantification: on P2 run SRR24915848 (S = 19,460,180) the procedure returned 5,000,000
  reads whose spot numbers matched the planned block positions exactly (0 mismatches). Only reads were fetched;
  no QC, DIAMOND or marker output was produced, and the reads were deleted.

### P1 (254 runs): unchanged (head, first 5 M R1 reads)
- P1 is the positive control / calibration (amendment v1.1), and the two genes that define the primary outcome,
  acsB and mcrA, passed Gate 1b in both samples. Changing P1 now would mix protocols within P1 (107 runs done).
- [NiFe]-hydrogenase results from P1 carry the Gate 1b caveat in the manuscript.
- **Pre-specified sensitivity (descriptive, not a gate):** 5 P1 runs, drawn now with `random.Random(20261116)` from
  the sorted P1 queue - SRR12529109, SRR12529110, SRR12529203, SRR12529313, SRR28204455 - are re-quantified with
  multi-offset sampling. Per gene group we report head vs multi-offset copies per genome (ratio, Poisson CI). These
  numbers are reported whatever they show and do not alter the P1 analysis.

## 3. Data source change (no effect on which reads are used)
- Default source is now the NCBI SRA copy on AWS (public Open Data bucket), read remotely with `fastq-dump`
  (sra-tools 3.1.1, separate environment `software_libs/conda_envs/sratools`), which transfers only the requested
  spots. Reason: after a connection upgrade, ENA delivers ~1 MB/s per file while AWS delivers ~8 MB/s.
- For head sampling the reads are identical to the ENA head: (i) block concatenation was verified to give spot
  numbers 1..5,000,000 in order (SRR12529199); (ii) a run quantified from the SRA first 5 M spots reproduced the
  ENA-based marker table exactly (55 / 55 features, identical read counts and qc_reads; SRR12529199).
- P1 falls back to the ENA head download if AWS fails (same reads; flagged `AWS_FAIL_ENA_FALLBACK`). P2 has no fallback.
- The source and sampling mode of every run are logged in `analysis/logs/source.tsv`.
- Google Colab was tested for one run (above) and is not part of the pipeline; the Colab-specific script
  (`02a_fetch_aws.sh`) is removed.

## 4. Pipeline fix recorded here for completeness (before this amendment, no effect on results)
- ENA byte-range resume (commit 560fe68): `curl --retry` could duplicate bytes after a dropped connection,
  corrupting the partial gzip; the read-count check then extended the download to the whole file. Affected runs
  were in flight only (seven partial files deleted). A completed run required all 5 M reads to decompress without
  error before the corruption point, so finished samples are unaffected.

## 5. Still planned, not implemented (unchanged from v1.1)
In-silico marker specificity benchmark; acsB phylogenetic placement; dated-tree sensitivity; Methods rewrite.

## Addendum (2026-11-17, during P2; no effect on results)
P2 started 2026-11-17 07:49. The job was stopped twice by the system memory guard (systemd-oomd) after 40 / 930 runs:
three concurrent DIAMOND searches at `-b 1` need ~5.4 GB each. DIAMOND block size changed from `-b 1` to `-b 0.4`
(`HS_DMND_BLOCK`). Test on P2 run SRR25820354 (same query, all other parameters equal): peak RAM 5.44 GB (`-b 1`) vs
2.77 GB (`-b 0.4`), wall time 7:39 for both, and **identical hit files** (31,775 hits, identical md5 of sorted output).
The reference database (99 MB) fits in one block at either setting, so each query is searched against the whole
reference in one pass. Thread count did not change memory (5.36 GB at 3 threads). P1 (run at `-b 1`) is not re-run.
The 40 P2 runs finished before the change used `-b 1`.
- Further fetch change (2026-11-17, no effect on reads used): remote `fastq-dump` spent ~4-5 min per process fetching
  the original spot-name column (11 spots took 236 s; the same rows without names took 11 s with `vdb-dump`), so each
  P2 run took ~45 min to fetch. The fetcher now makes one `vdb-dump` call per run for all 10 blocks
  (`-I -C READ,QUALITY,READ_LEN`; read 1 = first segment, converted by `vdb_tab_to_r1fq.py`, which also checks that
  every planned row arrives once, in order). Check on SRR32584038, rows 7,236,266-7,256,265 (20,000 reads): bases and
  qualities identical to `fastq-dump --split-files` read 1 (identical md5). Only the FASTQ header differs (SRA row
  number instead of the original spot name). Per-block `fastq-dump` remains as automatic fallback.
