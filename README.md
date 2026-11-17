# HydroSink — gut hydrogen-sink genes across vertebrates (reproducibility package)

**Status (2026-11-17): data processed, all pre-specified analyses run, manuscript draft v1 in `manuscript/manuscript.md`.**
Every number in the manuscript comes from a file in `analysis/results/` named next to it.

## What is where
| Path | Contents |
|---|---|
| `research_plan.md` | Hypotheses (P1–P4), inclusion criteria, analysis plan, decision log |
| `analysis/config/config.yaml` | Pre-registered parameters (sampling, thresholds, outcomes) |
| `analysis/PREREG_AMENDMENT_v1.1.md`, `v1.2.md` | Amendments committed before P2 data existed (git tags `prereg-v1`, `-v1.1`, `-v1.2`) |
| `analysis/inventory/` | ENA screening outputs, run lists (`p1_runs.tsv`, `p2_runs.tsv`), host traits, host tree, curation |
| `analysis/db/` | Marker database (Greening-lab v3 + UniRef90 RpoB; md5 in config) |
| `analysis/results/per_sample/`, `hits/` | Per-run marker tables, fastp reports, DIAMOND best hits (SHA-256: `p1_manifest.sha256`, `p2_manifest.sha256`) |
| `analysis/results/*.tsv` | All analysis outputs (P1, P2, P3, descriptives); sensitivity sets in sub-folders |
| `analysis/bench/` | In-silico marker specificity benchmark (protocol, panel, scripts); results in `analysis/results/bench/` |
| `analysis/reporting/` | PRISMA 2020 and STORMS v1.03 checklists |
| `figures/` | Figures 1–4, S1–S2 (PDF + PNG) |
| `literature.csv` | References with verification status (only VERIFIED entries are cited) |

## Reproduce the statistics (from per-sample outputs)
```bash
bash run_all.sh
```
Checks file integrity, rebuilds the matrices and re-runs every model, sensitivity analysis, table and figure.
About 2–3 h, dominated by the MCMCglmm chains. The dated-tree step needs internet access to the TimeTree API.

## Full pipeline (from public reads)
1. `analysis/scripts/01_build_inventory.py`, `01b_host_traits.py`, `01c_host_tree.py`, `01d_curate_metadata.py`: screening, traits, tree, curation.
2. `analysis/scripts/run_queue.sh <list> <parallel> 5000000 <threads>` runs `02_fetch_quantify.sh` per run:
   - fetch: `fetch_aws_blocks.sh`; P1 used `HS_SAMPLING=head`, P2 used `HS_SAMPLING=offsets`;
   - then fastp, DIAMOND and `aggregate_markers.py`.
3. `bash run_all.sh`
4. Benchmark: `analysis/bench/bench_fetch.py`, then `bench_all.sh` (simulation, search and scoring).

Tools: DIAMOND 2.1.10, fastp 0.23.4, Python 3.11, sra-tools 3.1.1, R 4.3.3 (metafor 4.8-0, lme4 1.1-37, MCMCglmm 2.36, ape 5.8-1), wgsim, NCBI datasets.
The environment is in `analysis/env/environment.yml`; cloud VM setup is in `analysis/scripts/gcp/vm_setup.sh`.

## Honest status notes
- **Not registered** in an external registry. The git history (tags) is the time-stamped record.
- **Post-hoc analyses**, all labelled as such in the manuscript:
  - no-mouse set;
  - *acsB* taxonomic attribution, which replaced the planned read placement because reads were not retained;
  - benchmark executed after the main analyses.
- **Not done:** per-study risk-of-bias scoring and certainty-of-evidence grading (see the PRISMA checklist).
