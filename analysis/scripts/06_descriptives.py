#!/usr/bin/env python3
"""Descriptive tables from executed outputs (no modelling).

Outputs (analysis/results/):
  prisma_final.tsv        PRISMA flow from inventory to analysed sets (final state of the run logs)
  table1_p2.tsv           P2 samples / species / studies by fermentation site and vertebrate class, depth, detection
  table1_p1.tsv           P1 animals and samples by study, species and compartment
  detection_by_site.tsv   share of P2 samples with >= 2 reads per marker, by fermentation site
"""
import pathlib
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]; A = ROOT / "analysis"; R = A / "results"
p2q = [l.strip() for l in open(A / "inventory/p2_queue.txt") if l.strip()]
p1q = [l.strip() for l in open(A / "inventory/p1_queue.txt") if l.strip()]
p1 = pd.read_csv(R / "p1_matrix.tsv", sep="\t"); p2 = pd.read_csv(R / "p2_matrix.tsv", sep="\t")
aset = pd.read_csv(R / "p2_analysis_set.tsv", sep="\t") if (R / "p2_analysis_set.tsv").exists() else None
fails = pd.read_csv(A / "logs/failures.tsv", sep="\t", header=None, usecols=[0, 1], names=["run", "reason"])
nq = set(p2q) - set(p2.run_accession)
reasons = fails[fails.run.isin(nq)].drop_duplicates("run", keep="last").reason.value_counts().to_dict()

rows = [("P2 runs selected (inventory, 10 per species cap)", len(p2q), ""),
        ("P2 runs quantified (marker table produced)", len(p2), f"not quantified: {len(nq)} {reasons}"),
        ("P2 runs passing QC (rpoB reads >= 200)", int(p2.qc_pass.sum()), f"failed QC: {int((~p2.qc_pass).sum())}"),
        ("P2 runs in phylogenetic models (host has tree tip)", len(aset) if aset is not None else "", "genus-only hosts without a tip excluded"),
        ("P2 species / studies in models", f"{aset.species.nunique()} / {aset.study_accession.nunique()}" if aset is not None else "", ""),
        ("P1 runs selected", len(p1q), ""), ("P1 runs quantified", len(p1), ""),
        ("P1 runs passing QC", int(p1.qc_pass.sum()), ""),
        ("P1 animals with both chambers (QC pass)",
         int(p1[p1.qc_pass].groupby(["study", "animal"]).chamber.nunique().eq(2).sum()), "")]
pd.DataFrame(rows, columns=["stage", "n", "note"]).to_csv(R / "prisma_final.tsv", sep="\t", index=False)

d = aset if aset is not None else p2[p2.qc_pass]
g = d.groupby(["fermentation_site", "vclass"])
t = pd.DataFrame({"samples": g.size(), "species": g.species.nunique(), "studies": g.study_accession.nunique(),
                  "median_qc_reads_M": (g.qc_reads.median() / 1e6).round(2), "median_rpoB_reads": g.rpoB_n.median(),
                  "wild": g.captivity.apply(lambda s: (s == "wild").sum()),
                  "captive_domestic_lab": g.captivity.apply(lambda s: s.isin(["captive", "domestic", "laboratory"]).sum()),
                  "captivity_unknown": g.captivity.apply(lambda s: (s == "unknown").sum())}).reset_index()
tot = pd.DataFrame([{"fermentation_site": "ALL", "vclass": "ALL", "samples": len(d), "species": d.species.nunique(),
                     "studies": d.study_accession.nunique(), "median_qc_reads_M": round(d.qc_reads.median() / 1e6, 2),
                     "median_rpoB_reads": d.rpoB_n.median(), "wild": (d.captivity == "wild").sum(),
                     "captive_domestic_lab": d.captivity.isin(["captive", "domestic", "laboratory"]).sum(),
                     "captivity_unknown": (d.captivity == "unknown").sum()}])
pd.concat([t, tot]).to_csv(R / "table1_p2.tsv", sep="\t", index=False)

q1 = p1[p1.qc_pass]
t1 = q1.groupby(["study", "host_species", "compartment"]).agg(samples=("run_accession", "size"), animals=("animal", "nunique")).reset_index()
t1.to_csv(R / "table1_p1.tsv", sep="\t", index=False)

genes = ["mcrA", "acsB", "dsrA", "asrA", "nrfA", "narG", "napA", "frdA", "NiFe", "FeFe"]
det = d.groupby("fermentation_site")[[f"{x}_n" for x in genes]].apply(lambda x: (x >= 2).mean().round(3))
det.columns = genes; det.insert(0, "n", d.groupby("fermentation_site").size())
det.to_csv(R / "detection_by_site.tsv", sep="\t")
print(pd.read_csv(R / "prisma_final.tsv", sep="\t").to_string(index=False)); print(det.to_string())
