#!/usr/bin/env python3
"""Gate 1a positive control: reproduce the Li et al. 2025 (Microbiome 13:28) direction —
acsB enriched in goat cecum vs rumen (~12x reported), mcrA enriched in rumen.
Pass criteria (pre-specified in research_plan.md): acsB cecum>rumen in >= 12/16 goats and
paired Wilcoxon p < 0.05; mcrA rumen>cecum median direction.
Output: analysis/results/gate1a_positive_control.tsv, figures/S1_positive_control.(png|pdf)
"""
import pandas as pd, numpy as np, pathlib, glob
from scipy.stats import wilcoxon
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[2]
li = pd.read_csv(ROOT / "analysis/inventory/li2025_runs.tsv", sep="\t")
li["comp"] = np.where(li.sample_alias.str.startswith("Rum"), "rumen", "cecum")
li["goat"] = li.sample_alias.str.replace(r"^(Rum|Ce)", "", regex=True)
rows = []
for r in li.itertuples():
    f = ROOT / f"analysis/results/per_sample/{r.run_accession}.markers.tsv"
    if not f.exists():
        continue
    m = pd.read_csv(f, sep="\t")
    g = m[m.level == "gene"].set_index("feature")
    get = lambda k, c: float(g.loc[k, c]) if k in g.index else 0.0
    rows.append(dict(run=r.run_accession, goat=r.goat, comp=r.comp, rpoB_reads=get("rpoB", "reads"),
                     acsB_reads=get("acsB", "reads"), mcrA_reads=get("mcrA", "reads"),
                     acsB=get("acsB", "cpg"), mcrA=get("mcrA", "cpg"), FeFe=get("FeFe", "cpg"), NiFe=get("NiFe", "cpg")))
d = pd.DataFrame(rows)
d.to_csv(ROOT / "analysis/results/gate1a_per_sample.tsv", sep="\t", index=False)
w = d.pivot(index="goat", columns="comp", values=["acsB", "mcrA"]).dropna()
out = []
for gene in ["acsB", "mcrA"]:
    c, ru = w[(gene, "cecum")], w[(gene, "rumen")]
    pc = 0.5 * min(x for x in pd.concat([c, ru]) if x > 0)
    lfc = np.log2((c + pc) / (ru + pc))
    p = wilcoxon(lfc).pvalue if len(lfc) >= 5 else np.nan
    out.append(dict(gene=gene, n_pairs=len(lfc), n_cecum_higher=int((lfc > 0).sum()),
                    median_fold_cecum_over_rumen=float(2 ** np.median(lfc)),
                    ratio_of_means=float(c.mean() / ru.mean()) if ru.mean() > 0 else np.inf,
                    wilcoxon_p=p))
res = pd.DataFrame(out)
res.to_csv(ROOT / "analysis/results/gate1a_positive_control.tsv", sep="\t", index=False)
print(d.groupby("comp")[["rpoB_reads", "acsB_reads", "mcrA_reads", "acsB", "mcrA", "FeFe", "NiFe"]].median().to_string())
print(res.to_string())
a = res.set_index("gene")
passed = (a.loc["acsB", "n_cecum_higher"] >= 0.75 * a.loc["acsB", "n_pairs"]) and a.loc["acsB", "wilcoxon_p"] < 0.05 \
         and a.loc["mcrA", "median_fold_cecum_over_rumen"] < 1
print("GATE 1a:", "PASS" if passed else "FAIL")

fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.2))
for i, gene in enumerate(["acsB", "mcrA"]):
    for goat, row in w.iterrows():
        ax[i].plot([0, 1], [row[(gene, "rumen")], row[(gene, "cecum")]], color="#888", lw=0.8, marker="o", ms=3)
    ax[i].set_xticks([0, 1]); ax[i].set_xticklabels(["Rumen", "Cecum"]); ax[i].set_yscale("symlog", linthresh=1e-3)
    ax[i].set_title(f"{gene} (copies per genome)", fontsize=9); ax[i].set_xlim(-0.3, 1.3)
fig.suptitle("Positive control: Li et al. 2025 goats (n=%d pairs)" % len(w), fontsize=9)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(ROOT / f"figures/S1_positive_control.{ext}", dpi=300)
