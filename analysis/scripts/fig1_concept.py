"""Fig. 1 - conceptual design schematic (no data)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
fig, ax = plt.subplots(figsize=(10, 5.2)); ax.set_xlim(0, 10); ax.set_ylim(0, 5.2); ax.axis("off")
def box(x, y, w, h, t, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", fc=fc, ec="#333", lw=0.8))
    ax.text(x + w/2, y + h/2, t, ha="center", va="center", fontsize=8.5, wrap=True)
def arr(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=10, lw=0.8, color="#333"))
ax.text(0.1, 5.0, "a  Hypothesis", fontsize=10, weight="bold")
box(0.2, 3.6, 2.0, 0.9, "Fermentative H$_2$\nproduction\n([FeFe]-hydrogenases)", "#e8e8e8")
box(2.9, 4.2, 2.2, 0.6, "Foregut: methanogenesis\n(mcrA) dominates", "#cfe0f3")
box(2.9, 3.3, 2.2, 0.6, "Hindgut: acetogenesis (acsB)\n+ respiratory sinks", "#f6dcc2")
arr(2.25, 4.05, 2.85, 4.5); arr(2.25, 4.05, 2.85, 3.6)
ax.text(0.2, 2.95, "P3: negative covariation mcrA vs acsB (CLR, within study)", fontsize=8, style="italic")
ax.text(5.5, 5.0, "b  Design", fontsize=10, weight="bold")
box(5.6, 4.1, 4.2, 0.6, "P1 (primary): paired foregut vs hindgut, same animal\nXie 2021 PRJNA657455; Li 2025 PRJNA1068803", "#cfe0f3")
box(5.6, 3.25, 4.2, 0.6, "P2: cross-host phylogenetic GLMM\nfermentation site + phylogeny + diet + captivity + technical", "#f6dcc2")
box(5.6, 2.4, 4.2, 0.6, "Validation: marsupial hold-out\nBowerman 2025 PRJNA1167890", "#e3efd9")
ax.text(0.1, 1.95, "c  Pipeline", fontsize=10, weight="bold")
steps = ["ENA/Metalog\ninventory", "QC + subsample\n(fastp, seqkit)", "DIAMOND vs curated\nH$_2$-marker DB", "Genome-equivalent\nnormalisation", "brms / metafor\nmeta-analysis"]
for i, s in enumerate(steps):
    box(0.2 + i*1.98, 0.6, 1.7, 0.9, s, "#f2f2f2")
    if i: arr(0.2 + i*1.98 - 0.25, 1.05, 0.2 + i*1.98 - 0.02, 1.05)
ax.text(0.2, 0.15, "Conceptual schematic only; contains no data.", fontsize=7, color="#666")
plt.savefig("figures/Fig1_concept_design.png", dpi=300, bbox_inches="tight"); plt.savefig("figures/Fig1_concept_design.pdf", bbox_inches="tight")
