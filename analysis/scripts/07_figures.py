#!/usr/bin/env python3
"""Manuscript figures, drawn only from analysis/results/*.tsv (each figure is skipped if its inputs do not exist).

Fig. 2  P1 within-animal hindgut - foregut contrast (pooled outcomes; per-stratum ALR_acet)
Fig. 3  P2 cross-host: per-species mean ALR_acet by fermentation site; posterior contrasts
Fig. 4  P3 within-unit proportionality rho(mcrA, acsB)
Fig. S1 marker detection by fermentation site (P2)
Output: figures/fig*.pdf + .png (300 dpi)
"""
import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[2]; R = ROOT / "analysis/results"; F = ROOT / "figures"; F.mkdir(exist_ok=True)
# categorical palette (validated: dataviz validate_palette.js, light surface; categories are always text-labelled)
SITE_COL = {"foregut": "#2a78d6", "hindgut": "#eb6834", "avian_caecal": "#1baf7a", "simple": "#eda100"}
SITE_LAB = {"foregut": "Foregut", "hindgut": "Hindgut", "avian_caecal": "Avian caecal", "simple": "Simple gut"}
INK, MUTED, GRID = "#1f1f1e", "#6b6b68", "#e4e4e1"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                     "font.family": "DejaVu Sans"})
LAB = {"alr_acet": "ALR acet (acsB/mcrA) †", "alr_resp": "ALR resp (respiratory/mcrA)", "acsB": "acsB", "mcrA": "mcrA",
       "resp": "Respiratory sum", "dsrA": "dsrA", "asrA": "asrA", "nrfA": "nrfA", "narG": "narG", "napA": "napA",
       "frdA": "frdA", "NiFe1": "[NiFe] group 1", "FeFeB": "[FeFe] group B"}


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(F / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig); print("wrote", name)


def vline0(ax):
    ax.axvline(0, color=MUTED, lw=0.8, zorder=0)
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)


def fig2():
    if not (R / "p1_pooled.tsv").exists(): return
    t = pd.read_csv(R / "p1_pooled.tsv", sep="\t"); es = pd.read_csv(R / "p1_stratum_effects.tsv", sep="\t")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={"width_ratios": [1.15, 1]})
    t = t.iloc[::-1].reset_index(drop=True); y = np.arange(len(t))
    a.hlines(y, t.pi_lb, t.pi_ub, color=GRID, lw=5, label="95% prediction interval")
    a.hlines(y, t.ci_lb, t.ci_ub, color=INK, lw=1.6, label="95% CI")
    a.plot(t.estimate, y, "o", ms=5, color=INK, mec="white", mew=1)
    a.set_yticks(y, [LAB.get(o, o) for o in t.analysis]); vline0(a)
    a.set_xlabel("Hindgut − foregut (ln units)"); a.set_title("a  Pooled within-animal contrast (8 strata, 53 animals)", loc="left", fontsize=8.5)
    a.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=7)
    e = es[es.outcome == "alr_acet"].copy(); e["se"] = np.sqrt(e.vi)
    e["lab"] = e.stratum.str.replace("_PRJNA", " PRJNA").str.replace("|", ": ", regex=False)
    e = e.sort_values("yi"); yy = np.arange(len(e))
    b.hlines(yy, e.yi - 1.96 * e.se, e.yi + 1.96 * e.se, color=MUTED, lw=1.4)
    b.plot(e.yi, yy, "s", ms=5, color=SITE_COL["hindgut"], mec="white", mew=1)
    for yi_, (_, r) in zip(yy, e.iterrows()): b.text(r.yi + 1.96 * r.se + 0.1, yi_, f"n={int(r.n)}", va="center", fontsize=6.5, color=MUTED)
    p = t[t.analysis == "alr_acet"].iloc[0]
    b.axvspan(p.ci_lb, p.ci_ub, color=GRID, zorder=0); b.axvline(p.estimate, color=INK, lw=1)
    b.set_yticks(yy, e.lab, fontsize=6.5); vline0(b); b.set_xlabel("Hindgut − foregut ALR acet (ln units)")
    b.set_title("b  ALR acet per study × host species", loc="left", fontsize=8.5)
    fig.text(0.01, -0.09, "† primary outcome. Grey band in b: pooled estimate 95% CI. Source: p1_pooled.tsv, p1_stratum_effects.tsv",
             fontsize=6.5, color=MUTED)
    fig.tight_layout(); save(fig, "fig2_p1_within_animal")


def fig3():
    if not (R / "p2_species_means.tsv").exists(): return
    sp = pd.read_csv(R / "p2_species_means.tsv", sep="\t"); con = pd.read_csv(R / "p2_contrasts.tsv", sep="\t")
    order = ["foregut", "hindgut", "avian_caecal", "simple"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw={"width_ratios": [1.2, 1]})
    rng = np.random.default_rng(20261115)
    for i, s in enumerate(order):
        v = sp[sp.fermentation_site == s]
        x = i + rng.uniform(-0.18, 0.18, len(v))
        a.scatter(x, v.alr_acet, s=np.clip(v.n, 1, 10) * 4 + 6, color=SITE_COL[s], alpha=0.75, edgecolor="white", lw=0.6)
        a.hlines(v.alr_acet.median(), i - 0.3, i + 0.3, color=INK, lw=1.6)
    a.set_xticks(range(4), [f"{SITE_LAB[s]}\n({(sp.fermentation_site == s).sum()} spp.)" for s in order])
    a.axhline(0, color=MUTED, lw=0.8, zorder=0); a.grid(axis="y", color=GRID, lw=0.6); a.set_axisbelow(True)
    a.set_ylabel("Species mean ALR acet, ln(acsB/mcrA)"); a.set_title("a  Faecal ALR acet by fermentation site", loc="left", fontsize=8.5)
    c = con[con.outcome == "alr_acet"].iloc[::-1].reset_index(drop=True); y = np.arange(len(c))
    b.hlines(y, c.l95, c.u95, color=INK, lw=1.6); b.plot(c.post_mean, y, "o", ms=5, color=INK, mec="white", mew=1)
    b.set_yticks(y, [x.replace("avian_caecal", "avian caecal") for x in c.contrast]); vline0(b)
    b.set_xlabel("Posterior difference in ALR acet (95% HPD)")
    b.set_title("b  Phylogenetic mixed model contrasts", loc="left", fontsize=8.5)
    fig.text(0.01, -0.03, "a: points = host species (size ∝ samples, max 10), bar = median. b: MCMCglmm with phylogeny, species and study "
             "random effects; diet, captivity, depth, read length fixed. Source: p2_species_means.tsv, p2_contrasts.tsv",
             fontsize=6.5, color=MUTED, wrap=True)
    fig.tight_layout(); save(fig, "fig3_p2_cross_host")


def fig4():
    if not (R / "p3_units.tsv").exists(): return
    u = pd.read_csv(R / "p3_units.tsv", sep="\t"); p = pd.read_csv(R / "p3_pooled.tsv", sep="\t").iloc[0]
    u["se"] = np.sqrt(u.vz); u = u.sort_values("rho_p").reset_index(drop=True); y = np.arange(len(u))
    fig, a = plt.subplots(figsize=(4.6, 5.6))
    lo, hi = np.tanh(u.z - 1.96 * u.se), np.tanh(u.z + 1.96 * u.se)
    col = [SITE_COL["foregut"] if x.startswith("P1") else MUTED for x in u.unit]
    a.hlines(y, lo, hi, color=col, lw=1.2); a.scatter(u.rho_p, y, s=12 + u.n / 3, color=col, edgecolor="white", lw=0.6, zorder=3)
    a.set_yticks(y, u.unit.str.replace("P1:", "P1 ").str.replace("P2:", ""), fontsize=6)
    a.axvspan(p.rho_ci_lb, p.rho_ci_ub, color=GRID, zorder=0); a.axvline(p.rho_pooled, color=INK, lw=1); vline0(a)
    a.set_xlabel("Proportionality ρp(mcrA, acsB), 95% CI"); a.set_xlim(-1, 1); a.set_xticks([-1, -0.5, 0, 0.5, 1])
    a.set_title(f"Pooled ρp = {p.rho_pooled:.2f} (95% CI {p.rho_ci_lb:.2f} to {p.rho_ci_ub:.2f}; {int(p.k)} units)", loc="left", fontsize=8)
    fig.text(0.01, -0.01, "Blue: P1 study × chamber units; grey: P2 studies with ≥10 samples. A trade-off predicts ρp < 0. "
             "Source: p3_units.tsv, p3_pooled.tsv", fontsize=6.3, color=MUTED, wrap=True)
    fig.tight_layout(); save(fig, "fig4_p3_tradeoff")


def figS1():
    if not (R / "detection_by_site.tsv").exists(): return
    d = pd.read_csv(R / "detection_by_site.tsv", sep="\t", index_col=0)
    order = ["foregut", "hindgut", "avian_caecal", "simple"]; genes = [c for c in d.columns if c != "n"]
    fig, a = plt.subplots(figsize=(7.2, 2.6)); w = 0.2
    for i, s in enumerate(order):
        a.bar(np.arange(len(genes)) + (i - 1.5) * w, d.loc[s, genes], w - 0.02, color=SITE_COL[s],
              label=f"{SITE_LAB[s]} (n={int(d.loc[s, 'n'])})")
    a.set_xticks(range(len(genes)), genes); a.set_ylabel("Samples with ≥2 reads"); a.set_ylim(0, 1.05)
    a.grid(axis="y", color=GRID, lw=0.6); a.set_axisbelow(True); a.legend(frameon=False, ncol=4, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(); save(fig, "figS1_detection_by_site")


for f in (fig2, fig3, fig4, figS1):
    f()
