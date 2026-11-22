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
import matplotlib.patches

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
        fig.savefig(F / f"{name}.{ext}", dpi=300, bbox_inches="tight",
                    metadata={"CreationDate": None} if ext == "pdf" else None)   # no timestamp: byte-reproducible
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
    groups = pd.read_csv(R / "p2_species_means.tsv", sep="\t")
    groups["alr_sum"] = groups.alr_acet * groups.n
    sp = groups.groupby(["species", "fermentation_site"], as_index=False)[["n", "alr_sum"]].sum()
    sp["alr_acet"] = sp.alr_sum / sp.n
    con = pd.read_csv(R / "p2_contrasts.tsv", sep="\t")
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
    a.set_ylabel("Species mean ALR acet (adjusted log ratio)"); a.set_title("a  Faecal ALR acet by fermentation site", loc="left", fontsize=8.5)
    c = con[con.outcome == "alr_acet"].iloc[::-1].reset_index(drop=True); y = np.arange(len(c))
    b.hlines(y, c.l95, c.u95, color=INK, lw=1.6); b.plot(c.post_mean, y, "o", ms=5, color=INK, mec="white", mew=1)
    b.set_yticks(y, [x.replace("avian_caecal", "avian caecal") for x in c.contrast]); vline0(b)
    b.set_xlabel("Posterior difference in ALR acet (95% HPD)")
    b.set_title("b  Phylogenetic mixed model contrasts", loc="left", fontsize=8.5)
    fig.text(0.01, -0.03, "a: one point per host species, with area increasing with sample count (max 10). Bars show medians. b: MCMCglmm with phylogeny, species and study "
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
    fig.text(0.01, -0.01, "Blue: P1 study × chamber units; grey: P2 studies with ≥10 samples. Intervals use the bootstrap variance on the Fisher z scale. "
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


def fig1_design():
    """Study design schematic (executed design; no data values except set sizes from prisma_final.tsv)."""
    pf = pd.read_csv(R / "prisma_final.tsv", sep="\t") if (R / "prisma_final.tsv").exists() else None
    n = (lambda s: pf.loc[pf.stage.str.startswith(s), "n"].iloc[0]) if pf is not None else (lambda s: "?")
    fig, ax = plt.subplots(figsize=(7.2, 4.6)); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 64)
    def box(x, y, w, h, txt, fc, fs=6.6):
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=MUTED, lw=0.6))
        ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, color=INK)
    def arrow(x0, y0, x1, y1): ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.7))
    ax.text(0, 62, "a  Hypothesis", fontsize=8.5, weight="bold", color=INK); ax.text(50, 62, "b  Design", fontsize=8.5, weight="bold", color=INK)
    box(1, 44, 15, 10, "Fermentative H$_2$\nproduction", "#ececea")
    box(22, 50, 24, 8, "Foregut: lower predicted\nacsB/mcrA marker ratio", "#d6e6f8", 6.0); box(22, 39, 24, 8, "Hindgut: higher predicted\nacsB/mcrA marker ratio", "#fbe1d6", 6.0)
    arrow(16.5, 50, 21.5, 54); arrow(16.5, 48, 21.5, 43)
    ax.text(1, 35, "P3: negative mcrA–acsB proportionality\nwithin study units\nP4 (exploratory): more respiratory markers\nwhere mcrA is lower",
            fontsize=5.8, color=INK, style="italic", va="top")
    box(50, 50.5, 49, 8, f"P1 positive control: foregut vs hindgut, same animal\nLi 2025 + Xie 2021: {n('P1 runs passing QC')} runs, {n('P1 animals')} ruminants\nrandom-effects meta-analysis", "#d6e6f8", 5.6)
    box(50, 40, 49, 8.5, f"P2: faeces across vertebrates: {n('P2 runs in phylogenetic models')} samples,\n{str(n('P2 species')).replace(' / ', ' species, ')} studies; MCMCglmm with fermentation site,\nphylogeny, species, study, diet, captivity, technical", "#fbe1d6", 5.6)
    box(50, 30, 49, 8, "P3: within-study proportionality ρp(mcrA, acsB), pooled\nChecks: sensitivity sets, identity ±10, acsB attribution,\ndated tree, in-silico benchmark", "#dcf1e8", 5.6)
    ax.text(0, 24, "c  Pipeline", fontsize=8.5, weight="bold", color=INK)
    steps = ["ENA screening\n(PRISMA)", "5 M R1 reads\nfrom SRA (P1 head;\nP2 10 offsets)", "fastp QC\n(Q20, ≥100 bp)",
             "DIAMOND blastx vs\nGreening DB v3\n+ RpoB", "Marker RPKM /\nbacterial rpoB RPKM", "metafor,\nMCMCglmm,\nproportionality"]
    w, gap = 14.0, 2.6
    for i, t in enumerate(steps):
        x = i * (w + gap); box(x + 0.5, 6, w, 11, t, "#f4f4f2", 5.6)
        if i < len(steps) - 1: arrow(x + w + 0.7, 11.5, x + w + gap + 0.3, 11.5)
    ax.text(0.5, 0.5, "Schematic of the executed design; set sizes from analysis/results/prisma_final.tsv.", fontsize=6, color=MUTED)
    save(fig, "fig1_study_design")


def figS2_prisma():
    if not (ROOT / "analysis/prisma_log.tsv").exists() or not (R / "prisma_final.tsv").exists(): return
    pl = pd.read_csv(ROOT / "analysis/prisma_log.tsv", sep="\t"); pf = pd.read_csv(R / "prisma_final.tsv", sep="\t")
    n = lambda s: int(pf.loc[pf.stage.str.startswith(s), "n"].iloc[0])
    main = [(f"Runs identified (ENA, 7 gut-metagenome taxa)\nn = {pl.n_after[0]:,}", None)]
    for _, r in pl.iloc[1:9].iterrows():
        main.append((f"{r.stage.split(' (')[0].replace('host_tax_id present', 'structured host taxon present').replace('read_count', 'read count')}\nn = {int(r.n_after):,}", f"excluded n = {int(r.excluded):,}"))
    main += [(f"P2 runs quantified\nn = {n('P2 runs quantified')}", "no reads passing QC (76-bp reads) n = 4"),
             (f"P2 runs passing QC (rpoB reads ≥ 200)\nn = {n('P2 runs passing QC')}", f"failed QC n = {n('P2 runs quantified') - n('P2 runs passing QC')}"),
             (f"P2 analysis set (host on tree)\nn = {n('P2 runs in phylogenetic models')} samples, "
              f"{pf.loc[pf.stage.str.startswith('P2 species'), 'n'].iloc[0].replace(' / ', ' species / ')} studies", f"genus-only hosts without tree tip n = {n('P2 runs passing QC') - n('P2 runs in phylogenetic models')}")]
    fig, ax = plt.subplots(figsize=(7.2, 9.4)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, len(main) + 2.2)
    for i, (txt, exc) in enumerate(main):
        y = len(main) + 1.2 - i
        ax.text(3.2, y, txt, ha="center", va="center", fontsize=7, color=INK, bbox=dict(boxstyle="round,pad=0.35", fc="#f4f4f2", ec=MUTED, lw=0.6))
        if i < len(main) - 1: ax.annotate("", xy=(3.2, y - 0.62), xytext=(3.2, y - 0.38), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.7))
        if exc: ax.text(6.4, y + 0.5, exc, ha="left", va="center", fontsize=6.5, color=MUTED)
    p1 = pl.iloc[-1]
    ax.text(3.2, 0.25, f"P1 (paired compartments, Li 2025 + Xie 2021): {int(p1.n_before)} runs → {int(p1.n_after)} fermentation-chamber runs → "
            f"{n('P1 runs passing QC')} pass QC; {n('P1 animals')} animals with both chambers", ha="center", fontsize=6.5, color=INK, wrap=True)
    save(fig, "figS2_prisma_flow")


for f in (fig1_design, fig2, fig3, fig4, figS1, figS2_prisma):
    f()
