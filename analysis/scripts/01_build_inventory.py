#!/usr/bin/env python3
"""Build the HydroSink run inventory with PRISMA-style screening counts.

Inputs (analysis/inventory/raw/): ena_tax*.tsv (ENA Portal API search, date in search_date.txt),
host_traits.tsv (01b_host_traits.py). Seed paired-compartment studies are added explicitly.

Sets produced:
  P1  : paired foregut/hindgut compartments from the same animals (Li 2025 PRJNA1068803; Xie 2021
        PRJNA657455 rumen vs cecum/colon/rectum) -> analysis/inventory/p1_runs.tsv
  P2  : one faecal run per sample, vertebrate non-human/mouse/rat hosts, <= CAP runs per host
        species spread across studies (round-robin, seed fixed) -> analysis/inventory/p2_runs.tsv
Every exclusion step is logged in analysis/prisma_log.tsv.
"""
import pandas as pd, numpy as np, glob, re, pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
INV = ROOT / "inventory"; RAW = INV / "raw"
SEED = 20261115; CAP = 10; MIN_READS = 5_000_000; MIN_LEN = 100
EXCL_KW = re.compile(r"antibiot|in ?vitro|germ.?free|gnotobio|enrichment|culture|3-?NOP|nitrooxypropanol|bromoform|"
                     r"asparagopsis|monensin|nitrate supplement|methane inhibit|bromochloromethane|virome|vlp|"
                     r"metatranscript|rna-?seq|amplicon|16s|mock|negative control|blank|extraction control|"
                     r"fmt|transplant|spike|isolate", re.I)
FECES = re.compile(r"fec|faec|feces|faeces|stool|scat|dropping|excreta|dung|manure", re.I)
SITE_PATTERNS = [("rumen", r"rumen|ruminal"), ("reticulum", r"reticul"), ("omasum", r"omasum"),
                 ("abomasum", r"abomas"), ("forestomach", r"forestomach|foregut"), ("cecum", r"c(a)?ec(um|al|a)\b|caec"),
                 ("colon", r"colon|colonic"), ("rectum", r"rect"), ("ileum", r"ileum|ileal"),
                 ("jejunum", r"jejun"), ("duodenum", r"duoden")]
DOMESTIC = {"Bos taurus", "Bos indicus", "Bubalus bubalis", "Ovis aries", "Capra hircus", "Sus scrofa domesticus",
            "Equus caballus", "Equus asinus", "Gallus gallus", "Meleagris gallopavo", "Anas platyrhynchos",
            "Canis lupus familiaris", "Felis catus", "Oryctolagus cuniculus", "Cavia porcellus", "Bos grunniens",
            "Camelus dromedarius", "Camelus bactrianus", "Anser cygnoides", "Coturnix japonica", "Columba livia"}
log = []


def step(name, before, after, note=""):
    log.append({"stage": name, "n_before": before, "n_after": after, "excluded": before - after, "note": note})


def sample_type(r):
    txt = " ".join(str(r.get(c, "")) for c in ["host_body_site", "isolation_source", "sample_title", "scientific_name"])
    if str(r.get("scientific_name", "")) == "feces metagenome" or FECES.search(txt):
        return "feces"
    for site, pat in SITE_PATTERNS:
        if re.search(pat, txt, re.I):
            return site
    return "unspecified"


def main():
    d = pd.concat([pd.read_csv(f, sep="\t", dtype=str) for f in sorted(glob.glob(str(RAW / "ena_tax*.tsv")))])
    n0 = len(d); d = d.drop_duplicates("run_accession")
    step("identified: ENA WGS/METAGENOMIC/ILLUMINA runs, 7 gut taxa", n0, len(d), "duplicates across taxa removed")
    traits = pd.read_csv(INV / "host_traits.tsv", sep="\t", dtype=str)
    n = len(d); d = d[d.host_tax_id.notna()]
    step("host_tax_id present (free-text-only hosts: mostly domestic pig/cattle/chicken/dog, not recovered)", n, len(d))
    n = len(d); d = d.merge(traits, on="host_tax_id", how="inner")
    step("vertebrate host other than human/mouse/rat (ENA lineage)", n, len(d))
    d = d[~d.host_tax_id.isin(["9606", "10090", "10116"]) & (d.fermentation_site != "excluded_human")]
    txt = d[["sample_title", "isolation_source", "host_body_site"]].fillna("").agg(" ".join, axis=1)
    n = len(d); d = d[~txt.str.contains(EXCL_KW)]
    step("no exclusion keywords (treatment/in-vitro/non-WGS/controls)", n, len(d))
    d["read_count_n"] = pd.to_numeric(d.read_count, errors="coerce")
    d["base_count_n"] = pd.to_numeric(d.base_count, errors="coerce")
    d["read_len"] = d.base_count_n / d.read_count_n / np.where(d.library_layout.eq("PAIRED"), 2, 1)
    n = len(d); d = d[d.read_count_n >= MIN_READS]
    step(f"read_count >= {MIN_READS:,}", n, len(d))
    n = len(d); d = d[d.read_len >= MIN_LEN]
    step(f"mean read length >= {MIN_LEN} bp", n, len(d))
    d = d.sort_values("read_count_n", ascending=False)
    n = len(d); d = d.drop_duplicates("sample_accession")
    step("one run per BioSample (deepest kept)", n, len(d))
    d["sample_type"] = d.apply(sample_type, axis=1)
    d["captivity"] = np.where(d.binomial.isin(DOMESTIC) | d.sci_name.isin(DOMESTIC), "domestic",
                     np.where(txt.reindex(d.index).str.contains(r"\bzoo\b|captiv", case=False), "captive",
                     np.where(txt.reindex(d.index).str.contains(r"\bwild\b|free.?rang", case=False), "wild", "unknown")))
    d.to_csv(INV / "eligible_runs.tsv", sep="\t", index=False)

    # ---------- P2: faecal cross-host set ----------
    f = d[d.sample_type == "feces"].copy()
    step("P2: faecal samples", len(d), len(f))
    f["species"] = f.binomial
    rng = np.random.default_rng(SEED)
    picks = []
    for sp, g in f.groupby("species"):
        by_study = {s: gg.sample(frac=1, random_state=int(rng.integers(1e9))) for s, gg in g.groupby("study_accession")}
        order = sorted(by_study, key=lambda s: rng.random())
        chosen, i = [], 0
        while len(chosen) < CAP and any(len(by_study[s]) > i for s in order):
            for s in order:
                if len(by_study[s]) > i and len(chosen) < CAP:
                    chosen.append(by_study[s].iloc[i])
            i += 1
        picks += chosen
    p2 = pd.DataFrame(picks)
    step(f"P2: cap {CAP} runs/host species, round-robin across studies", len(f), len(p2))
    p2.to_csv(INV / "p2_runs.tsv", sep="\t", index=False)
    print(f"P2: {len(p2)} runs, {p2.species.nunique()} species, {p2.study_accession.nunique()} studies")
    print(p2.groupby(["vclass", "fermentation_site"]).agg(runs=("run_accession", "size"), species=("species", "nunique")).to_string())

    # ---------- P1: paired compartments ----------
    li = pd.read_csv(INV / "li2025_runs.tsv", sep="\t", dtype=str)
    li["study"] = "Li2025_PRJNA1068803"; li["host_species"] = "Capra hircus"
    li["compartment"] = np.where(li.sample_alias.str.startswith("Rum"), "rumen", "cecum")
    li["animal"] = "goat" + li.sample_alias.str.replace(r"^(Rum|Ce)", "", regex=True)
    xie = pd.read_csv(INV / "raw" / "xie2021_runs.tsv", sep="\t", dtype=str)
    m = xie.sample_alias.str.extract(r"^(?P<sp>.+)_(?P<reg>[A-Za-z]+)-(?P<idx>\d+)$")
    xie = pd.concat([xie, m], axis=1)
    xie["study"] = "Xie2021_PRJNA657455"; xie["host_species"] = xie.host
    xie["compartment"] = xie.reg.str.lower(); xie["animal"] = xie.sp + "_" + xie.idx
    keep = ["run_accession", "study", "host_species", "compartment", "animal", "sample_alias"]
    p1 = pd.concat([li[keep], xie[keep]])
    n = len(p1)
    p1["chamber"] = p1.compartment.map({"rumen": "foregut", "reticulum": "foregut", "omasum": "foregut",
                                        "cecum": "hindgut", "colon": "hindgut", "rectum": "hindgut"})
    p1 = p1[p1.chamber.notna()]
    step("P1: seed paired-compartment runs (Li 2025 + Xie 2021); fermentation chambers only "
         "(abomasum + small intestine excluded)", n, len(p1), "seed studies, added directly")
    p1.to_csv(INV / "p1_runs.tsv", sep="\t", index=False)
    print(f"P1: {len(p1)} runs, {p1.animal.nunique()} animals")
    pd.DataFrame(log).to_csv(ROOT / "prisma_log.tsv", sep="\t", index=False)
    print(pd.DataFrame(log).to_string())


if __name__ == "__main__":
    main()
