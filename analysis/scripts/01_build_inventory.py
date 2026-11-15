#!/usr/bin/env python3
"""Build a machine-readable run-level inventory from seed accessions (datasets.csv) and,
optionally, a Metalog animal metadata export. Applies the pre-registered inclusion rules and
writes PRISMA-style counts.

Requires network access to the ENA Portal API. Not executed in the authoring environment.

Usage:
  python 01_build_inventory.py --seeds ../../datasets.csv --config ../config/config.yaml \
      [--metalog metalog_animal.tsv] --out ../inventory.tsv --prisma ../prisma_log.tsv
"""
import argparse, io, re, sys, time
import pandas as pd
import requests
import yaml

ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = ("run_accession,sample_accession,secondary_sample_accession,study_accession,"
          "scientific_name,host_tax_id,host_scientific_name,instrument_platform,instrument_model,"
          "library_strategy,library_source,library_layout,read_count,base_count,nominal_length,"
          "sample_title,sample_description,collection_date,country,fastq_ftp,fastq_md5")
ACC_RE = re.compile(r"^(PRJ[EDN][A-Z]?\d+|[ESD]RP\d+)$")


def fetch_runs(acc: str, retries: int = 3) -> pd.DataFrame:
    params = {"accession": acc, "result": "read_run", "fields": FIELDS, "format": "tsv"}
    for i in range(retries):
        r = requests.get(ENA, params=params, timeout=60)
        if r.status_code == 200 and r.text.strip():
            return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str)
        time.sleep(2 ** i)
    print(f"[warn] no runs returned for {acc} (HTTP {r.status_code})", file=sys.stderr)
    return pd.DataFrame()


def apply_filters(df: pd.DataFrame, cfg: dict, log: list) -> pd.DataFrame:
    inc = cfg["inclusion"]
    def step(name, mask):
        nonlocal df
        before = len(df)
        df = df[mask(df)].copy()
        log.append({"step": name, "n_before": before, "n_after": len(df), "excluded": before - len(df)})
    step("library_strategy==WGS", lambda d: d["library_strategy"].str.upper().eq(inc["library_strategy"]))
    step("library_source==METAGENOMIC", lambda d: d["library_source"].str.upper().eq(inc["library_source"]))
    step("platform Illumina", lambda d: d["instrument_platform"].str.upper().isin(inc["platforms"]))
    step("not human/mouse/rat host", lambda d: ~d["host_tax_id"].fillna("").isin([str(x) for x in inc["exclude_host_taxids"]]))
    kw = re.compile("|".join(re.escape(k) for k in inc["exclude_keywords"]), re.I)
    step("no exclusion keywords", lambda d: ~(d["sample_title"].fillna("") + " " + d["sample_description"].fillna("")).str.contains(kw))
    rc = pd.to_numeric(df["read_count"], errors="coerce")
    df["read_count_num"] = rc
    step("read_count >= min (pre-QC proxy)", lambda d: d["read_count_num"].fillna(0) >= inc["min_read_pairs_after_qc"])
    step("deduplicate sample_accession (keep deepest run)",
         lambda d: ~d.sort_values("read_count_num", ascending=False).duplicated("sample_accession", keep="first").reindex(d.index))
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--metalog")
    ap.add_argument("--out", required=True)
    ap.add_argument("--prisma", required=True)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    seeds = pd.read_csv(a.seeds, dtype=str)
    accs = sorted({x.strip() for cell in seeds["accession"].dropna() for x in re.split(r"[+;,]", cell) if ACC_RE.match(x.strip())})
    if a.metalog:
        ml = pd.read_csv(a.metalog, sep="\t", dtype=str)
        col = next(c for c in ml.columns if c.lower() in {"study_accession", "study", "bioproject"})
        accs = sorted(set(accs) | set(ml[col].dropna().unique()))
    frames = []
    for acc in accs:
        d = fetch_runs(acc)
        if len(d):
            d["seed_accession"] = acc
            frames.append(d)
    runs = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=FIELDS.split(","))
    log = [{"step": "records identified (runs)", "n_before": None, "n_after": len(runs), "excluded": 0}]
    inv = apply_filters(runs, cfg, log)
    if a.metalog:
        key = "sample_accession" if "sample_accession" in ml.columns else ml.columns[0]
        inv = inv.merge(ml, how="left", left_on="sample_accession", right_on=key, suffixes=("", "_metalog"))
    inv.to_csv(a.out, sep="\t", index=False)
    pd.DataFrame(log).to_csv(a.prisma, sep="\t", index=False)
    print(f"{len(accs)} accessions -> {len(runs)} runs -> {len(inv)} included runs")


if __name__ == "__main__":
    main()
