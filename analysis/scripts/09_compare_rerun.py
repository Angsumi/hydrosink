#!/usr/bin/env python3
"""Reproducibility check after run_all.sh: compare every tracked analysis/results/**/*.tsv with its committed version
(git HEAD). Reports identical files, and for differing files the maximum absolute and relative numeric difference.
Output: analysis/results/rerun_comparison.tsv"""
import io, pathlib, subprocess
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
files = subprocess.run(["git", "ls-files", "analysis/results"], cwd=ROOT, capture_output=True, text=True).stdout.split()
files = [f for f in files if f.endswith(".tsv") and "/per_sample/" not in f and "/bench/" not in f]
rows = []
for f in sorted(files):
    new = (ROOT / f).read_bytes() if (ROOT / f).exists() else None
    old = subprocess.run(["git", "show", f"HEAD:{f}"], cwd=ROOT, capture_output=True).stdout
    if new is None: rows.append(dict(file=f, status="MISSING")); continue
    if new == old: rows.append(dict(file=f, status="identical")); continue
    try:
        a = pd.read_csv(io.BytesIO(old), sep="\t"); b = pd.read_csv(io.BytesIO(new), sep="\t")
        if a.shape != b.shape or list(a.columns) != list(b.columns):
            rows.append(dict(file=f, status="DIFFERENT_SHAPE", detail=f"{a.shape} vs {b.shape}")); continue
        na, nb = a.select_dtypes("number"), b[a.select_dtypes("number").columns]
        same_text = a.select_dtypes(exclude="number").fillna("").equals(b.select_dtypes(exclude="number").fillna(""))
        d = (na - nb).abs().to_numpy(); r = d / np.maximum(na.abs().to_numpy(), 1e-12)
        rows.append(dict(file=f, status="numeric_diff" if same_text else "TEXT_DIFF",
                         max_abs_diff=np.nanmax(d) if d.size else 0, max_rel_diff=np.nanmax(r) if r.size else 0))
    except Exception as e:
        rows.append(dict(file=f, status="UNPARSEABLE_DIFF", detail=str(e)[:80]))
out = pd.DataFrame(rows); out.to_csv(ROOT / "analysis/results/rerun_comparison.tsv", sep="\t", index=False)
print(out.status.value_counts().to_string()); print(out[out.status != "identical"].to_string(index=False))
