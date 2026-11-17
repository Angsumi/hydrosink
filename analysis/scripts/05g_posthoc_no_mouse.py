#!/usr/bin/env python3
"""POST-HOC set (not pre-registered; see analysis/results/posthoc_no_mouse/README.txt): P2 without Mus musculus
(entered via taxids 10092 / 10088 despite the 10090 exclusion) and with the single diet-less Elephas maximus sample set
to herbivore, as its other records. Writes analysis/results/posthoc_no_mouse/{p1_matrix,p2_matrix}.tsv."""
import pathlib, shutil
import pandas as pd
R = pathlib.Path(__file__).resolve().parents[2] / "analysis/results"; O = R / "posthoc_no_mouse"; O.mkdir(exist_ok=True)
a = pd.read_csv(R / "p2_matrix.tsv", sep="\t", dtype={"diet_class": str})
n0 = len(a); a = a[a.species != "Mus musculus"].copy(); a.loc[a.species == "Elephas maximus", "diet_class"] = "herbivore"
a.to_csv(O / "p2_matrix.tsv", sep="\t", index=False); shutil.copy(R / "p1_matrix.tsv", O / "p1_matrix.tsv")
print("removed", n0 - len(a))
