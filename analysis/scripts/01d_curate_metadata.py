#!/usr/bin/env python3
"""[prereg-v1.1] Outcome-blind metadata curation for the P2 inventory (written 2026-11-16, before any P2 marker data
were assembled or inspected). Sources: ENA study title/description (raw/studies/p2_study_meta.tsv) and BioSample
attribute values (raw/studies/p2_sample_attrs.xml). Rules are keyword-based and applied per sample (sample text) then
per study (study text); manual per-study decisions are in config/captivity_manual.tsv with an evidence quote.
Output: inventory/p2_curation.tsv (run_accession, captivity_v11, captivity_basis, species_v11, species_basis)."""
import pandas as pd, re, pathlib, xml.etree.ElementTree as ET
ROOT = pathlib.Path(__file__).resolve().parents[2]; I = ROOT / "analysis/inventory"
d = pd.read_csv(I / "p2_runs.tsv", sep="\t", dtype=str)
st = pd.read_csv(I / "raw/studies/p2_study_meta.tsv", sep="\t", dtype=str).fillna("")
st["text"] = st.study_title + " " + st.study_description
stext = dict(zip(st.study_accession, st.text))
samp = {}
try:
    raw = (I / "raw/studies/p2_sample_attrs.xml").read_text()
    for blk in re.findall(r"<SAMPLE .*?</SAMPLE>", raw, flags=re.S):
        acc = re.search(r'accession="([^"]+)"', blk).group(1)
        samp[acc] = " ".join(re.findall(r"<VALUE>([^<]*)</VALUE>", blk)) + " " + " ".join(re.findall(r"<TITLE>([^<]*)</TITLE>", blk))
except FileNotFoundError:
    pass
RULES = [  # order matters: first match wins
    ("wild", r"\bwild\b(?![- ]type)|free[- ]?(ranging|living|roaming)|wild[- ]caught|natural habitat|field[- ]collected|nature reserve|national park"),
    ("captive", r"\bzoo\b|zoological|captiv|sanctuary|breeding (center|centre|base)|aquarium|rescue cent"),
    ("domestic", r"\bfarm(ed)?\b|livestock|dairy|feedlot|\bpig farm|poultry|broiler|laying hen|\blayers?\b|\bpets?\b|household|companion animal"),
    ("laboratory", r"laborator|\blab\b|germ[- ]free|specific[- ]pathogen|\bSPF\b|C57BL|BALB|Sprague|Wistar|inbred"),
]
def classify(t):
    for lab, rx in RULES:
        m = re.search(rx, t, flags=re.I)
        if m: return lab, m.group(0)
    return None, None
man = {}
mf = ROOT / "analysis/config/captivity_manual.tsv"
if mf.exists():
    mm = pd.read_csv(mf, sep="\t", dtype=str); man = {(a if pd.isna(sp) or sp == "" else a + "|" + sp): (c, e) for a, sp, c, e in mm[["study_accession", "species", "captivity", "evidence"]].itertuples(index=False)}
cap, basis = [], []
for _, r in d.iterrows():
    key = r.study_accession + "|" + str(r.species)
    if key in man or r.study_accession in man:
        c, e = man.get(key, man.get(r.study_accession)); cap.append(c); basis.append(f"manual:{e}"); continue
    c, e = classify(samp.get(r.sample_accession, ""))
    if c: cap.append(c); basis.append(f"sample:'{e}'"); continue
    c, e = classify(stext.get(r.study_accession, ""))
    # study-level 'wild' never overrides a domestic species (e.g. horses in an equid survey)
    if c == "wild" and r.captivity == "domestic": c = None
    if c: cap.append(c); basis.append(f"study:'{e}'"); continue
    if r.captivity == "domestic": cap.append("domestic"); basis.append("v1:domestic species list"); continue
    cap.append("unknown"); basis.append("none")
d["captivity_v11"], d["captivity_basis"] = cap, basis
# species labels: genus-only hosts resolved per study from text, else kept at genus with a flag
SP = pd.read_csv(ROOT / "analysis/config/species_manual.tsv", sep="\t", dtype=str) if (ROOT / "analysis/config/species_manual.tsv").exists() else pd.DataFrame(columns=["study_accession", "species", "species_v11", "evidence"])
spm = {(a, b): (c, e) for a, b, c, e in SP[["study_accession", "species", "species_v11", "evidence"]].itertuples(index=False)}
d["species_v11"] = [spm.get((a, b), (b, None))[0] for a, b in zip(d.study_accession, d.species)]
d["species_basis"] = [("manual:" + spm[(a, b)][1]) if (a, b) in spm else ("binomial" if re.match(r"^[A-Z][a-z]+ [a-z]+$", str(b)) and not str(b).endswith(" sp.") else "genus_only_unresolved")
                      for a, b in zip(d.study_accession, d.species)]
out = d[["run_accession", "study_accession", "species", "captivity", "captivity_v11", "captivity_basis", "species_v11", "species_basis"]]
out.to_csv(I / "p2_curation.tsv", sep="\t", index=False)
print("captivity v1 :", d.captivity.value_counts().to_dict()); print("captivity v11:", d.captivity_v11.value_counts().to_dict())
print("species unresolved:", (d.species_basis == "genus_only_unresolved").sum())
