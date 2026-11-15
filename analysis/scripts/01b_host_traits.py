#!/usr/bin/env python3
"""Assign host traits BEFORE any marker data are inspected (blinded predictor assignment).

fermentation_site (primary predictor), assigned by lineage rules following
Stevens & Hume 1998 Physiol Rev 78:393-427 (doi:10.1152/physrev.1998.78.2.393):
  foregut       : Ruminantia, Camelidae, Hippopotamidae, Tayassuidae, Macropodidae/Potoroidae,
                  Colobinae, Bradypodidae/Choloepodidae (sloths), Opisthocomus (hoatzin)
  hindgut       : Equidae, Rhinocerotidae, Tapiridae, Elephantidae, Sirenia, Lagomorpha, Rodentia,
                  Vombatidae, Phascolarctidae, Phalangeriformes (non-macropod diprotodonts),
                  Hyracoidea, Suidae, Primates (non-colobine), herbivorous Testudines (Testudinidae),
                  Iguanidae
  avian_caecal  : Galliformes, Anseriformes, Struthioniformes/Rheiformes/Casuariiformes/Apterygiformes/
                  Tinamiformes, Gruiformes, Otidiformes, Columbiformes? (no: small caeca -> simple)
  simple        : everything else among vertebrates (Carnivora incl. giant panda, Chiroptera,
                  Eulipotyphla, Cetacea, Pholidota, Xenarthra other than sloths, most birds,
                  snakes, non-iguanid lizards, crocodiles, non-tortoise turtles, amphibians, fish)
Any taxon matching no rule at a rank finer than its class is flagged 'rule_default' for audit.

diet: EltonTraits 1.0 (Wilman et al. 2014 Ecology 95:2027; figshare 10.6084/m9.figshare.3559887)
  plant% = Fruit+Nect+Seed+PlantO ; class herbivore (>=70), carnivore (animal >=70), omnivore otherwise.
Output: analysis/inventory/host_traits.tsv
"""
import pandas as pd, pathlib, re

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "inventory" / "raw"

FOREGUT = ["Ruminantia", "Camelidae", "Hippopotamidae", "Tayassuidae", "Macropodidae", "Potoroidae",
           "Colobinae", "Bradypodidae", "Choloepodidae", "Megalonychidae", "Opisthocomus"]
HINDGUT = ["Equidae", "Rhinocerotidae", "Tapiridae", "Elephantidae", "Sirenia", "Lagomorpha", "Rodentia",
           "Vombatidae", "Phascolarctidae", "Phalangeriformes", "Diprotodontia", "Hyracoidea", "Suidae", "Primates",
           "Testudinidae", "Iguanidae"]
AVIAN_CAECAL = ["Galliformes", "Anseriformes", "Struthioniformes", "Rheiformes", "Casuariiformes",
                "Apterygiformes", "Tinamiformes", "Gruiformes", "Otidiformes"]


def site(lineage, name):
    tokens = [t.strip() for t in lineage.split(";")] + [name.split()[0] if name else ""]
    has = lambda lst: any(x in tokens for x in lst)
    if has(FOREGUT):
        return "foregut", "rule"
    if "Phalangeriformes" in tokens and ("Macropodidae" in tokens or "Potoroidae" in tokens):
        return "foregut", "rule"
    if has(HINDGUT):
        if "Hominidae" in tokens and name.startswith("Homo"):
            return "excluded_human", "rule"
        return "hindgut", "rule"
    if "Aves" in tokens:
        return ("avian_caecal", "rule") if has(AVIAN_CAECAL) else ("simple", "rule")
    if "Mammalia" in tokens:
        known_simple = ["Carnivora", "Chiroptera", "Eulipotyphla", "Cetacea", "Pholidota", "Cingulata",
                        "Vermilingua", "Dasyuromorphia", "Didelphimorphia", "Monotremata", "Peramelemorphia",
                        "Scandentia", "Afrosoricida", "Macroscelidea", "Tubulidentata"]
        return ("simple", "rule") if has(known_simple) else ("simple", "rule_default")
    return "simple", "rule"


def main():
    h = pd.read_csv(RAW / "host_taxa_lineage.tsv", sep="\t", dtype=str)
    h = h[h.lineage.fillna("").str.contains("Vertebrata")].copy()
    h["vclass"] = h.lineage.apply(lambda l: next((c for c in ["Mammalia", "Aves", "Lepidosauria", "Testudines",
                                   "Crocodylia", "Amphibia", "Actinopterygii", "Chondrichthyes"] if c in l), "other"))
    h[["fermentation_site", "site_basis"]] = h.apply(lambda r: pd.Series(site(r.lineage, r.sci_name or "")), axis=1)
    # species binomial for trait matching
    h["binomial"] = h.sci_name.fillna("").apply(lambda s: " ".join(s.split()[:2]) if len(s.split()) >= 2 else s)
    mam = pd.read_csv(RAW / "MamFuncDat.txt", sep="\t", encoding="latin-1", dtype=str)
    bird = pd.read_csv(RAW / "BirdFuncDat.txt", sep="\t", encoding="latin-1", dtype=str)
    def diet_tab(df, plant_cols, animal_cols):
        d = df[["Scientific"] + plant_cols + animal_cols].dropna(subset=["Scientific"]).copy()
        for c in plant_cols + animal_cols:
            d[c] = pd.to_numeric(d[c], errors="coerce")
        d["plant_pct"] = d[plant_cols].sum(axis=1); d["animal_pct"] = d[animal_cols].sum(axis=1)
        return d[["Scientific", "plant_pct", "animal_pct"]]
    md = diet_tab(mam, ["Diet-Fruit", "Diet-Nect", "Diet-Seed", "Diet-PlantO"],
                  ["Diet-Inv", "Diet-Vend", "Diet-Vect", "Diet-Vfish", "Diet-Vunk", "Diet-Scav"])
    bd = diet_tab(bird, ["Diet-Fruit", "Diet-Nect", "Diet-Seed", "Diet-PlantO"],
                  ["Diet-Inv", "Diet-Vend", "Diet-Vect", "Diet-Vfish", "Diet-Vunk", "Diet-Scav"])
    diet = pd.concat([md, bd]).drop_duplicates("Scientific")
    # domestic forms: map to wild ancestor names used by EltonTraits
    alias = {"Bos taurus": "Bos taurus", "Sus scrofa domesticus": "Sus scrofa", "Canis lupus familiaris": "Canis lupus",
             "Felis catus": "Felis silvestris", "Capra hircus": "Capra hircus", "Ovis aries": "Ovis aries",
             "Equus caballus": "Equus caballus", "Gallus gallus": "Gallus gallus", "Bos indicus": "Bos taurus"}
    h["diet_key"] = h.binomial.map(lambda s: alias.get(s, s))
    h = h.merge(diet, left_on="diet_key", right_on="Scientific", how="left").drop(columns=["Scientific"])
    h["diet_class"] = h.apply(lambda r: "NA" if pd.isna(r.plant_pct) else
                              ("herbivore" if r.plant_pct >= 70 else "carnivore" if r.animal_pct >= 70 else "omnivore"), axis=1)
    out = ROOT / "inventory" / "host_traits.tsv"
    h.drop(columns=["lineage"]).to_csv(out, sep="\t", index=False)
    print(h.groupby(["vclass", "fermentation_site"]).size().to_string())
    print("diet matched:", (h.diet_class != "NA").sum(), "/", len(h))
    print("rule_default (audit):", (h.site_basis == "rule_default").sum())


if __name__ == "__main__":
    main()
