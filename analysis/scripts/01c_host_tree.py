#!/usr/bin/env python3
"""Host phylogeny for P2 species: Open Tree of Life synthetic tree (TNRS match -> induced subtree).
Branch lengths are added later in R (Grafen's method, ape::compute.brlen) because the OToL
synthesis tree is topology-only. Output: analysis/inventory/host_tree.nwk, host_ott_map.tsv
"""
import pandas as pd, requests, pathlib, re, json

ROOT = pathlib.Path(__file__).resolve().parents[1]
API = "https://api.opentreeoflife.org/v3"
p2 = pd.read_csv(ROOT / "inventory/p2_runs.tsv", sep="\t", dtype=str)
names = sorted(set(p2.species.dropna()))
# domestic subspecies / forms -> species name for matching
SYN = {"Mops pumilus": "Chaerephon pumilus"}
def _clean(n):
    n = SYN.get(n, n)
    n = re.sub(r"\s+sp\.$", "", n)          # genus-level hosts -> genus node
    return re.sub(r"\s+(domesticus|familiaris|taurus taurus)$", "", n)
clean = {n: _clean(n) for n in names}
rows = []
for i in range(0, len(names), 200):
    batch = [clean[n] for n in names[i:i + 200]]
    r = requests.post(f"{API}/tnrs/match_names", json={"names": batch, "do_approximate_matching": False}, timeout=120).json()
    for res in r["results"]:
        m = res["matches"][0]["taxon"] if res["matches"] else None
        rows.append(dict(query=res["name"], ott_id=m["ott_id"] if m else None, ott_name=m["unique_name"] if m else None))
mp = pd.DataFrame(rows).drop_duplicates("query")
mp = pd.DataFrame({"species": names, "query": [clean[n] for n in names]}).merge(mp, on="query", how="left")
mp.to_csv(ROOT / "inventory/host_ott_map.tsv", sep="\t", index=False)
ids = sorted({int(x) for x in mp.ott_id.dropna()})
r = requests.post(f"{API}/tree_of_life/induced_subtree", json={"ott_ids": ids, "label_format": "id"}, timeout=300)
j = r.json()
if "newick" not in j:  # some ids may be pruned/broken: drop them and retry
    bad = set()
    for k in ("unknown", "broken"):
        bad |= {int(re.sub(r"\D", "", x)) for x in (j.get(k) or {})}
    ids = [x for x in ids if x not in bad]
    j = requests.post(f"{API}/tree_of_life/induced_subtree", json={"ott_ids": ids, "label_format": "id"}, timeout=300).json()
open(ROOT / "inventory/host_tree.nwk", "w").write(j["newick"])
print(f"species {len(names)}; matched {mp.ott_id.notna().sum()}; in tree {len(ids)}")
print("unmatched:", list(mp[mp.ott_id.isna()].species))
