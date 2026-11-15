#!/usr/bin/env python3
"""Build the HydroSink DIAMOND reference + label table.

Sources (downloaded by `03_fetch_db.sh`, checksums in analysis/db/SOURCES.tsv):
  * Greening lab metabolic marker gene databases, compiled v3 (figshare 14431208,
    file 27599390, md5 986f248875b26ca378b59e58517665e6, CC-BY-4.0). All 51 genes are
    kept so that near-homologues (SdhA vs FrdA, CooS vs AcsB, ...) compete for best hits.
  * UniRef50 clusters named exactly "DNA-directed RNA polymerase subunit beta"
    (Bacteria + Archaea; UniProt release recorded in SOURCES.tsv) = single-copy normaliser.

Per-gene identity thresholds follow the Greening-lab read-level protocol as published in
Nat Commun 2025 (PMC12749492): query cover 80 %, identity 50 % except rho 40 %; nuoF,
group-4 [NiFe], mmoX, [FeFe], coxL, amoA, nxrA, rbcL 60 %; psaA 80 %; psbA, isoA, atpA,
aro, ygfK 70 %; hbsT 75 %. rpoB normaliser: 50 % (our choice; logged in decision log).

Output: analysis/db/hydrosink.faa, hydrosink.labels.tsv  (seq_id, gene, guild, subgroup, min_id)
"""
import re, sys, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "db" / "raw"
GDIR = RAW / "greening" / "Greening lab metabolic marker gene databases"

# filename keyword -> gene symbol (order matters: first match wins)
GENE_RULES = [
    ("NiFe_hydrogenase", "NiFe"), ("FeFe_hydrogenase", "FeFe"), ("Fe hydrogenase", "Fe-only"),
    ("AcsB", "acsB"), ("McrA", "mcrA"), ("Malonyl-CoA reductase", "mcr_malonyl"),
    ("DsrA", "dsrA"), ("AsrA", "asrA"), ("NrfA", "nrfA"), ("NarG", "narG"), ("NapA", "napA"),
    ("Fumarate_reductase", "frdA"), ("SdhA", "sdhA"), ("CooS", "cooS"), ("CoxL", "coxL"),
    ("FdhA", "fdhA"), ("NuoF", "nuoF"), ("RbcL", "rbcL"), ("PsaA", "psaA"), ("PsbA", "psbA"),
    ("Isoprene", "isoA"), ("Selenate", "ygfK"), ("Arsenite", "aro"), ("AtpA", "atpA"),
    ("HbsT", "hbsT"), ("HbsC", "hbsC"), ("rhodopsin", "rho"), ("MmoA", "mmoX"),
    ("PmoA", "pmoA"), ("amoA", "amoA"), ("NxrA", "nxrA"),
]
ID60 = {"nuoF", "mmoX", "FeFe", "coxL", "amoA", "nxrA", "rbcL"}
ID70 = {"psbA", "isoA", "atpA", "aro", "ygfK"}
GUILD = {"mcrA": "methanogenesis", "acsB": "acetogenesis", "dsrA": "sulfite_reduction",
         "asrA": "sulfite_reduction", "nrfA": "nitrite_ammonification", "narG": "nitrate_reduction",
         "napA": "nitrate_reduction", "frdA": "fumarate_reduction", "NiFe": "hydrogenase",
         "FeFe": "hydrogenase", "Fe-only": "hydrogenase", "rpoB": "normaliser"}


def gene_of(fname):
    for key, g in GENE_RULES:
        if key.lower() in fname.lower():
            return g
    return re.sub(r"[^A-Za-z0-9]+", "_", fname.rsplit(".", 1)[0])[:30]


def min_id(gene, subgroup):
    if gene == "rho":
        return 40
    if gene == "psaA":
        return 80
    if gene == "hbsT":
        return 75
    if gene in ID70:
        return 70
    if gene in ID60:
        return 60
    if gene == "NiFe" and subgroup.startswith("4"):
        return 60
    return 50


def read_fasta(p):
    hdr, seq = None, []
    for line in open(p, errors="replace"):
        line = line.rstrip()
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(seq)
            hdr, seq = line[1:], []
        elif line:
            seq.append(line)
    if hdr is not None:
        yield hdr, "".join(seq)


def main():
    out_faa = open(ROOT / "db" / "hydrosink.faa", "w")
    out_lab = open(ROOT / "db" / "hydrosink.labels.tsv", "w")
    out_lab.write("seq_id\tgene\tguild\tsubgroup\tmin_id\tlen_aa\tsource\n")
    n = 0
    files = sorted(p for p in GDIR.iterdir() if p.is_file())
    for f in files:
        gene = gene_of(f.name)
        for hdr, seq in read_fasta(f):
            seq = seq.replace("*", "").replace("-", "").upper()
            if len(seq) < 50:
                continue
            sub = ""
            m = re.search(r"\[(?:NiFe|FeFe)\]\s*Group\s*([0-9A-Za-z]+)", hdr)
            if m:
                sub = m.group(1)
            n += 1
            sid = f"HS{n:06d}"
            out_faa.write(f">{sid}\n{seq}\n")
            out_lab.write(f"{sid}\t{gene}\t{GUILD.get(gene, 'other')}\t{sub}\t{min_id(gene, sub)}\t{len(seq)}\t"
                          f"greening_v3:{hdr.split()[0]}\n")
    # Full-length (1000-1600 aa) bacterial UniRef90 RpoB = normaliser; RpoC (beta') kept as decoy.
    # (UniRef50 'rpoB' representatives were mostly fragments, median 135 aa -> rejected.)
    for hdr, seq in read_fasta(RAW / "uniref90_rpoB_full.faa"):
        name = re.sub(r"^\S+\s+", "", hdr).split(" n=")[0]
        if name == "DNA-directed RNA polymerase subunit beta":
            gene, guild = "rpoB", "normaliser"
        elif name == "DNA-directed RNA polymerase subunit beta'":
            gene, guild = "rpoC_decoy", "other"
        else:
            continue
        n += 1
        sid = f"HS{n:06d}"
        out_faa.write(f">{sid}\n{seq}\n")
        out_lab.write(f"{sid}\t{gene}\t{guild}\t\t50\t{len(seq)}\t{hdr.split()[0]}\n")
    out_faa.close(); out_lab.close()
    print(f"wrote {n} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
