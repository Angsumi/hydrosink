# Building the H₂-sink marker reference (protocol; not yet executed)

Goal: a small, curated protein reference in which every sequence carries a **guild label**. Reads are assigned only when they pass both a global and a guild-specific threshold.

1. **Seed sequences**
   - Hydrogenases: download the classified catalytic-subunit sequences from HydDB (Søndergaard et al. 2016, *Sci Rep*; citation to verify). Keep the classes listed in `config.yaml` (`hydrogenases`). Label them as uptake (group 1 [NiFe]) or bidirectional/evolving (group 3, group 4, [FeFe] A/B).
   - Terminal reductases and marker enzymes: for each KO in `config.yaml`, pull UniRef90 members of the KO (KEGG `link/uniprot`, then UniRef90 representatives), and add proteins from GTDB r220 representative genomes annotated with KofamScan (adaptive score threshold).
   - Gut-specific diversity: run KofamScan on the SPIRE protein catalogues of D01 (Xie 2021) and D11 (Youngblut 2020) and add the hits that pass KO thresholds. This captures uncultured gut lineages missing from UniRef.
2. **Disambiguation (the critical step)**
   - *acsB* (K14138): keep only the bacterial Wood–Ljungdahl acetyl-CoA synthase β clade. Remove archaeal CdhC and sulfate-reducer homologues by phylogenetic placement (MAFFT → IQ-TREE reference tree → EPA-ng). Label sequences from known acetogen genera as high-confidence acetogenic.
   - *dsrAB*: separate the reductive (sulfate-reducing) from the oxidative (sulfur-oxidising) type by placement on the DsrAB reference tree (Müller et al. 2015; to verify).
   - *frdA* vs *sdhA*: place on the FrdA/SdhA tree; keep the anaerobic fumarate-reductase clade only.
   - *mcrA*: include anaerobic methanotroph (ANME) homologues so they can be flagged, though they are not expected in guts.
3. **Clustering**: `mmseqs easy-cluster --min-seq-id 0.9 -c 0.8`; record the guild label of each representative.
4. **Index**: `diamond makedb --in hydrosink_markers.faa -d hydrosink_markers`
5. **Benchmark (required before use)**: simulate reads with InSilicoSeq from (a) genomes carrying each marker and (b) negative-control genomes carrying near-homologues. Report per-guild precision and recall at 40/50/60% identity; choose thresholds that give precision ≥ 0.95. Write the results to `analysis/marker_benchmark.tsv`.
6. **Positive control**: on D02 (Li 2025), confirm acsB enrichment in the cecum relative to the rumen, in the direction and roughly the size reported (~12×). If it does not reproduce, stop and revise the database before running the meta-analysis.
