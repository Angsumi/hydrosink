POST-HOC (not pre-registered), added 2026-11-17 after the pre-registered P2 results were seen.
Reason: the ENA screen excluded mouse by taxid 10090 only; 20 Mus musculus samples entered P2 under taxids 10092
(M. m. domesticus) and 10088 (genus Mus, resolved by curation), 10 of them from a DSS-colitis model (PRJDB41576).
Also, hosts resolved from genus to species had no EltonTraits diet (Elephas maximus: 1 sample) -> set to herbivore.
This set removes all Mus musculus and harmonises that diet label. Run: HS_RESULTS=analysis/results/posthoc_no_mouse Rscript analysis/scripts/05d_sensitivity_sets.R
