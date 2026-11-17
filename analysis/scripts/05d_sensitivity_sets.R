#!/usr/bin/env Rscript
# HydroSink — runs the pre-registered models on an alternative results directory (sensitivity data sets, e.g.
# identity threshold -10 / +10, analysis/results/sens_id_m10). Uses the functions of 05_meta_analysis.R unchanged:
#   P1 (full), P3 (full), and for P2 the primary outcome only (alr_acet: null vs full model, DIC, contrasts).
# Usage: HS_RESULTS=analysis/results/sens_id_m10 Rscript 05d_sensitivity_sets.R
commandArgs <- function(...) "none"            # load the function definitions without running the main analyses
source("analysis/scripts/05_meta_analysis.R")
rm(commandArgs)
run_p1(); run_p3()
d <- fread(file.path(R, "p2_matrix.tsv"))[qc_pass == TRUE]
d <- add_ratios(d)
d[is.na(diet_class) | diet_class == "NA", diet_class := "unknown"]
d[, fermentation_site := relevel(factor(fermentation_site), ref = "simple")]
d[, diet_class := relevel(factor(diet_class), ref = "omnivore")]
d[, captivity := relevel(factor(captivity), ref = "unknown")]
d[, log_rpoB := log10(rpoB_n)]; d[, read_len_c := as.numeric(read_len) - 150]
pt <- prep_tree(d); d <- pt$d; Ainv <- inverseA(pt$tree, nodes = "TIPS")$Ainv
cov_null <- "diet_class + captivity + log_rpoB + read_len_c"
m0 <- fit_mcmc(as.formula(paste("alr_acet ~", cov_null)), d, Ainv)
m1 <- fit_mcmc(as.formula(paste("alr_acet ~ fermentation_site +", cov_null)), d, Ainv)
out <- cbind(n = nrow(d), DIC_null = m0$DIC, DIC_full = m1$DIC, dDIC = m0$DIC - m1$DIC,
             rbind(contrast(m1, "foregut", "hindgut", "foregut - hindgut"), contrast(m1, "hindgut", "simple", "hindgut - simple"),
                   contrast(m1, "foregut", "simple", "foregut - simple"), contrast(m1, "avian_caecal", "simple", "avian_caecal - simple")))
fwrite(out, file.path(R, "p2_primary.tsv"), sep = "\t"); print(out)
