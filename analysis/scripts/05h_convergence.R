#!/usr/bin/env Rscript
# HydroSink — POST-HOC (not pre-registered), added 2026-11-22 in response to internal review.
# Multi-chain convergence check of the P2 MCMCglmm models. Refits every P2 model (null + full for the four Gaussian
# outcomes, and the categorical mcrA-detection model) as independent chains with distinct seeds, using the functions,
# priors, iterations, burn-in and thinning of 05_meta_analysis.R unchanged. Reports Gelman-Rubin PSRF (point and upper
# 95% CI), multivariate PSRF, pooled effective sample size, and estimates pooled over chains.
# Inputs : <HS_RESULTS>/p2_matrix.tsv (802-sample primary set: analysis/results; 782-sample set: analysis/results/posthoc_no_mouse)
# Outputs: analysis/results/convergence/<label>/{diagnostics,dic_by_chain,dic,contrasts,fixed_effects,variance_partition}.tsv
# Usage  : HS_RESULTS=analysis/results HS_LABEL=n802 [HS_CHAINS=4] [HS_CORES=16] Rscript 05h_convergence.R
commandArgs <- function(...) "none"            # load the function definitions without running the main analyses
source("analysis/scripts/05_meta_analysis.R")
rm(commandArgs)
suppressPackageStartupMessages({library(coda); library(parallel)})
LABEL <- Sys.getenv("HS_LABEL", "n802"); NCH <- as.integer(Sys.getenv("HS_CHAINS", "4"))
CORES <- as.integer(Sys.getenv("HS_CORES", "16"))
OUT <- file.path(Sys.getenv("HS_CONV_OUT", "analysis/results/convergence"), LABEL); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEEDS <- 20261122L + 1000L * seq_len(NCH)

# --- analysis set: identical construction to run_p2() ---
d <- fread(file.path(R, "p2_matrix.tsv"))[qc_pass == TRUE]
d <- add_ratios(d)
d[is.na(diet_class) | diet_class == "NA", diet_class := "unknown"]
d[, fermentation_site := relevel(factor(fermentation_site), ref = "simple")]
d[, diet_class := relevel(factor(diet_class), ref = "omnivore")]
d[, captivity := relevel(factor(captivity), ref = "unknown")]
d[, log_rpoB := log10(rpoB_n)]; d[, read_len_c := as.numeric(read_len) - 150]
pt <- prep_tree(d); d <- pt$d; Ainv <- inverseA(pt$tree, nodes = "TIPS")$Ainv
cov_null <- "diet_class + captivity + log_rpoB + read_len_c"
cat(sprintf("[%s] n = %d samples, %d species, %d studies; NITT = %d; chains = %d\n", LABEL, nrow(d),
            uniqueN(d$species), uniqueN(d$study_accession), NITT, NCH))

jobs <- list()
for (y in c("alr_acet", "alr_resp", "l_acsB", "l_mcrA")) {
  jobs[[paste0(y, "|null")]] <- list(y = y, fx = paste(y, "~", cov_null), fam = "gaussian")
  jobs[[paste0(y, "|full")]] <- list(y = y, fx = paste(y, "~ fermentation_site +", cov_null), fam = "gaussian")
}
jobs[["mcrA_present|full"]] <- list(y = "mcrA_present", fx = paste("mcrA_present ~ fermentation_site +", cov_null), fam = "categorical")
if (nzchar(Sys.getenv("HS_JOBS"))) jobs <- jobs[names(jobs) %in% strsplit(Sys.getenv("HS_JOBS"), ",")[[1]]]  # e.g. "mcrA_present|full"
grid <- expand.grid(job = names(jobs), chain = seq_len(NCH), stringsAsFactors = FALSE)

t0 <- Sys.time()
fits <- mclapply(seq_len(nrow(grid)), function(i) {
  j <- jobs[[grid$job[i]]]; set.seed(SEEDS[grid$chain[i]])
  m <- fit_mcmc(as.formula(j$fx), d, Ainv, family = j$fam)
  list(Sol = m$Sol, VCV = m$VCV, DIC = m$DIC)
}, mc.cores = CORES, mc.preschedule = FALSE)
bad <- vapply(fits, function(f) inherits(f, "try-error") || is.null(f$Sol), logical(1))
if (any(bad)) stop("failed fits: ", paste(grid$job[bad], grid$chain[bad], collapse = "; "))
cat(sprintf("[%s] %d fits in %.1f min\n", LABEL, length(fits), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(list(grid = grid, fits = fits, seeds = SEEDS, nitt = NITT), file.path(OUT, "chains.rds"))

# derived quantities per draw (contrasts, variance shares), so they are diagnosed like parameters
derived <- function(f, full, cat) {
  S <- f$Sol; V <- f$VCV; tot <- rowSums(V)
  out <- cbind(S, V)
  if (!cat) out <- cbind(out, phylo_H2 = V[, "phylo"] / tot, species_share = V[, "species"] / tot, study_share = V[, "study_accession"] / tot)
  if (full) {
    g <- function(l) if (l == "simple") 0 else S[, paste0("fermentation_site", l)]
    out <- cbind(out, `foregut - hindgut` = g("foregut") - g("hindgut"), `hindgut - simple` = g("hindgut") - g("simple"),
                 `foregut - simple` = g("foregut") - g("simple"), `avian_caecal - simple` = g("avian_caecal") - g("simple"))
  }
  mcmc(out, start = 1, thin = 1)
}

diag_rows <- list(); dic_rows <- list(); con_rows <- list(); fx_rows <- list(); vp_rows <- list()
for (jn in names(jobs)) {
  j <- jobs[[jn]]; full <- grepl("\\|full$", jn); cat_ <- j$fam == "categorical"; idx <- which(grid$job == jn)
  ml <- mcmc.list(lapply(fits[idx], derived, full = full, cat = cat_))
  pooled <- as.matrix(ml)
  # drop parameters that are constant (e.g. fixed residual of the categorical model) before Gelman-Rubin
  keep <- apply(pooled, 2, var) > 0
  gd <- gelman.diag(ml[, keep], autoburnin = FALSE, multivariate = FALSE)$psrf
  fcols <- intersect(colnames(pooled)[keep], colnames(fits[[idx[1]]]$Sol))
  mpsrf <- tryCatch(gelman.diag(ml[, fcols], autoburnin = FALSE)$mpsrf,
                    error = function(e) NA_real_)
  ess <- effectiveSize(ml[, keep])
  diag_rows[[jn]] <- data.table(outcome = j$y, model = sub(".*\\|", "", jn), parameter = rownames(gd),
                                psrf = gd[, 1], psrf_upper = gd[, 2], ess = ess[rownames(gd)], mpsrf_fixed = mpsrf)
  dic_rows[[jn]] <- data.table(outcome = j$y, model = sub(".*\\|", "", jn), chain = grid$chain[idx],
                               DIC = vapply(fits[idx], `[[`, numeric(1), "DIC"))
  if (full) {
    S <- pooled
    fx_rows[[jn]] <- data.table(outcome = j$y, term = colnames(fits[[idx[1]]]$Sol),
                                post_mean = colMeans(S[, colnames(fits[[idx[1]]]$Sol), drop = FALSE]),
                                l95 = HPDinterval(mcmc(S[, colnames(fits[[idx[1]]]$Sol), drop = FALSE]))[, 1],
                                u95 = HPDinterval(mcmc(S[, colnames(fits[[idx[1]]]$Sol), drop = FALSE]))[, 2])
    if (!cat_) {
      cn <- c("foregut - hindgut", "hindgut - simple", "foregut - simple", "avian_caecal - simple")
      con_rows[[jn]] <- rbindlist(lapply(cn, function(k) { x <- S[, k]
        data.table(outcome = j$y, contrast = k, post_mean = mean(x), l95 = HPDinterval(mcmc(x))[1], u95 = HPDinterval(mcmc(x))[2],
                   pMCMC = 2 * min(mean(x > 0), mean(x < 0)),
                   chain_means = paste(sprintf("%.3f", vapply(ml, function(c) mean(c[, k]), numeric(1))), collapse = ";")) }))
      vp_rows[[jn]] <- rbindlist(lapply(c("phylo_H2", "species_share", "study_share"), function(k) { x <- mcmc(S[, k])
        data.table(outcome = j$y, component = k, post_mode = posterior.mode(x), l95 = HPDinterval(x)[1], u95 = HPDinterval(x)[2],
                   chain_modes = paste(sprintf("%.3f", vapply(ml, function(c) posterior.mode(mcmc(c[, k])), numeric(1))), collapse = ";")) }))
    }
  }
}
diag <- rbindlist(diag_rows); dicc <- rbindlist(dic_rows)
dic <- dcast(dicc, outcome + chain ~ model, value.var = "DIC")
if (!"null" %in% names(dic)) dic[, null := NA_real_]
dic[, dDIC := null - full]
dic_sum <- dic[, .(DIC_null_mean = mean(null), DIC_full_mean = mean(full), dDIC_mean = mean(dDIC),
                   dDIC_min = min(dDIC), dDIC_max = max(dDIC), dDIC_by_chain = paste(sprintf("%.2f", dDIC), collapse = ";")), by = outcome]
meta <- data.table(label = LABEL, n = nrow(d), species = uniqueN(d$species), studies = uniqueN(d$study_accession),
                   chains = NCH, nitt = NITT, burnin = round(NITT * 3 / 13), thin = max(1, round(NITT / 1300)),
                   draws_per_chain = nrow(fits[[1]]$Sol), seeds = paste(SEEDS, collapse = ";"),
                   max_psrf = max(diag$psrf), max_psrf_upper = max(diag$psrf_upper), min_ess = min(diag$ess),
                   worst_psrf_param = diag[which.max(psrf), paste(outcome, model, parameter)],
                   worst_ess_param = diag[which.min(ess), paste(outcome, model, parameter)])
fwrite(meta, file.path(OUT, "summary.tsv"), sep = "\t")
fwrite(diag, file.path(OUT, "diagnostics.tsv"), sep = "\t")
fwrite(dic, file.path(OUT, "dic_by_chain.tsv"), sep = "\t")
fwrite(dic_sum, file.path(OUT, "dic.tsv"), sep = "\t")
if (length(con_rows)) fwrite(rbindlist(con_rows), file.path(OUT, "contrasts.tsv"), sep = "\t")
fwrite(rbindlist(fx_rows), file.path(OUT, "fixed_effects.tsv"), sep = "\t")
if (length(vp_rows)) fwrite(rbindlist(vp_rows), file.path(OUT, "variance_partition.tsv"), sep = "\t")
print(meta); print(dic_sum); if (length(con_rows)) print(rbindlist(con_rows)[outcome == "alr_acet"]); print(rbindlist(vp_rows))
print(diag[order(-psrf)][1:10]); print(diag[order(ess)][1:10])
