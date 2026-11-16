#!/usr/bin/env Rscript
# HydroSink — pre-specified statistical analyses (written before P1/P2 marker data were inspected).
# Inputs : analysis/results/p1_matrix.tsv, p2_matrix.tsv (05a_assemble.py); analysis/inventory/host_tree.nwk,
#          host_ott_map.tsv, host_traits.tsv (lineage for vclass/order)
# Outputs: analysis/results/*.tsv (every number quoted in the manuscript comes from these files)
# Usage  : Rscript 05_meta_analysis.R [p1|p2|p3|all]  (default all)
suppressPackageStartupMessages({library(data.table); library(metafor); library(lme4); library(MCMCglmm); library(ape)})
set.seed(20261115)
args <- commandArgs(TRUE); what <- if (length(args)) args[1] else "all"
R <- Sys.getenv("HS_RESULTS", "analysis/results")      # override only for code testing on synthetic data
NITT <- as.integer(Sys.getenv("HS_NITT", "130000"))    # override only for code testing
hm <- function(x) 0.5 * min(x[x > 0], na.rm = TRUE)          # pseudocount: half min non-zero
MK <- c("mcrA", "acsB", "dsrA", "asrA", "nrfA", "narG", "napA", "frdA", "NiFe", "FeFe", "rpoB")
add_ratios <- function(d) {
  cA <- hm(d$acsB_cpg); cM <- hm(d$mcrA_cpg); cR <- hm(d$resp_cpg)
  d[, alr_acet := log((acsB_cpg + cA) / (mcrA_cpg + cM))]
  d[, alr_resp := log((resp_cpg + cR) / (mcrA_cpg + cM))]
  d[, l_acsB := log(acsB_cpg + cA)][, l_mcrA := log(mcrA_cpg + cM)][, l_resp := log(resp_cpg + cR)]
  for (g in c("dsrA", "asrA", "nrfA", "narG", "napA", "frdA", "NiFe1", "FeFeB")) {
    cg <- hm(d[[paste0(g, "_cpg")]]); if (!is.finite(cg)) cg <- 1e-6
    d[, (paste0("l_", g)) := log(get(paste0(g, "_cpg")) + cg)]
  }
  d[, mcrA_present := as.integer(mcrA_n >= 2)]
  d[]
}
pool_tab <- function(r, label) data.table(analysis = label, k = r$k, estimate = r$b[1], se = r$se, ci_lb = r$ci.lb,
  ci_ub = r$ci.ub, pi_lb = predict(r)$pi.lb, pi_ub = predict(r)$pi.ub, p = r$pval, tau2 = r$tau2, I2 = r$I2, QEp = r$QEp)

# ============================ P1: within-animal foregut vs hindgut ============================
run_p1 <- function() {
  d <- fread(file.path(R, "p1_matrix.tsv"))[qc_pass == TRUE]
  d <- add_ratios(d)
  d[, stratum := paste(study, host_species, sep = "|")]
  d[, chamber := factor(chamber, levels = c("foregut", "hindgut"))]
  outcomes <- c(alr_acet = "alr_acet", alr_resp = "alr_resp", acsB = "l_acsB", mcrA = "l_mcrA", resp = "l_resp",
                dsrA = "l_dsrA", asrA = "l_asrA", nrfA = "l_nrfA", narG = "l_narG", napA = "l_napA", frdA = "l_frdA",
                NiFe1 = "l_NiFe1", FeFeB = "l_FeFeB")
  # stage 1: per animal mean(hindgut) - mean(foregut); per stratum mean and variance
  per_animal <- rbindlist(lapply(names(outcomes), function(o) {
    v <- outcomes[[o]]
    a <- d[, .(fg = mean(get(v)[chamber == "foregut"]), hg = mean(get(v)[chamber == "hindgut"])), by = .(stratum, animal)]
    a[is.finite(fg) & is.finite(hg)][, `:=`(outcome = o, diff = hg - fg)]
  }))
  fwrite(per_animal, file.path(R, "p1_per_animal.tsv"), sep = "\t")
  es <- per_animal[, .(yi = mean(diff), vi = var(diff) / .N, n = .N, n_pos = sum(diff > 0)), by = .(outcome, stratum)]
  fwrite(es, file.path(R, "p1_stratum_effects.tsv"), sep = "\t")
  # stage 2: REML + Hartung-Knapp; prediction interval; BH across marker outcomes (not the primary)
  tab <- rbindlist(lapply(unique(es$outcome), function(o) {
    e <- es[outcome == o]; r <- rma(yi, vi, data = e, method = "REML", test = "knha")
    cbind(pool_tab(r, o), n_animals = sum(e$n), animals_positive = sum(e$n_pos))
  }))
  tab[, role := fifelse(analysis == "alr_acet", "primary", "secondary")]
  tab[role == "secondary", q_BH := p.adjust(p, "BH")]
  fwrite(tab, file.path(R, "p1_pooled.tsv"), sep = "\t")
  # leave-one-stratum-out for the primary outcome
  e <- es[outcome == "alr_acet"]; r <- rma(yi, vi, data = e, method = "REML", test = "knha")
  l1 <- as.data.table(leave1out(r)); l1[, left_out := e$stratum]
  fwrite(l1, file.path(R, "p1_leave1out.tsv"), sep = "\t")
  # one-stage mixed model cross-check
  m <- lmer(alr_acet ~ chamber + (1 | animal) + (1 + chamber | stratum), data = d,
            control = lmerControl(optimizer = "bobyqa"))
  cf <- summary(m)$coefficients
  ci <- confint(m, parm = "chamberhindgut", method = "Wald")
  fwrite(data.table(term = rownames(cf), cf, ci_lb = c(NA, ci[1]), ci_ub = c(NA, ci[2])),
         file.path(R, "p1_mixed_model.tsv"), sep = "\t")
  # sensitivity: rumen-vs-cecum only; rumen-vs-rectum (faecal proxy)
  sens <- rbindlist(lapply(list(rumen_cecum = c("rumen", "cecum"), rumen_rectum = c("rumen", "rectum")), function(cc) {
    dd <- d[compartment %in% cc]
    a <- dd[, .(diff = mean(alr_acet[compartment == cc[2]]) - mean(alr_acet[compartment == cc[1]])), by = .(stratum, animal)][is.finite(diff)]
    e <- a[, .(yi = mean(diff), vi = var(diff) / .N), by = stratum]
    pool_tab(rma(yi, vi, data = e, method = "REML", test = "knha"), paste(cc, collapse = "_vs_"))
  }))
  fwrite(sens, file.path(R, "p1_sensitivity.tsv"), sep = "\t")
  # [prereg-v1.1] P1 has only TWO source studies (Li 2025, Xie 2021). Strata (study x host species) are not
  # independent replication. Report (i) each study separately and (ii) a three-level model with strata nested
  # in study; the stratum-level pooling above is descriptive only.
  e <- es[outcome == "alr_acet"]; e[, study := sub("\\|.*", "", stratum)]
  per_study <- rbindlist(lapply(unique(e$study), function(s) {
    x <- e[study == s]
    if (nrow(x) == 1) data.table(analysis = paste0("study:", s), k = 1, estimate = x$yi, se = sqrt(x$vi),
      ci_lb = x$yi - 1.96 * sqrt(x$vi), ci_ub = x$yi + 1.96 * sqrt(x$vi), p = 2 * pnorm(-abs(x$yi / sqrt(x$vi))))
    else { r <- rma(yi, vi, data = x, method = "REML", test = "knha")
      data.table(analysis = paste0("study:", s), k = r$k, estimate = r$b[1], se = r$se, ci_lb = r$ci.lb, ci_ub = r$ci.ub, p = r$pval) }
  }), fill = TRUE)
  m3 <- rma.mv(yi, vi, random = ~ 1 | study / stratum, data = e, method = "REML", test = "t")
  per_study <- rbind(per_study, data.table(analysis = "three-level (stratum in study)", k = m3$k, estimate = m3$b[1],
    se = m3$se, ci_lb = m3$ci.lb, ci_ub = m3$ci.ub, p = m3$pval, sigma2_study = m3$sigma2[1], sigma2_stratum = m3$sigma2[2]), fill = TRUE)
  fwrite(per_study, file.path(R, "p1_per_study.tsv"), sep = "\t")
  print(tab[, .(analysis, k, estimate, ci_lb, ci_ub, pi_lb, pi_ub, p, I2, animals_positive, n_animals)])
  print(per_study)
}

# ============================ P2: across hosts (faeces), phylogenetic GLMM ============================
prep_tree <- function(d) {
  tr <- read.tree("analysis/inventory/host_tree.nwk")
  mp <- fread("analysis/inventory/host_ott_map.tsv")
  mp[, tip := paste0("ott", ott_id)]
  # genus-level / internal-node matches: tips labelled 'mrcaottX...ottY' are matched by suffix search
  miss <- setdiff(mp$tip, tr$tip.label)
  for (t in miss) { hit <- grep(paste0(t, "($|[^0-9])"), tr$tip.label, value = TRUE); if (length(hit)) mp[tip == t, tip := hit[1]] }
  tr$node.label <- NULL
  tr <- collapse.singles(tr); tr <- multi2di(tr, random = FALSE)
  tr <- compute.brlen(tr, method = "Grafen")
  d <- merge(d, mp[, .(species, phylo = tip)], by = "species")
  d <- d[phylo %in% tr$tip.label]
  tr <- keep.tip(tr, unique(d$phylo))
  list(d = d, tree = tr)
}
pri3 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
  G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
fit_mcmc <- function(fx, d, Ainv, family = "gaussian", nitt = NITT, burnin = round(NITT * 3 / 13), thin = max(1, round(NITT / 1300))) {
  pr <- pri3; if (family == "categorical") pr$R <- list(V = 1, fix = 1)
  MCMCglmm(fx, random = ~ phylo + species + study_accession, ginverse = list(phylo = Ainv), family = family,
           data = as.data.frame(d), prior = pr, nitt = nitt, burnin = burnin, thin = thin, verbose = FALSE, pr = FALSE)
}
summ_fixed <- function(m, label) {
  s <- summary(m)$solutions
  data.table(model = label, term = rownames(s), post_mean = s[, 1], l95 = s[, 2], u95 = s[, 3], pMCMC = s[, 5])
}
phylo_signal <- function(m) {
  V <- m$VCV; tot <- rowSums(V)
  h <- V[, "phylo"] / tot
  data.table(component = c("phylo_H2", "species", "study"),
             post_mode = c(posterior.mode(mcmc(h)), posterior.mode(mcmc(V[, "species"] / tot)), posterior.mode(mcmc(V[, "study_accession"] / tot))),
             l95 = c(HPDinterval(mcmc(h))[1], HPDinterval(mcmc(V[, "species"] / tot))[1], HPDinterval(mcmc(V[, "study_accession"] / tot))[1]),
             u95 = c(HPDinterval(mcmc(h))[2], HPDinterval(mcmc(V[, "species"] / tot))[2], HPDinterval(mcmc(V[, "study_accession"] / tot))[2]))
}
contrast <- function(m, a, b, label) {   # posterior of level a minus level b (reference = simple)
  S <- m$Sol; get <- function(l) if (l == "simple") 0 else S[, paste0("fermentation_site", l)]
  x <- get(a) - get(b)
  data.table(contrast = label, post_mean = mean(x), l95 = HPDinterval(mcmc(x))[1], u95 = HPDinterval(mcmc(x))[2],
             pMCMC = 2 * min(mean(x > 0), mean(x < 0)))
}
run_p2 <- function() {
  d <- fread(file.path(R, "p2_matrix.tsv"))[qc_pass == TRUE]
  d <- add_ratios(d)
  d[is.na(diet_class) | diet_class == "NA", diet_class := "unknown"]
  d[, fermentation_site := relevel(factor(fermentation_site), ref = "simple")]
  d[, diet_class := relevel(factor(diet_class), ref = "omnivore")]
  d[, captivity := relevel(factor(captivity), ref = "unknown")]
  d[, log_rpoB := log10(rpoB_n)]; d[, read_len_c := as.numeric(read_len) - 150]
  pt <- prep_tree(d); d <- pt$d; Ainv <- inverseA(pt$tree, nodes = "TIPS")$Ainv
  fwrite(d, file.path(R, "p2_analysis_set.tsv"), sep = "\t")
  cov_null <- "diet_class + captivity + log_rpoB + read_len_c"
  res <- list(); sig <- list(); dic <- list(); con <- list()
  for (y in c("alr_acet", "alr_resp", "l_acsB", "l_mcrA")) {
    m0 <- fit_mcmc(as.formula(paste(y, "~", cov_null)), d, Ainv)
    m1 <- fit_mcmc(as.formula(paste(y, "~ fermentation_site +", cov_null)), d, Ainv)
    res[[y]] <- summ_fixed(m1, y)
    sig[[y]] <- cbind(outcome = y, phylo_signal(m1))
    dic[[y]] <- data.table(outcome = y, DIC_null = m0$DIC, DIC_full = m1$DIC, dDIC = m0$DIC - m1$DIC)
    con[[y]] <- cbind(outcome = y, rbind(contrast(m1, "foregut", "hindgut", "foregut - hindgut"),
                                          contrast(m1, "hindgut", "simple", "hindgut - simple"),
                                          contrast(m1, "foregut", "simple", "foregut - simple"),
                                          contrast(m1, "avian_caecal", "simple", "avian_caecal - simple")))
    if (y == "alr_acet") saveRDS(m1, file.path(R, "p2_alr_acet_full.rds"))
  }
  # hurdle component: mcrA detected (>=2 reads)
  mh <- fit_mcmc(as.formula(paste("mcrA_present ~ fermentation_site +", cov_null)), d, Ainv, family = "categorical")
  res[["mcrA_present"]] <- summ_fixed(mh, "mcrA_present")
  fwrite(rbindlist(res), file.path(R, "p2_fixed_effects.tsv"), sep = "\t")
  fwrite(rbindlist(sig), file.path(R, "p2_variance_partition.tsv"), sep = "\t")
  fwrite(rbindlist(dic), file.path(R, "p2_model_comparison.tsv"), sep = "\t")
  fwrite(rbindlist(con), file.path(R, "p2_contrasts.tsv"), sep = "\t")
  # descriptive: per-species means (for figures)
  sp <- d[, .(n = .N, alr_acet = mean(alr_acet), alr_resp = mean(alr_resp), mcrA_det = mean(mcrA_present)),
          by = .(species, phylo, vclass, fermentation_site, diet_class)]
  fwrite(sp, file.path(R, "p2_species_means.tsv"), sep = "\t")
  # sensitivity analyses on the primary outcome
  sets <- list(mammals_only = d[vclass == "Mammalia"], non_domestic = d[captivity != "domestic"],
               cap3_per_species = d[, .SD[sample(.N, min(.N, 3))], by = species],
               rpoB_ge_1000 = d[rpoB_n >= 1000], herbivores_only = d[diet_class == "herbivore"])
  sens <- rbindlist(lapply(names(sets), function(nm) {
    s <- copy(sets[[nm]]); s[, fermentation_site := droplevels(fermentation_site)]
    s[, diet_class := droplevels(diet_class)]; s[, captivity := droplevels(captivity)]
    tr <- keep.tip(pt$tree, unique(s$phylo)); Ai <- inverseA(tr, nodes = "TIPS")$Ainv
    cv <- if (nm == "herbivores_only") "captivity + log_rpoB + read_len_c" else cov_null
    if (nm == "non_domestic") cv <- "diet_class + log_rpoB + read_len_c"
    m <- fit_mcmc(as.formula(paste("alr_acet ~ fermentation_site +", cv)), s, Ai)
    lv <- levels(s$fermentation_site)
    out <- list(); if (all(c("foregut", "hindgut") %in% lv)) out[[1]] <- contrast(m, "foregut", "hindgut", "foregut - hindgut")
    if ("hindgut" %in% lv) out[[2]] <- contrast(m, "hindgut", "simple", "hindgut - simple")
    cbind(set = nm, n = nrow(s), species = uniqueN(s$species), rbindlist(out))
  }), fill = TRUE)
  fwrite(sens, file.path(R, "p2_sensitivity.tsv"), sep = "\t")

  # ---------------- [prereg-v1.1] additional pre-specified analyses (added before any P2 data were analysed) ----------------
  # (a) known-captivity subset and platform covariate (study-level technical confounding)
  d[, platform := fcase(grepl("NovaSeq", instrument_model), "NovaSeq", grepl("HiSeq", instrument_model), "HiSeq",
                        grepl("BGI|DNBSEQ|MGISEQ", instrument_model), "DNBSEQ", default = "other")]
  extra <- list(captivity_known = list(s = d[captivity != "unknown"], cv = cov_null),
                platform_cov = list(s = d, cv = paste(cov_null, "+ platform")))
  sens2 <- rbindlist(lapply(names(extra), function(nm) {
    s <- copy(extra[[nm]]$s); for (f in c("fermentation_site", "diet_class", "captivity")) s[, (f) := droplevels(get(f))]
    tr <- keep.tip(pt$tree, unique(s$phylo)); Ai <- inverseA(tr, nodes = "TIPS")$Ainv
    m <- fit_mcmc(as.formula(paste("alr_acet ~ fermentation_site +", extra[[nm]]$cv)), s, Ai)
    lv <- levels(s$fermentation_site); out <- list()
    if (all(c("foregut", "hindgut") %in% lv)) out[[1]] <- contrast(m, "foregut", "hindgut", "foregut - hindgut")
    if ("hindgut" %in% lv) out[[2]] <- contrast(m, "hindgut", "simple", "hindgut - simple")
    cbind(set = nm, n = nrow(s), species = uniqueN(s$species), studies = uniqueN(s$study_accession), rbindlist(out))
  }), fill = TRUE)
  fwrite(sens2, file.path(R, "p2_sensitivity_v1_1.tsv"), sep = "\t")

  # (b) within-study design check: only studies sampling >= 2 fermentation categories; study as FIXED effect,
  #     so the category contrast is estimated purely within studies (no between-study comparison).
  ms <- d[, .(nc = uniqueN(fermentation_site)), by = study_accession][nc >= 2]$study_accession
  w <- copy(d[study_accession %in% ms]); w[, fermentation_site := droplevels(fermentation_site)]
  ws <- data.table(n_studies = length(ms), n_samples = nrow(w), n_species = uniqueN(w$species))
  if (length(ms) >= 2 && nrow(w) >= 10) {
    lmw <- lmer(alr_acet ~ fermentation_site + log_rpoB + factor(study_accession) + (1 | species), data = w)
    cf <- summary(lmw)$coefficients; cf <- cf[grepl("fermentation_site", rownames(cf)), , drop = FALSE]
    ws <- cbind(ws, data.table(term = rownames(cf), estimate = cf[, 1], se = cf[, 2], t = cf[, 3]))
  }
  fwrite(ws, file.path(R, "p2_within_study.tsv"), sep = "\t")

  # (c) marsupial holdout (out-of-sample validation): fit on non-marsupial hosts; predict marsupial samples from the
  #     fixed effects only (phylogeny contributes nothing for an unseen clade); compare with observed values.
  lin <- fread("analysis/inventory/raw/host_taxa_lineage.tsv", select = c("host_tax_id", "lineage"))
  d[, marsupial := host_tax_id %in% lin[grepl("Metatheria", lineage)]$host_tax_id]
  ho <- d[marsupial == TRUE]; tr0 <- d[marsupial == FALSE]
  if (nrow(ho) >= 5 && uniqueN(ho$fermentation_site) >= 2) {
    for (f in c("fermentation_site", "diet_class", "captivity")) tr0[, (f) := droplevels(get(f))]
    Ai <- inverseA(keep.tip(pt$tree, unique(tr0$phylo)), nodes = "TIPS")$Ainv
    fx <- as.formula(paste("alr_acet ~ fermentation_site +", cov_null))
    m <- fit_mcmc(fx, tr0, Ai)
    hd <- copy(ho); for (f in c("fermentation_site", "diet_class", "captivity")) hd[, (f) := factor(as.character(get(f)), levels = levels(tr0[[f]]))]
    hd <- hd[complete.cases(hd[, .(fermentation_site, diet_class, captivity)])]
    X <- model.matrix(update(fx, NULL ~ .), data = hd)
    B <- m$Sol[, colnames(X), drop = FALSE]; P <- X %*% t(B)
    # predictive draws for a new clade/species/study: add phylo + species + study + residual variance
    Pn <- P + matrix(rnorm(length(P), 0, rep(sqrt(rowSums(m$VCV)), each = nrow(P))), nrow(P))
    hd[, `:=`(pred = rowMeans(P), pred_l95 = apply(Pn, 1, quantile, 0.025), pred_u95 = apply(Pn, 1, quantile, 0.975))]
    hs <- hd[, .(n = .N, observed = mean(alr_acet), predicted = mean(pred)), by = .(species, fermentation_site)]
    fwrite(hs, file.path(R, "p2_holdout_marsupial_species.tsv"), sep = "\t")
    fwrite(data.table(n_samples = nrow(hd), n_species = nrow(hs),
                      spearman_species = if (nrow(hs) >= 3) cor(hs$observed, hs$predicted, method = "spearman") else NA,
                      mae = mean(abs(hd$alr_acet - hd$pred)), coverage95 = mean(hd$alr_acet >= hd$pred_l95 & hd$alr_acet <= hd$pred_u95),
                      obs_hind_minus_fore = hd[fermentation_site == "hindgut", mean(alr_acet)] - hd[fermentation_site == "foregut", mean(alr_acet)],
                      pred_hind_minus_fore = hd[fermentation_site == "hindgut", mean(pred)] - hd[fermentation_site == "foregut", mean(pred)]),
           file.path(R, "p2_holdout_marsupial.tsv"), sep = "\t")
  }
  print(rbindlist(dic)); print(rbindlist(con)); print(rbindlist(sig)); print(sens2); print(ws)
}

# ============================ P3: methanogen-acetogen trade-off (proportionality) ============================
rho_p <- function(X, a, b) {          # proportionality rho on CLR-transformed composition (Erb & Notredame 2016)
  L <- log(X); C <- L - rowMeans(L)
  1 - var(C[, a] - C[, b]) / (var(C[, a]) + var(C[, b]))
}
run_p3 <- function() {
  units <- list()
  p2 <- fread(file.path(R, "p2_matrix.tsv"))[qc_pass == TRUE]
  for (s in p2[, .N, by = study_accession][N >= 10]$study_accession) units[[paste0("P2:", s)]] <- p2[study_accession == s]
  p1 <- fread(file.path(R, "p1_matrix.tsv"))[qc_pass == TRUE]
  for (st in unique(p1$study)) for (ch in c("foregut", "hindgut")) {   # within chamber: avoids compartment-driven correlation
    x <- p1[study == st & chamber == ch]; if (nrow(x) >= 10) units[[paste0("P1:", st, ":", ch)]] <- x }
  tab <- rbindlist(lapply(names(units), function(u) {
    x <- units[[u]]; X <- as.matrix(x[, paste0(MK, "_cpg"), with = FALSE]); colnames(X) <- MK
    for (j in seq_len(ncol(X))) { c0 <- hm(X[, j]); if (!is.finite(c0)) c0 <- 1e-6; X[, j] <- X[, j] + c0 }
    r <- rho_p(X, "mcrA", "acsB")
    bs <- replicate(1000, { i <- sample(nrow(X), replace = TRUE); rho_p(X[i, ], "mcrA", "acsB") })
    sp <- cor(log(X[, "mcrA"]), log(X[, "acsB"]), method = "spearman")
    data.table(unit = u, n = nrow(X), mcrA_detect = mean(x$mcrA_n >= 2), rho_p = r,
               z = atanh(pmin(pmax(r, -0.999), 0.999)), vz = var(atanh(pmin(pmax(bs, -0.999), 0.999))), spearman_logcpg = sp)
  }))
  tab <- tab[is.finite(vz) & vz > 0]
  fwrite(tab, file.path(R, "p3_units.tsv"), sep = "\t")
  r <- rma(z, vz, data = tab, method = "REML", test = "knha")
  out <- pool_tab(r, "rho_p(mcrA,acsB) Fisher-z")
  out[, `:=`(rho_pooled = tanh(estimate), rho_ci_lb = tanh(ci_lb), rho_ci_ub = tanh(ci_ub))]
  # restricted to units where mcrA is detected in >= 50 % of samples (avoids pseudocount-dominated units)
  t2 <- tab[mcrA_detect >= 0.5]
  if (nrow(t2) >= 3) { r2 <- rma(z, vz, data = t2, method = "REML", test = "knha")
    o2 <- pool_tab(r2, "rho_p, units with mcrA detected >=50%"); o2[, `:=`(rho_pooled = tanh(estimate), rho_ci_lb = tanh(ci_lb), rho_ci_ub = tanh(ci_ub))]
    out <- rbind(out, o2) }
  fwrite(out, file.path(R, "p3_pooled.tsv"), sep = "\t")
  print(out)
}

if (what %in% c("p1", "all")) run_p1()
if (what %in% c("p2", "all")) run_p2()
if (what %in% c("p3", "all")) run_p3()
writeLines(capture.output(sessionInfo()), file.path(R, "sessionInfo.txt"))
