#!/usr/bin/env Rscript
# HydroSink - primary and secondary models. Not executed in the authoring environment (no R available).
# Inputs: analysis/markers_merged.tsv (one row per sample x marker, copies_per_GE)
#         analysis/inventory.tsv (sample metadata incl. host_species, study, animal_id, site,
#                                 fermentation_site, diet_class, captivity, sample_type, read_len, depth)
#         analysis/host_tree.nwk (time-calibrated host phylogeny; TimeTree/VertLife)
suppressPackageStartupMessages({
  library(data.table); library(brms); library(ape); library(metafor); library(loo); library(propr)
})
set.seed(20261115)
dir.create("analysis/results", showWarnings = FALSE)

m   <- fread("analysis/markers_merged.tsv")
inv <- fread("analysis/inventory.tsv")
wide <- dcast(m, run_accession ~ marker, value.var = "copies_per_GE", fill = 0, fun.aggregate = sum)
d <- merge(inv, wide, by = "run_accession")
hm <- function(x) 0.5 * min(x[x > 0])
d[, alr_acet := log((acsB + hm(acsB)) / (mcrA + hm(mcrA)))]
d[, resp := rowSums(.SD), .SDcols = intersect(c("nrfA", "dsrA", "frdA", "narG", "napA"), names(d))]
d[, alr_resp := log((resp + hm(resp)) / (mcrA + hm(mcrA)))]
d[, mcrA_present := as.integer(mcrA > 0)]
d[, log_depth := log10(depth)]

# ---- Host phylogeny: correlation matrix (Brownian) ----
tr <- read.tree("analysis/host_tree.nwk")
tr <- keep.tip(tr, intersect(tr$tip.label, unique(d$host_species)))
A  <- vcv.phylo(tr, corr = TRUE)
d  <- d[host_species %in% rownames(A)]
d[, phylo := host_species]

pri <- c(prior(normal(0, 1.5), class = "b"), prior(student_t(3, 0, 1), class = "sd"))

# ---- P1: within-host compartment contrast (paired; ruminant studies) ----
p1 <- d[!is.na(animal_id) & site_class %in% c("foregut", "hindgut")]
p1[, site_class := factor(site_class, levels = c("foregut", "hindgut"))]
fit_p1 <- brm(alr_acet ~ site_class + log_depth + (1 + site_class | study) + (1 | animal_id) + (1 | host_species),
              data = p1, prior = pri, chains = 4, iter = 4000, cores = 4, backend = "cmdstanr",
              control = list(adapt_delta = 0.95), file = "analysis/results/fit_p1")
# Two-stage cross-check (mirrors random_effects.py)
es <- p1[, .(diff = alr_acet[site_class == "hindgut"] - alr_acet[site_class == "foregut"]), by = .(study, animal_id)][
  , .(yi = mean(diff), vi = var(diff) / .N, n = .N), by = study]
rma_p1 <- rma(yi, vi, data = es, method = "REML", test = "knha")
saveRDS(list(rma = rma_p1, loo = leave1out(rma_p1)), "analysis/results/p1_two_stage.rds")

# ---- P2: across hosts, fermentation site beyond phylogeny ----
f_null <- bf(alr_acet ~ diet_class + captivity + sample_type + log_depth + read_len +
               (1 | gr(phylo, cov = A)) + (1 | host_species) + (1 | study))
f_full <- bf(alr_acet ~ fermentation_site + diet_class + captivity + sample_type + log_depth + read_len +
               (1 | gr(phylo, cov = A)) + (1 | host_species) + (1 | study))
fit0 <- brm(f_null, data = d, data2 = list(A = A), prior = pri, chains = 4, iter = 4000, cores = 4,
            backend = "cmdstanr", file = "analysis/results/fit_p2_null")
fit1 <- brm(f_full, data = d, data2 = list(A = A), prior = pri, chains = 4, iter = 4000, cores = 4,
            backend = "cmdstanr", file = "analysis/results/fit_p2_full")
cmp <- loo_compare(loo(fit0), loo(fit1))
lambda <- hypothesis(fit1, "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sd_host_species__Intercept^2 + sd_study__Intercept^2 + sigma^2) = 0", class = NULL)
# Zero-inflation of mcrA: hurdle component
fit_h <- brm(bf(mcrA_present ~ fermentation_site + diet_class + log_depth + (1 | gr(phylo, cov = A)) + (1 | study)),
             family = bernoulli(), data = d, data2 = list(A = A), prior = pri, chains = 4, iter = 3000,
             cores = 4, backend = "cmdstanr", file = "analysis/results/fit_p2_hurdle")
saveRDS(list(loo_compare = cmp, lambda = lambda), "analysis/results/p2_summary.rds")

# ---- P3: trade-off (proportionality on CLR within study after host centring) ----
mk <- intersect(c("mcrA", "acsB", "nrfA", "dsrA", "frdA", "narG", "napA"), names(d))
rho_by_study <- d[, {
  X <- as.matrix(.SD) + 1e-6
  pr <- propr(X, metric = "rho", ivar = "clr")
  R <- getMatrix(pr); .(rho_mcrA_acsB = R["mcrA", "acsB"], n = .N)
}, by = study, .SDcols = mk][n >= 10]
rho_by_study[, `:=`(zi = atanh(rho_mcrA_acsB), vi = 1 / (n - 3))]
rma_p3 <- rma(zi, vi, data = rho_by_study, method = "REML", test = "knha")
saveRDS(rma_p3, "analysis/results/p3_rho_meta.rds")

# ---- Multiplicity across marker families ----
fams <- intersect(c("acsB", "nrfA", "dsrA", "frdA", "narG", "napA"), names(p1))
fam_tab <- rbindlist(lapply(fams, function(g) {
  q <- copy(p1)[, y := log((get(g) + hm(get(g))) / (mcrA + hm(mcrA)))]
  e <- q[, .(diff = y[site_class == "hindgut"] - y[site_class == "foregut"]), by = .(study, animal_id)][
    , .(yi = mean(diff), vi = var(diff) / .N), by = study]
  r <- rma(yi, vi, data = e, method = "REML", test = "knha")
  data.table(marker = g, k = r$k, est = r$b[1], ci_lb = r$ci.lb, ci_ub = r$ci.ub, p = r$pval, I2 = r$I2)
}))
fam_tab[, q_BH := p.adjust(p, "BH")]
fwrite(fam_tab, "analysis/results/p1_by_marker_family.tsv", sep = "\t")

# ---- Sensitivity analyses ----
sens <- list(
  contents_only = d[sample_type == "contents"],
  wild_only     = d[captivity == "wild"],
  no_aves       = d[host_class != "Aves"],
  depth_2M      = d[subsampled_pairs == 2e6]
)
sens_fits <- lapply(names(sens), function(nm)
  update(fit1, newdata = sens[[nm]], file = paste0("analysis/results/fit_p2_", nm), recompile = FALSE))
names(sens_fits) <- names(sens)
loso <- lapply(unique(d$study), function(s)
  update(fit1, newdata = d[study != s], file = paste0("analysis/results/loso_", s), recompile = FALSE))

# ---- Out-of-sample validation: marsupial hold-out ----
train <- d[host_order != "Diprotodontia" & host_infraclass != "Marsupialia"]
test  <- d[host_infraclass == "Marsupialia"]
fit_tr <- update(fit1, newdata = train, file = "analysis/results/fit_p2_train", recompile = FALSE)
pred <- posterior_predict(fit_tr, newdata = test, allow_new_levels = TRUE)
val <- data.table(obs = test$alr_acet, pred = colMeans(pred))
saveRDS(list(cor = cor.test(val$obs, val$pred, method = "spearman"), rmse = sqrt(mean((val$obs - val$pred)^2))),
        "analysis/results/validation_marsupial.rds")
writeLines(capture.output(sessionInfo()), "analysis/results/sessionInfo.txt")
