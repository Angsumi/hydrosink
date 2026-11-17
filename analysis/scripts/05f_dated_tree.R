#!/usr/bin/env Rscript
# HydroSink — dated-tree sensitivity for P2 (pre-planned in amendment v1.1 §5; executed after the main analyses).
# Topology: the same Open Tree of Life tree as the main analysis (host_tree.nwk, analysis tips, singletons collapsed).
# Node ages: TimeTree 5 API (api.timetree.org, MRCA of up to 10 analysis species from each child clade; NCBI taxon ids).
# Polytomies are resolved first (multi2di). A node whose TimeTree age exceeds its parent's (topology/age conflict) or
# that has no age is interpolated halfway between its parent and its oldest dated descendant; finally every node is
# forced to be at most 98% of its parent's age (ultrametric, positive branch lengths). TimeTree ages of 0 count as missing. Then the primary P2 model is refitted
# (alr_acet null vs full, contrasts, variance partition) exactly as in 05_meta_analysis.R, replacing only the tree.
# Outputs: analysis/results/dated_tree/{host_tree_dated.nwk, node_ages.tsv, p2_primary.tsv, p2_variance_partition.tsv}
commandArgs <- function(...) "none"; source("analysis/scripts/05_meta_analysis.R"); rm(commandArgs)
O <- "analysis/results/dated_tree"; dir.create(O, showWarnings = FALSE, recursive = TRUE)
d <- fread("analysis/results/p2_analysis_set.tsv")
d[, fermentation_site := relevel(factor(fermentation_site), ref = "simple")]
d[is.na(diet_class) | diet_class == "", diet_class := "unknown"]
d[, diet_class := relevel(factor(diet_class), ref = "omnivore")]; d[, captivity := relevel(factor(captivity), ref = "unknown")]
# same topology as prep_tree(), before multi2di / Grafen
tr <- read.tree("analysis/inventory/host_tree.nwk"); mp <- fread("analysis/inventory/host_ott_map.tsv"); mp[, tip := paste0("ott", ott_id)]
for (t in setdiff(mp$tip, tr$tip.label)) { hit <- grep(paste0(t, "($|[^0-9])"), tr$tip.label, value = TRUE); if (length(hit)) mp[tip == t, tip := hit[1]] }
tr$node.label <- NULL; tr <- collapse.singles(tr); tr <- keep.tip(tr, unique(d$phylo)); tr <- collapse.singles(tr)
tr <- multi2di(tr, random = FALSE); tr$edge.length <- NULL   # resolve polytomies first; every node is then dated
taxid <- d[, .(taxid = host_tax_id[1]), by = phylo]; tx <- setNames(taxid$taxid, taxid$phylo)
nt <- Ntip(tr); nodes <- (nt + 1):(nt + tr$Nnode)
desc_tips <- function(n) { if (n <= nt) return(tr$tip.label[n]); tr$tip.label[unlist(phangorn_desc(n))] }
phangorn_desc <- function(n) { out <- c(); st <- n; while (length(st)) { x <- st[1]; st <- st[-1]; ch <- tr$edge[tr$edge[, 1] == x, 2]
  out <- c(out, ch[ch <= nt]); st <- c(st, ch[ch > nt]) }; list(out) }
CF <- file.path(O, "timetree_cache.rds"); cache <- if (file.exists(CF)) readRDS(CF) else list()
age <- rep(NA_real_, length(nodes)); nmiss <- integer(length(nodes)); nq <- integer(length(nodes))
for (i in seq_along(nodes)) {
  ch <- tr$edge[tr$edge[, 1] == nodes[i], 2]
  ids <- unique(unlist(lapply(ch, function(c) head(tx[desc_tips(c)], 10))))
  url <- paste0("https://api.timetree.org/mrca/id/", paste(ids, collapse = "+"), "/json")
  key <- paste(sort(ids), collapse = "+")
  if (!is.null(cache[[key]])) txt <- cache[[key]] else {
    txt <- tryCatch(paste(readLines(url, warn = FALSE), collapse = " "), error = function(e) NA_character_); Sys.sleep(0.3)
    if (!is.na(txt)) cache[[key]] <- txt }
  fld <- function(k) { m <- regmatches(txt, regexpr(paste0('"', k, '": *\\[[^]]*\\]'), txt)); if (!length(m)) return(integer(0))
    v <- gsub("[^0-9,]", "", sub('.*\\[', "", m)); if (v == "") integer(0) else as.integer(strsplit(v, ",")[[1]]) }
  ag <- if (is.na(txt)) NA else suppressWarnings(as.numeric(sub('.*"adjusted_age": *([0-9.]+).*', "\\1", txt)))
  if (!is.na(ag) && ag > 0 && length(fld("found_ids")) >= 2) age[i] <- ag   # age 0 (e.g. same species) -> missing
  nmiss[i] <- if (is.na(txt)) length(ids) else length(fld("missing_ids")); nq[i] <- length(ids)
}
raw <- age; saveRDS(cache, CF)
parent_of <- function(n) tr$edge[tr$edge[, 2] == n, 1]
root <- nt + 1
ord <- c(root, tr$edge[, 2][tr$edge[, 2] > nt])     # preorder (ape edge matrix is in preorder after read.tree/keep.tip)
A <- setNames(age, nodes); conflict <- setNames(rep(FALSE, length(nodes)), nodes)
for (n in ord) { if (n == root) next; p <- parent_of(n)
  if (!is.na(A[as.character(n)]) && !is.na(A[as.character(p)]) && A[as.character(n)] > A[as.character(p)]) {
    A[as.character(n)] <- NA; conflict[as.character(n)] <- TRUE } }   # conflicting age dropped, then interpolated
oldest_desc <- function(n) { ch <- tr$edge[tr$edge[, 1] == n, 2]; ch <- ch[ch > nt]; if (!length(ch)) return(0)
  v <- sapply(ch, function(c) if (!is.na(A[as.character(c)])) A[as.character(c)] else oldest_desc(c)); max(v) }
interp <- setNames(rep(FALSE, length(nodes)), nodes)
for (n in ord) { if (!is.na(A[as.character(n)])) next
  p <- if (n == root) NA else A[as.character(parent_of(n))]; lo <- oldest_desc(n)
  A[as.character(n)] <- if (is.na(p)) lo * 1.1 + 1 else (p + lo) / 2; interp[as.character(n)] <- TRUE }
for (n in ord) { if (n == root) next; p <- as.character(parent_of(n))   # strictly younger than parent (<= 98% of its age)
  if (A[as.character(n)] > 0.98 * A[p]) { A[as.character(n)] <- 0.98 * A[p]; conflict[as.character(n)] <- TRUE } }
stopifnot(all(A > 0))
H <- c(rep(0, nt), A)            # node heights (Myr) indexed by node number; tips at 0 -> ultrametric
tr$edge.length <- H[tr$edge[, 1]] - H[tr$edge[, 2]]; tr2 <- tr; stopifnot(is.ultrametric(tr2), all(tr2$edge.length > 0))
fwrite(data.table(node = nodes, n_ids_queried = nq, n_ids_missing = nmiss, timetree_age = raw, used_age = A,
                  conflict_set_to_parent = conflict, interpolated = interp), file.path(O, "node_ages.tsv"), sep = "\t")
write.tree(tr2, file.path(O, "host_tree_dated.nwk"))
Ainv <- inverseA(tr2, nodes = "TIPS", scale = TRUE)$Ainv
cov_null <- "diet_class + captivity + log_rpoB + read_len_c"
m0 <- fit_mcmc(as.formula(paste("alr_acet ~", cov_null)), d, Ainv)
m1 <- fit_mcmc(as.formula(paste("alr_acet ~ fermentation_site +", cov_null)), d, Ainv)
out <- cbind(n = nrow(d), DIC_null = m0$DIC, DIC_full = m1$DIC, dDIC = m0$DIC - m1$DIC,
             rbind(contrast(m1, "foregut", "hindgut", "foregut - hindgut"), contrast(m1, "hindgut", "simple", "hindgut - simple"),
                   contrast(m1, "foregut", "simple", "foregut - simple"), contrast(m1, "avian_caecal", "simple", "avian_caecal - simple")))
fwrite(out, file.path(O, "p2_primary.tsv"), sep = "\t"); fwrite(phylo_signal(m1), file.path(O, "p2_variance_partition.tsv"), sep = "\t")
cat(sprintf("nodes %d; dated by TimeTree %d; conflicts %d; interpolated %d\n", length(nodes), sum(!is.na(raw)), sum(conflict), sum(interp)))
print(out); print(phylo_signal(m1))
