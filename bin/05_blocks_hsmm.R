#!/usr/bin/env Rscript
## Stage 05 (R) -- large-scale methylation BLOCKS for a cross-array
## longitudinal EWAS. Replaces EWASGalaxy section 2.7 (cpgCollapse with a fixed
## 500 bp gap / 1500 bp width rule, then bump hunting with a >=250 kb loess
## window) with:
##   * open-sea clusters formed where probes are close AND co-methylated, so a
##     cluster is a property of the locus rather than of the panel
##   * cluster-level effects from the within-subject model, each with its own
##     standard error
##   * a 3-state HMM along the chromosome whose transition matrix is an
##     explicit function of genomic distance, A(d) = e^{-d/L} I + (1-e^{-d/L}) 1 pi'
##   * posterior decoding, giving a soft boundary and a per-block posterior in
##     place of a thresholded loess curve
##
## Port of bin/05_blocks_hsmm.py; the numerical core is bin/ewasml.R, proved
## equivalent function by function in tests/test_equivalence.R. Baum-Welch uses
## no RNG, so unlike stage 04 this stage agrees with the Python exactly --
## including the called blocks, not merely in distribution.
##
## Inputs (from stage 01): mval.f64 + mval_dims.json, pheno_used.csv,
## probe_annotation.csv, plus crossarray_probe_map.csv.gz for island context.

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

T0 <- Sys.time()
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)),
      file = stderr())
}

here <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "ewasml.R"))

# ---------------------------------------------------------------------------
opt_list <- list(
  make_option("--in-dir", type = "character",
              help = "directory holding mval.f64/mval_dims.json, pheno_used.csv, probe_annotation.csv"),
  make_option("--probe-map", type = "character",
              help = "crossarray_probe_map.csv.gz (probe, status, island context)"),
  make_option("--out-dir", type = "character"),
  make_option("--exposure", type = "character", default = "days_on_clozapine"),
  make_option("--exposure-scale", type = "double", default = 100),
  make_option("--subject", type = "character", default = "Subject_ID"),
  make_option("--array-col", type = "character", default = "Array_Type"),
  make_option("--covars", type = "character",
              default = "smoking_score,cd8t,cd4t,bcell,gran,mono,nk",
              help = "TIME-VARYING covariates only; time-invariant ones are removed by the within transform"),
  make_option("--max-gap", type = "integer", default = 1500L,
              help = "maximum gap when forming open-sea clusters (cpgCollapse used 500 bp gap / 1500 bp width) [default %default]"),
  make_option("--rho-min", type = "double", default = 0.20,
              help = "minimum within-subject co-methylation to join two open-sea probes [default %default]"),
  make_option("--length-scale", type = "double", default = 250000,
              help = "HSMM distance length scale L in bp; the 250 kb default matches the loess window it replaces [default %default]"),
  make_option("--min-post", type = "double", default = 0.80),
  make_option("--min-clusters", type = "integer", default = 3L),
  make_option("--fixed-collapse", action = "store_true", default = FALSE,
              help = "also report the legacy fixed-width collapse, for comparison of the resulting units"),
  make_option("--var-method", type = "character", default = "limma",
              help = "variance moderation: 'limma' (squeezeVar) or 'mom', the ported method-of-moments path the Python uses [default %default]"),
  make_option("--seed", type = "integer", default = 1L)
)
opt <- parse_args(OptionParser(option_list = opt_list))
names(opt) <- gsub("-", "_", names(opt))   # optparse keeps the hyphens
for (need in c("in_dir", "probe_map", "out_dir"))
  if (is.null(opt[[need]])) stop(sprintf("--%s is required", gsub("_", "-", need)))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- inputs ---------------------------------------------------------------
mv <- read_f64(file.path(opt$in_dir, "mval"))
M <- mv$mat
ph <- read.csv(file.path(opt$in_dir, "pheno_used.csv"), stringsAsFactors = FALSE)
ph[[opt$subject]] <- as.character(ph[[opt$subject]])
ann <- read.csv(file.path(opt$in_dir, "probe_annotation.csv"),
                stringsAsFactors = FALSE)
if (!is.null(mv$colnames)) {
  key <- if ("Sample_Name" %in% names(ph)) "Sample_Name" else names(ph)[1]
  ph <- ph[match(mv$colnames, as.character(ph[[key]])), , drop = FALSE]
}
names(ann)[1] <- "probe"
if (!is.null(mv$rownames))
  ann <- ann[match(mv$rownames, as.character(ann$probe)), , drop = FALSE]
chrcol <- names(ann)[tolower(names(ann)) %in% c("chr", "seqnames", "chromosome")][1]
poscol <- names(ann)[tolower(names(ann)) %in% c("pos", "mapinfo", "start", "position")][1]
if (is.na(chrcol) || is.na(poscol))
  stop("probe_annotation.csv needs a chromosome and a position column")
names(ann)[names(ann) == chrcol] <- "chr"
names(ann)[names(ann) == poscol] <- "pos"
ann$chr <- sub("^chr", "", as.character(ann$chr))

pm <- read.csv(opt$probe_map, stringsAsFactors = FALSE)
need_cols <- c("probe", "status", "Relation_to_UCSC_CpG_Island")
if (!all(need_cols %in% names(pm)))
  stop("probe map needs columns: ", paste(need_cols, collapse = ", "))
i <- match(ann$probe, pm$probe)
ann$status <- pm$status[i]
ann$island <- pm$Relation_to_UCSC_CpG_Island[i]

# Open sea is the content stage 05 is about, and the content that differs most
# between the two arrays; restricting to shared probes keeps a cluster's
# composition independent of which array a sample came from.
isl <- ifelse(is.na(ann$island), "OpenSea", ann$island)
sel <- isl == "OpenSea" & ann$chr %in% as.character(1:22) &
  !is.na(ann$status) & ann$status == "shared_450K_EPIC"
log_msg(sprintf("open-sea, shared, autosomal probes: %d of %d harmonised",
                sum(sel), nrow(ann)))
if (sum(sel) < 1000)
  stop("too few open-sea probes survived filtering to look for blocks")
ann <- ann[sel, , drop = FALSE]
M <- M[which(sel), , drop = FALSE]
ann$chr_num <- as.integer(ann$chr)
ord <- order(ann$chr_num, ann$pos)
ann <- ann[ord, , drop = FALSE]; rownames(ann) <- NULL
M <- M[ord, , drop = FALSE]

expo <- as.double(ph[[opt$exposure]]) / opt$exposure_scale
subj <- as.character(ph[[opt$subject]])
arr <- as.character(ph[[opt$array_col]])
cov_names <- strsplit(opt$covars, ",")[[1]]
cov_names <- cov_names[nzchar(cov_names) & cov_names %in% names(ph)]
cov <- if (length(cov_names)) as.matrix(sapply(ph[cov_names], as.double)) else NULL

# ---- data-driven open-sea clusters ----------------------------------------
resid <- demean_rows_by_group(M, subj)
cl <- comethylation_clusters(ann$chr_num, as.double(ann$pos), resid,
                             max_gap = opt$max_gap, rho_min = opt$rho_min)
rm(resid); invisible(gc(FALSE))
n_cl <- if (max(cl) >= 0) as.integer(max(cl) + 1L) else 0L
log_msg(sprintf("co-methylated open-sea clusters: %d covering %d probes",
                n_cl, sum(cl >= 0)))
if (n_cl < 200)
  stop("too few open-sea clusters to fit a block model; relax --rho-min or --max-gap")

# collapse each cluster to its mean M-value profile: the cpgCollapse step, but
# over data-defined clusters
keep <- which(cl >= 0)
cid <- cl[keep]
o2 <- order(cid)                      # stable in R by default
keep <- keep[o2]; cid <- cid[o2]
Mc <- matrix(0, n_cl, ncol(M))
cchr <- character(n_cl); cstart <- numeric(n_cl)
cend <- numeric(n_cl); cnp <- integer(n_cl)
runs <- split(keep, cid)
for (k in seq_len(n_cl)) {
  rows <- runs[[as.character(k - 1L)]]
  Mc[k, ] <- colMeans(M[rows, , drop = FALSE])
  cchr[k] <- ann$chr[rows[1]]
  cstart[k] <- min(ann$pos[rows]); cend[k] <- max(ann$pos[rows])
  cnp[k] <- length(rows)
}
cmid <- (cstart + cend) / 2
o <- order(as.integer(cchr), cmid)
Mc <- Mc[o, , drop = FALSE]
cchr <- cchr[o]; cstart <- cstart[o]; cend <- cend[o]; cnp <- cnp[o]; cmid <- cmid[o]
log_msg(sprintf("cluster widths: median %d bp, median %d probes/cluster",
                as.integer(median(cend - cstart + 1)), as.integer(median(cnp))))

# ---- cluster-level within-subject effects ---------------------------------
fit <- fit_within(Mc, expo, subj, cov, var_method = opt$var_method)
log_msg(sprintf("cluster effects: median |z|=%.3f, max |z|=%.2f",
                median(abs(fit$z)), max(abs(fit$z))))

# ---- distance-aware block HSMM --------------------------------------------
hs <- fit_block_hsmm(fit$beta, fit$se, cmid, cchr,
                     length_scale = opt$length_scale, seed = opt$seed)
log_msg(sprintf("HSMM converged in %d iterations, loglik=%.1f", hs$n_iter, hs$loglik))
labs <- state_labels(hs$neutral)
log_msg(sprintf("  state means (M-value per %g units): %s", opt$exposure_scale,
                paste(sprintf("%s=%.4f", labs, hs$mu), collapse = ", ")))
log_msg(sprintf("  state sd: %s; stationary pi: %s",
                paste(round(hs$sigma, 4), collapse = ", "),
                paste(round(hs$pi, 4), collapse = ", ")))
if (hs$neutral != 2L)
  # Every cluster effect fell on one side of zero, so the state pinned at
  # mu = 0 sorted to an end. Direction is labelled relative to that column,
  # not to the middle one: with neutral at an end only one direction can be
  # called, and the middle column is a real effect state that would go
  # unreported if it were treated as no change.
  log_msg(sprintf(paste("  NOTE: the zero-effect state is column %d, not the",
                        "middle one -- all cluster effects lie on one side of",
                        "zero, so only %smethylated blocks can be called"),
                  hs$neutral, if (hs$neutral == 1L) "hyper" else "hypo"))

bl <- call_blocks(cchr, cstart, cend, hs$posterior, min_post = opt$min_post,
                  min_clusters = opt$min_clusters, neutral = hs$neutral)
log_msg(sprintf("blocks called: %d", nrow(bl)))

clus <- data.frame(chr = cchr, start = as.integer(cstart), end = as.integer(cend),
                   n_probes = cnp, effect_M = fit$beta, se = fit$se, z = fit$z,
                   stringsAsFactors = FALSE)
for (k in 1:3) clus[[paste0("post_", labs[k])]] <- hs$posterior[, k]

# ---- cross-array check on the called blocks -------------------------------
per_arr <- list()
for (a in sort(unique(arr))) {
  s <- arr == a
  # Subject count alone does not make the within fit identified: after
  # differencing out subject means, the exposure and the time-varying
  # covariates must still leave residual degrees of freedom. Same guard as
  # stage 04; without it an unidentified arm returns numbers rather than
  # declining to fit.
  df_a <- within_df(expo[s], subj[s], if (!is.null(cov)) cov[s, , drop = FALSE] else NULL)
  if (length(unique(subj[s])) < 3 || df_a <= 0) {
    log_msg(sprintf("  %s: no within-subject information to fit (df=%d); skipped",
                    a, df_a))
    next
  }
  per_arr[[a]] <- fit_within(Mc[, s, drop = FALSE], expo[s], subj[s],
                             if (!is.null(cov)) cov[s, , drop = FALSE] else NULL,
                             var_method = opt$var_method)
  clus[[paste0("effect_", a)]] <- per_arr[[a]]$beta
}
arrs <- names(per_arr)
r_pearson <- NA_real_; sign_conc <- NA_real_
if (length(arrs) == 2 && nrow(bl)) {
  ka <- arrs[1]; kb <- arrs[2]
  ea <- numeric(nrow(bl)); eb <- numeric(nrow(bl))
  for (i in seq_len(nrow(bl))) {
    m <- which(clus$chr == bl$chr[i] & clus$start >= bl$start[i] &
               clus$end <= bl$end[i])
    wa <- 1 / per_arr[[ka]]$se[m]^2
    wb <- 1 / per_arr[[kb]]$se[m]^2
    ea[i] <- sum(wa * per_arr[[ka]]$beta[m]) / sum(wa)
    eb[i] <- sum(wb * per_arr[[kb]]$beta[m]) / sum(wb)
  }
  bl[[paste0("effect_", ka)]] <- ea
  bl[[paste0("effect_", kb)]] <- eb
  ok <- is.finite(ea) & is.finite(eb)
  if (sum(ok) > 2) r_pearson <- cor(ea[ok], eb[ok])
  if (length(ea)) sign_conc <- mean(sign(ea) == sign(eb))
  log_msg(sprintf("cross-array block effects: r=%.3f, sign concordance=%.1f%%",
                  r_pearson, 100 * sign_conc))
} else if (length(arrs) < 2) {
  log_msg("cross-array comparison not possible: fewer than two arrays fitted")
}

# ---- legacy fixed-width collapse, for comparison of the units -------------
legacy <- NULL
if (opt$fixed_collapse) {
  cn <- ann$chr_num; pos <- as.double(ann$pos)
  brk <- c(TRUE, cn[-1] != cn[-length(cn)] | diff(pos) > 500)
  gid <- cumsum(brk)
  agg <- data.frame(g = gid, pos = pos)
  w <- do.call(rbind, lapply(split(agg$pos, agg$g), function(v)
    data.frame(width = max(v) - min(v) + 1, count = length(v))))
  w <- w[w$count >= 2 & w$width <= 1500, , drop = FALSE]
  # median of an empty selection is NA, which jsonlite would write as the
  # string "NA"; report null instead, matching the Python driver.
  legacy <- list(n_units = nrow(w),
                 median_width_bp = if (nrow(w)) as.numeric(median(w$width)) else NULL,
                 median_probes = if (nrow(w)) as.numeric(median(w$count)) else NULL)
  log_msg(sprintf(paste("legacy fixed-width collapse would give %d units",
                        "(median %.0f bp) vs %d co-methylation clusters (median %.0f bp)"),
                  legacy$n_units,
                  if (is.null(legacy$median_width_bp)) NA_real_ else legacy$median_width_bp,
                  n_cl,
                  median(cend - cstart + 1)))
}

# ---- outputs --------------------------------------------------------------
write.csv(bl, file.path(opt$out_dir, "blocks_hsmm.csv"), row.names = FALSE)
gz <- gzfile(file.path(opt$out_dir, "openSea_cluster_effects.csv.gz"), "w")
write.csv(clus, gz, row.names = FALSE)
close(gz)
write(jsonlite::toJSON(list(
  mu = hs$mu, sigma = hs$sigma, pi = hs$pi, length_scale = hs$length_scale,
  loglik = hs$loglik, n_iter = hs$n_iter,
  # 0-based, so this field matches the Python driver's JSON
  neutral_state = hs$neutral - 1L,
  state_labels = labs, n_clusters = n_cl, n_openSea_probes = nrow(ann),
  n_blocks = nrow(bl), cross_array_r = r_pearson,
  cross_array_sign_concordance = sign_conc, legacy_fixed_collapse = legacy,
  covars_used = cov_names, var_method = opt$var_method,
  implementation = "R", args = opt,
  runtime_s = round(as.numeric(difftime(Sys.time(), T0, units = "secs")), 1)),
  # na = "null": a skipped array leaves cross_array_r as NA, and jsonlite
  # would otherwise write the string "NA" where the Python writes null
  auto_unbox = TRUE, digits = NA, null = "null", na = "null", pretty = 1),
  file.path(opt$out_dir, "hsmm_params.json"))
log_msg("done")
