#!/usr/bin/env Rscript
## Stage 04 (R) -- replaces EWASGalaxy section 2.6 (clusterMaker + loess +
## bumphunter permutation) with:
##   * data-driven co-methylation clustering
##   * weighted total-variation denoising with a distance-decayed penalty
##   * within-subject permutation FWER
##   * subject-level stability selection
##   * per-array replication of every called region
##
## Port of bin/04_dmr_ml.py; the numerical core is bin/ewasml.R, proved
## equivalent function by function in tests/test_equivalence.R. The two
## implementations differ only where the language forces it: R's RNG produces
## different permutation draws, so FWER p-values and stability frequencies
## agree in distribution, not to the digit.
##
## Inputs (from stage 01): mval.f64 + mval_dims.json, pheno_used.csv,
## probe_annotation.csv.

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
  make_option("--out-dir", type = "character"),
  make_option("--exposure", type = "character", default = "days_on_clozapine"),
  make_option("--exposure-scale", type = "double", default = 100,
              help = "divide the exposure by this so effects read per unit [default %default]"),
  make_option("--subject", type = "character", default = "Subject_ID"),
  make_option("--array-col", type = "character", default = "Array_Type"),
  make_option("--covars", type = "character",
              default = "smoking_score,cd8t,cd4t,bcell,gran,mono,nk",
              help = "TIME-VARYING covariates only; time-invariant ones are removed by the within transform"),
  make_option("--max-gap", type = "integer", default = 1000L),
  make_option("--rho-min", type = "double", default = 0.30),
  make_option("--decay-bp", type = "double", default = 1000),
  make_option("--lam-grid", type = "character", default = "0.05,0.1,0.25,0.5,1.0,2.0,4.0"),
  make_option("--cv-folds", type = "integer", default = 5L),
  make_option("--min-probes", type = "integer", default = 3L),
  make_option("--min-effect", type = "double", default = 0.05,
              help = "minimum |region effect| in M-value units [default %default]"),
  make_option("--n-perm", type = "integer", default = 200L),
  make_option("--n-boot", type = "integer", default = 50L),
  make_option("--perm-iter", type = "integer", default = 800L,
              help = "ADMM iterations for permutation replicates [default %default]"),
  make_option("--seed", type = "integer", default = 1L),
  make_option("--also-naive-perm", action = "store_true", default = FALSE,
              help = "additionally run the bumphunter-style free permutation, to quantify its mis-calibration"),
  make_option("--var-method", type = "character", default = "limma",
              help = "variance moderation: 'limma' (squeezeVar) or 'mom' (the ported method-of-moments) [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
# The flags are dashed to match 04_dmr_ml.py exactly, so either script can be
# dropped into the same Nextflow process; optparse keeps the dash in the name.
names(opt) <- gsub("-", "_", names(opt))
for (req in c("in_dir", "out_dir"))
  if (is.null(opt[[req]])) stop("--", gsub("_", "-", req), " is required")
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)

# ---------------------------------------------------------------------------
load_inputs <- function(opt) {
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
  if (!is.null(mv$rownames)) {
    ann <- ann[match(mv$rownames, as.character(ann[[1]])), , drop = FALSE]
  }
  names(ann)[1] <- "probe"
  chrcol <- names(ann)[tolower(names(ann)) %in% c("chr", "seqnames", "chromosome")][1]
  poscol <- names(ann)[tolower(names(ann)) %in% c("pos", "mapinfo", "start", "position")][1]
  if (is.na(chrcol) || is.na(poscol))
    stop("probe_annotation.csv needs a chromosome and a position column")
  names(ann)[names(ann) == chrcol] <- "chr"
  names(ann)[names(ann) == poscol] <- "pos"
  # genomic order, autosomes only (sex chromosomes need a sex-stratified model)
  ann$chr <- sub("^chr", "", as.character(ann$chr))
  aut <- ann$chr %in% as.character(1:22)
  ann <- ann[aut, , drop = FALSE]
  M <- M[which(aut), , drop = FALSE]
  ann$chr_num <- as.integer(ann$chr)
  ord <- order(ann$chr_num, ann$pos)
  ann <- ann[ord, , drop = FALSE]
  rownames(ann) <- NULL
  M <- M[ord, , drop = FALSE]
  log_msg(sprintf("loaded %d autosomal probes x %d samples", nrow(M), ncol(M)))
  list(M = M, ann = ann, ph = ph)
}

build_design <- function(ph, opt) {
  expo <- as.double(ph[[opt$exposure]]) / opt$exposure_scale
  subj <- as.character(ph[[opt$subject]])
  arr <- as.character(ph[[opt$array_col]])
  cov_names <- strsplit(opt$covars, ",")[[1]]
  cov_names <- cov_names[nzchar(cov_names) & cov_names %in% names(ph)]
  cov <- if (length(cov_names))
    as.matrix(sapply(ph[cov_names], as.double)) else NULL
  list(expo = expo, subj = subj, arr = arr, cov = cov, cov_names = cov_names)
}

#' Chromosome-wise TV denoising of the precision-weighted effect track.
tv_track <- function(ann, beta, se, lam0, decay_bp, n_iter) {
  w <- ifelse(se > 0, 1 / se^2, 0)
  theta <- numeric(length(beta))
  for (c in sort(unique(ann$chr_num))) {
    idx <- which(ann$chr_num == c)
    lam <- distance_penalty(as.numeric(ann$pos[idx]), lam0 = lam0,
                            decay_bp = decay_bp, w = w[idx])
    theta[idx] <- tv_denoise(beta[idx], w[idx], lam, n_iter = n_iter)
  }
  list(theta = theta, w = w)
}

#' Choose the smoothness by leaving out whole SUBJECTS.
#'
#' Held-out subjects (not held-out samples) are required: a subject's visits
#' are the unit the within-subject effect is estimated from, so splitting them
#' across train and test would leak the quantity being validated. Score =
#' precision-weighted squared error of the training track against the held-out
#' subjects' own within-subject effects.
cv_lambda <- function(M, ann, expo, subj, cov, opt) {
  uniq <- sort(unique(subj))
  # deterministic, language-independent split: see fold_assign in ewasml.R
  folds <- fold_assign(subj, opt$cv_folds, seed = opt$seed)
  grid <- as.numeric(strsplit(opt$lam_grid, ",")[[1]])
  rows <- list()
  n_skipped <- 0L
  for (lam0 in grid) {
    tot <- 0; ntot <- 0L
    for (f in folds) {
      te <- subj %in% f
      tr <- !te
      if (length(unique(subj[te])) < 2 || length(unique(subj[tr])) < 4) {
        n_skipped <<- n_skipped + 1L; next
      }
      # Subject counts alone do not make a fold estimable: a fold of mostly
      # single-visit subjects has no residual df once subject intercepts are
      # removed. Skip it rather than fail the stage.
      dfte <- within_df(expo[te], subj[te],
                        if (!is.null(cov)) cov[te, , drop = FALSE] else NULL)
      dftr <- within_df(expo[tr], subj[tr],
                        if (!is.null(cov)) cov[tr, , drop = FALSE] else NULL)
      if (dfte <= 0 || dftr <= 0) { n_skipped <<- n_skipped + 1L; next }
      ftr <- fit_within(M[, tr, drop = FALSE], expo[tr], subj[tr],
                        if (!is.null(cov)) cov[tr, , drop = FALSE] else NULL,
                        var_method = opt$var_method)
      fte <- fit_within(M[, te, drop = FALSE], expo[te], subj[te],
                        if (!is.null(cov)) cov[te, , drop = FALSE] else NULL,
                        var_method = opt$var_method)
      th <- tv_track(ann, ftr$beta, ftr$se, lam0, opt$decay_bp, opt$perm_iter)$theta
      wte <- ifelse(fte$se > 0, 1 / fte$se^2, 0)
      tot <- tot + sum(wte * (fte$beta - th)^2)
      ntot <- ntot + sum(wte > 0)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      lam0 = lam0, wsse = tot, n = ntot,
      score = if (ntot > 0) tot / ntot else NA_real_)
    log_msg(sprintf("  CV lam0=%-6g score=%.6g", lam0,
                    rows[[length(rows)]]$score))
  }
  cv <- do.call(rbind, rows)
  if (n_skipped)
    log_msg(sprintf("  %d/%d fold-evaluations skipped as not identified ",
                    n_skipped, length(grid) * length(folds)),
            "(too few subjects or no residual df)")
  if (!any(is.finite(cv$score)))
    stop(sprintf(paste("no cross-validation fold was estimable: with %d subjects",
                       "and %d folds every split left too little within-subject",
                       "variation. Use fewer folds, or drop single-visit subjects."),
                 length(uniq), opt$cv_folds))
  list(best = cv$lam0[which.min(cv$score)], cv = cv)
}

# ---------------------------------------------------------------------------
inp <- load_inputs(opt)
M <- inp$M; ann <- inp$ann; ph <- inp$ph
des <- build_design(ph, opt)
expo <- des$expo; subj <- des$subj; arr <- des$arr; cov <- des$cov

n_vis <- table(subj)
log_msg(sprintf("%d subjects, visits/subject min=%d median=%d max=%d",
                length(n_vis), min(n_vis), as.integer(median(n_vis)), max(n_vis)))
log_msg("array split: ", paste(names(table(arr)), table(arr), sep = "=",
                               collapse = ", "))
if (any(n_vis < 2))
  log_msg(sprintf("note: %d subject(s) with a single visit contribute nothing ",
                  sum(n_vis < 2)),
          "to a within-subject contrast and are carried but not informative")

## --- 1. observed probe-level fit -------------------------------------------
fit <- fit_within(M, expo, subj, cov, var_method = opt$var_method)
log_msg(sprintf("probe fit: df=%d, median |z|=%.3f, max |z|=%.2f",
                fit$df, median(abs(fit$z)), max(abs(fit$z))))

## --- 2. co-methylation clusters --------------------------------------------
resid <- demean_rows_by_group(M, subj)
Xw <- demean_by_group(cbind(expo, if (!is.null(cov)) cov), subj)
Xw <- Xw[, apply(abs(Xw), 2, max) > 1e-9, drop = FALSE]
resid <- resid - (resid %*% Xw %*% pinv(crossprod(Xw))) %*% t(Xw)
cl <- comethylation_clusters(ann$chr_num, as.numeric(ann$pos), resid,
                             max_gap = opt$max_gap, rho_min = opt$rho_min)
n_cl <- if (max(cl) >= 0) max(cl) + 1L else 0L
log_msg(sprintf("co-methylation clusters: %d clusters covering %d probes (%.1f%%)",
                n_cl, sum(cl >= 0), 100 * mean(cl >= 0)))
rm(resid); invisible(gc())

## --- 3. smoothness by held-out-subject CV ----------------------------------
log_msg("cross-validating the TV smoothness over held-out subjects")
cvr <- cv_lambda(M, ann, expo, subj, cov, opt)
lam0 <- cvr$best
write.csv(cvr$cv, file.path(opt$out_dir, "lambda_cv.csv"), row.names = FALSE)
log_msg(sprintf("selected lam0=%g", lam0))

## --- 4. observed regions ---------------------------------------------------
tt <- tv_track(ann, fit$beta, fit$se, lam0, opt$decay_bp, n_iter = 3000L)
theta <- tt$theta
cl_arg <- ifelse(cl >= 0, cl, -seq_along(cl))
obs <- call_regions(ann$chr, ann$pos, fit$beta, fit$se, theta,
                    min_probes = opt$min_probes, min_effect = opt$min_effect,
                    cluster = cl_arg)
reg <- obs$table
log_msg(sprintf("observed regions: %d; max |z| = %.2f", nrow(reg), obs$max_abs_z))
if (!nrow(reg)) {
  write.csv(reg, file.path(opt$out_dir, "dmr_ml.csv"), row.names = FALSE)
  write_json(list(n_regions = 0L, lam0 = lam0),
             file.path(opt$out_dir, "run_config.json"), auto_unbox = TRUE,
             pretty = TRUE)
  log_msg("no regions passed the effect-size floor; stopping")
  quit(status = 0L)
}

## --- 5. permutation FWER ---------------------------------------------------
null_max <- function(scheme, n_perm) {
  if (n_perm <= 0) return(numeric(0))
  fn <- if (scheme == "within") within_subject_permutation else naive_permutation
  out <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    pe <- fn(expo, subj)
    f <- fit_within(M, pe, subj, cov, var_method = opt$var_method)
    th <- tv_track(ann, f$beta, f$se, lam0, opt$decay_bp, opt$perm_iter)$theta
    r <- call_regions(ann$chr, ann$pos, f$beta, f$se, th,
                      min_probes = opt$min_probes, min_effect = opt$min_effect)
    out[b] <- r$max_abs_z
    if (b %% 25 == 0)
      log_msg(sprintf("  %s permutation %d/%d, running max|z| 95th pct = %.2f",
                      scheme, b, n_perm, quantile(out[seq_len(b)], 0.95)))
  }
  out
}

log_msg(sprintf("within-subject permutation null, %d replicates", opt$n_perm))
null_w <- null_max("within", opt$n_perm)
reg$p_fwer_within <- vapply(reg$z, function(z)
  (1 + sum(null_w >= abs(z))) / (1 + length(null_w)), numeric(1))
# Monte Carlo standard error of the permutation p-value. Reported because the
# FWER decision at 0.05 is a step function of the null tail: when many
# correlated regions have similar |z|, a p-value whose MCSE straddles the
# threshold is not evidence about the region, only about n_perm.
mcse <- function(p, n_perm) if (n_perm > 0) sqrt(p * (1 - p) / n_perm) else 0 * p
reg$p_fwer_within_mcse <- mcse(reg$p_fwer_within, length(null_w))
n_borderline <- sum(reg$p_fwer_within - 2 * reg$p_fwer_within_mcse < 0.05 &
                    reg$p_fwer_within + 2 * reg$p_fwer_within_mcse > 0.05)
if (n_borderline)
  log_msg(sprintf(paste("%d region(s) have a FWER p-value within 2 MC standard",
                        "errors of 0.05; raise --n-perm above %d to resolve them"),
                  n_borderline, opt$n_perm))
nulls <- data.frame(within = null_w)
if (opt$also_naive_perm) {
  log_msg(sprintf("naive (bumphunter-style) permutation null, %d replicates",
                  opt$n_perm))
  null_n <- null_max("naive", opt$n_perm)
  reg$p_fwer_naive <- vapply(reg$z, function(z)
    (1 + sum(null_n >= abs(z))) / (1 + length(null_n)), numeric(1))
  reg$p_fwer_naive_mcse <- mcse(reg$p_fwer_naive, length(null_n))
  nulls$naive <- null_n
  if (length(null_w) && length(null_n))
    log_msg(sprintf("null max|z| 95th pct: within=%.2f naive=%.2f",
                    quantile(null_w, 0.95), quantile(null_n, 0.95)))
}
write.csv(nulls, file.path(opt$out_dir, "permutation_null.csv"), row.names = FALSE)

## --- 6. stability selection over subject subsamples ------------------------
log_msg(sprintf("stability selection, %d subject subsamples", opt$n_boot))
keys <- sprintf("%s:%d-%d", reg$chr, reg$start, reg$end)
fit_subset <- function(sub) {
  sel <- subj %in% sub
  if (length(unique(subj[sel])) < 3 ||
      within_df(expo[sel], subj[sel],
                if (!is.null(cov)) cov[sel, , drop = FALSE] else NULL) <= 0)
    return(character(0))
  f <- fit_within(M[, sel, drop = FALSE], expo[sel], subj[sel],
                  if (!is.null(cov)) cov[sel, , drop = FALSE] else NULL,
                  var_method = opt$var_method)
  th <- tv_track(ann, f$beta, f$se, lam0, opt$decay_bp, opt$perm_iter)$theta
  r <- call_regions(ann$chr, ann$pos, f$beta, f$se, th,
                    min_probes = opt$min_probes, min_effect = opt$min_effect)
  if (!nrow(r$table)) return(character(0))
  found <- character(0)
  for (i in seq_len(nrow(reg))) {
    # count as re-discovered when the two intervals overlap
    hit <- r$table$chr == reg$chr[i] & r$table$start <= reg$end[i] &
      reg$start[i] <= r$table$end
    if (any(hit)) found <- c(found, keys[i])
  }
  found
}
ss <- stability_selection(fit_subset, subj, n_boot = opt$n_boot, frac = 0.5,
                          seed = opt$seed + 7L)
reg$stability <- ifelse(keys %in% names(ss$freq), ss$freq[keys], 0)
reg$stability[is.na(reg$stability)] <- 0
log_msg(sprintf("stability: median=%.2f, %d regions selected in >=60%% of subsamples",
                median(reg$stability), sum(reg$stability >= 0.6)))

## --- 7. cross-array generalisation -----------------------------------------
log_msg("cross-array generalisation")
r_pearson <- NA_real_; sign_conc <- NA_real_
arrs_used <- character(0)
for (a in sort(unique(arr))) {
  sel <- arr == a
  if (length(unique(subj[sel])) < 3 ||
      within_df(expo[sel], subj[sel],
                if (!is.null(cov)) cov[sel, , drop = FALSE] else NULL) <= 0) {
    log_msg(sprintf("  %s: not enough within-subject information to fit; skipped", a))
    next
  }
  f <- fit_within(M[, sel, drop = FALSE], expo[sel], subj[sel],
                  if (!is.null(cov)) cov[sel, , drop = FALSE] else NULL,
                  var_method = opt$var_method)
  pe <- numeric(nrow(reg)); se_ <- numeric(nrow(reg))
  for (i in seq_len(nrow(reg))) {
    idx <- (reg$probe_start[i] + 1L):(reg$probe_end[i] + 1L)   # stored 0-based
    ws <- ifelse(f$se[idx] > 0, 1 / f$se[idx]^2, 0)
    pe[i] <- if (sum(ws) > 0) sum(ws * f$beta[idx]) / sum(ws) else NA_real_
    se_[i] <- if (sum(ws) > 0) sqrt(1 / sum(ws)) else NA_real_
  }
  reg[[paste0("effect_", a)]] <- pe
  reg[[paste0("se_", a)]] <- se_
  arrs_used <- c(arrs_used, a)
}
if (length(arrs_used) == 2) {
  x <- reg[[paste0("effect_", arrs_used[1])]]
  y <- reg[[paste0("effect_", arrs_used[2])]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) > 2) r_pearson <- cor(x[ok], y[ok])
  if (any(ok)) sign_conc <- mean(sign(x[ok]) == sign(y[ok]))
  reg$sign_concordant <- is.finite(x) & is.finite(y) & sign(x) == sign(y)
  log_msg(sprintf("cross-array region effects: r=%.3f, sign concordance=%.1f%% over %d regions",
                  r_pearson, 100 * sign_conc, sum(ok)))
}

reg <- reg[order(reg$p_fwer_within), , drop = FALSE]
write.csv(reg, file.path(opt$out_dir, "dmr_ml.csv"), row.names = FALSE)

cfg <- opt
cfg$help <- NULL
cfg <- c(cfg, list(
  implementation = "R", ewasml_r = TRUE,
  n_probes = nrow(M), n_samples = ncol(M), n_subjects = length(unique(subj)),
  covars_used = des$cov_names, lam0_selected = lam0, n_clusters = n_cl,
  n_regions = nrow(reg), n_fwer_05 = sum(reg$p_fwer_within <= 0.05),
  n_fwer_borderline = n_borderline,
  cross_array_r = r_pearson, cross_array_sign_concordance = sign_conc,
  # --n-perm 0 is a legitimate mode (smoothness diagnostics without the
  # expensive null), so the null may be empty
  null_within_q95 = if (length(null_w)) quantile(null_w, 0.95, names = FALSE) else NULL,
  probe_fit_df = fit$df, seed = opt$seed, var_method = opt$var_method,
  runtime_s = round(as.numeric(difftime(Sys.time(), T0, units = "secs")), 1)))
if (opt$also_naive_perm && !is.null(nulls$naive) && length(nulls$naive))
  cfg$null_naive_q95 <- quantile(nulls$naive, 0.95, names = FALSE)
write_json(cfg, file.path(opt$out_dir, "run_config.json"), auto_unbox = TRUE,
           pretty = TRUE, na = "null")
log_msg(sprintf("wrote %d regions; %d at FWER<=0.05", nrow(reg),
                sum(reg$p_fwer_within <= 0.05)))
