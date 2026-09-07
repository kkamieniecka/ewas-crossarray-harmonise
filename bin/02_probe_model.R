#!/usr/bin/env Rscript
# 02_probe_model.R -- probe-level longitudinal effect track.
#
# This is the stage that feeds BOTH the legacy bumphunter/blockFinder baseline
# and the ML region finders, so all methods are compared on identical input.
#
# Why the design is what it is:
#   Array type is perfectly confounded with cohort and Sentrix chip (no chip
#   carries samples from both arrays), so no additive batch correction can
#   separate "array" from "cohort". What rescues the design is that the effect
#   of interest -- change in methylation with days on clozapine -- is WITHIN
#   subject, and every subject sits entirely on one array. Fitting a
#   subject-level intercept therefore absorbs array, cohort and chip together,
#   and the exposure coefficient is estimated from within-subject contrasts
#   only. Array type is left in the model solely as a variance stratum.
#
# Outputs (all on the harmonised probe space):
#   probe_stats_joint.csv  coefficient / SE / t / p per probe, joint fit
#   probe_stats_450K.csv   same, 450K subjects only  (cross-array validation)
#   probe_stats_EPIC.csv   same, EPIC subjects only
#   resid.rds              covariate-adjusted residual matrix (co-methylation)
#
# Usage:
#   Rscript 02_probe_model.R --harmonised results/harmonised/harmonised.rds \
#     --out_dir results/model --exposure days_on_clozapine

suppressPackageStartupMessages({
  library(optparse); library(limma); library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--harmonised", type = "character"),
  make_option("--out_dir",    type = "character", default = "results/model"),
  make_option("--exposure",   type = "character", default = "days_on_clozapine"),
  make_option("--covars",     type = "character",
              default = "age,sex,smoking_score,cd8t,cd4t,bcell,mono,nk"),
  make_option("--time_scale", type = "character", default = "per100days",
              help = "raw | per100days | log1p"),
  make_option("--seed",       type = "integer", default = 20240717)
)))
set.seed(opt$seed)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

h <- readRDS(opt$harmonised)
M <- h$mval; ph <- as.data.frame(h$pheno); anno <- h$anno
stopifnot(identical(colnames(M), ph$Sample_Name))
covars <- trimws(strsplit(opt$covars, ",")[[1]])
covars <- covars[covars %in% names(ph)]

# gran is dropped: the six Houseman fractions sum to ~1, so including all six
# with an intercept is rank-deficient. gran is the largest fraction and is
# used as the reference.
covars <- setdiff(covars, "gran")

scale_time <- function(x) switch(opt$time_scale,
  raw = x, per100days = x / 100, log1p = log1p(x),
  stop("unknown time_scale"))
ph$.expo <- scale_time(as.numeric(ph[[opt$exposure]]))

fit_block <- function(Msub, phsub, label) {
  # drop covariates that are constant or all-NA within this subset
  cv <- covars[vapply(covars, function(c)
    length(unique(na.omit(phsub[[c]]))) > 1, logical(1))]
  form <- as.formula(paste("~ .expo +", paste(cv, collapse = " + ")))
  # model.matrix() silently drops incomplete rows, and its rownames are row
  # numbers rather than sample names, so the mask has to be built explicitly
  # and applied by position to keep M, the phenotypes and the design aligned.
  ok   <- which(stats::complete.cases(
    as.data.frame(phsub)[, c(".expo", cv), drop = FALSE]))
  phsub <- phsub[ok, , drop = FALSE]
  Msub <- Msub[, ok, drop = FALSE]
  des  <- model.matrix(form, data = phsub)
  stopifnot(nrow(des) == ncol(Msub))
  log(label, ": ", ncol(Msub), " samples, ",
      length(unique(phsub$Subject_ID)), " subjects, ",
      ncol(des), " design columns")

  # duplicateCorrelation estimates the within-subject (intra-block)
  # correlation and limma then fits a generalised least squares model, which
  # is the repeated-measures analogue of the per-probe regression used in
  # minfi section 2.6. Weighting by array stratum lets the two cohorts have
  # different residual variance.
  aw  <- arrayWeights(Msub, design = des)
  dc  <- duplicateCorrelation(Msub, design = des, block = phsub$Subject_ID,
                              weights = aw)
  log(label, ": consensus within-subject correlation = ",
      signif(dc$consensus.correlation, 4))
  fit <- lmFit(Msub, design = des, block = phsub$Subject_ID,
               correlation = dc$consensus.correlation, weights = aw)
  fit <- eBayes(fit, robust = TRUE)
  tt  <- topTable(fit, coef = ".expo", number = Inf, sort.by = "none")
  out <- data.table(probe = rownames(tt),
                    beta_M = tt$logFC,
                    se_M   = tt$logFC / tt$t,
                    t      = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val)
  attr(out, "consensus_cor") <- dc$consensus.correlation
  attr(out, "n_samples") <- ncol(Msub)
  attr(out, "n_subjects") <- length(unique(phsub$Subject_ID))
  list(stats = out, design = des, samples = phsub$Sample_Name, weights = aw)
}

# ---- joint fit -------------------------------------------------------------
# Array_Type enters as a fixed effect so that its (unidentifiable-from-cohort)
# mean shift is removed; the exposure coefficient remains within-subject.
ph$Array_Type <- factor(ph$Array_Type)
covars_joint <- c(covars, "Array_Type")
covars_keep <- covars
covars <- covars_joint
joint <- fit_block(M, ph, "joint")
covars <- covars_keep

res_joint <- merge(anno, joint$stats, by = "probe")
setorder(res_joint, chr, pos)
fwrite(res_joint, file.path(opt$out_dir, "probe_stats_joint.csv"))

# ---- per-array fits, for cross-array replication ---------------------------
per_arr <- list()
for (arr in levels(ph$Array_Type)) {
  idx <- which(ph$Array_Type == arr)
  f <- fit_block(M[, idx, drop = FALSE], ph[idx, , drop = FALSE], arr)
  per_arr[[arr]] <- f$stats
  o <- merge(anno, f$stats, by = "probe"); setorder(o, chr, pos)
  fwrite(o, file.path(opt$out_dir, paste0("probe_stats_", arr, ".csv")))
}

# ---- residual matrix for data-driven co-methylation clustering -------------
# Residuals are taken from the covariate-only model (exposure retained), so the
# correlation structure reflects genuine co-methylation plus shared technical
# variation, not the effect being tested.
des_cov <- joint$design[, setdiff(colnames(joint$design), ".expo"), drop = FALSE]
Mk <- M[, joint$samples, drop = FALSE]
# lm.fit()$coefficients is (n_covariates x n_probes), so the fitted values are
# des_cov %*% coefficients = (n_samples x n_probes) and must be transposed back
# to the probes-by-samples orientation of M.
fitted_cov <- des_cov %*% lm.fit(des_cov, t(Mk))$coefficients
resid <- Mk - t(fitted_cov)
stopifnot(dim(resid) == dim(Mk))
saveRDS(list(resid = resid, samples = joint$samples,
             subject = ph$Subject_ID[match(joint$samples, ph$Sample_Name)],
             array   = as.character(ph$Array_Type[match(joint$samples, ph$Sample_Name)])),
        file.path(opt$out_dir, "resid.rds"))

meta <- data.table(
  fit = c("joint", levels(ph$Array_Type)),
  n_samples  = c(attr(joint$stats, "n_samples"),
                 vapply(per_arr, function(x) attr(x, "n_samples"), numeric(1))),
  n_subjects = c(attr(joint$stats, "n_subjects"),
                 vapply(per_arr, function(x) attr(x, "n_subjects"), numeric(1))),
  within_subject_cor = c(attr(joint$stats, "consensus_cor"),
                 vapply(per_arr, function(x) attr(x, "consensus_cor"), numeric(1))),
  n_probes = c(nrow(joint$stats), vapply(per_arr, nrow, numeric(1))),
  n_fdr05  = c(sum(joint$stats$fdr < 0.05),
               vapply(per_arr, function(x) sum(x$fdr < 0.05), numeric(1))))
fwrite(meta, file.path(opt$out_dir, "model_summary.csv"))
print(meta)
log("done")
