#!/usr/bin/env Rscript
# 03_baseline_bumphunter.R -- the LEGACY baseline, run on exactly the same
# harmonised matrix as the replacement, so the comparison is not confounded by
# preprocessing.
#
# This reproduces what the EWASGalaxy suite does today, i.e. sections 2.6 and
# 2.7 of Aryee et al. 2014 as minfi implements them:
#
#   2.6  clusterMaker(maxGap = 300) -> bumphunter(smooth = TRUE, loess,
#        B permutations of the design column of interest)
#   2.7  cpgCollapse(open-sea, 500 bp gap / 1500 bp width) -> blockFinder
#        (bump hunting with a 250 kb loess window)
#
# Two deliberate properties of this baseline, because they ARE the problem:
#   * the design matrix has no subject term, so between-subject variation --
#     which here is also between-array, between-cohort and between-chip
#     variation -- lands in the residual;
#   * bumphunter's permutation reassigns the exposure column freely across all
#     samples, which breaks the longitudinal pairing.
# Both are what the replacement fixes; this script quantifies the cost.
#
# Usage
#   Rscript 03_baseline_bumphunter.R --in-dir results/harmonised \
#       --out-dir results/baseline --n-perm 200

suppressPackageStartupMessages({
  library(optparse); library(minfi); library(bumphunter)
  library(data.table); library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--in-dir", type = "character"),
  make_option("--out-dir", type = "character"),
  make_option("--exposure", type = "character", default = "days_on_clozapine"),
  make_option("--exposure-scale", type = "double", default = 100),
  make_option("--covars", type = "character",
              default = "age,sex,smoking_score,cd8t,cd4t,bcell,gran,mono,nk"),
  make_option("--adjust-array", action = "store_true", default = FALSE,
              help = "add Array_Type to the design; note it is collinear with cohort and chip"),
  make_option("--max-gap", type = "integer", default = 300L),
  make_option("--cutoff", type = "double", default = 0.10),
  # blockFinder operates on collapsed open-sea units, whose effect sizes are
  # much smaller than probe-level ones, so it needs its own cutoff. At the
  # probe-level default, bumphunter finds no bump inside blockFinder and minfi
  # then fails on an atomic result table rather than returning an empty one.
  make_option("--block-cutoff", type = "double", default = 0.02),
  make_option("--n-perm", type = "integer", default = 200L),
  make_option("--block-window", type = "double", default = 250000),
  make_option("--cores", type = "integer", default = 4L),
  make_option("--seed", type = "integer", default = 1L),
  make_option("--skip-blocks", action = "store_true", default = FALSE)
)))
dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)
t0 <- Sys.time()
log <- function(...) cat(sprintf("[%6.1fs] ", as.numeric(difftime(Sys.time(), t0, units = "secs"))), ..., "\n", sep = "")

h <- readRDS(file.path(opt$`in-dir`, "harmonised.rds"))
M <- h$mval; pheno <- h$pheno; ann <- h$anno
stopifnot(identical(colnames(M), as.character(pheno$Sample_Name)))

# autosomes, genomic order -- matches the python stage's probe universe
ann$chr_s <- sub("^chr", "", as.character(ann$chr))
keep <- ann$chr_s %in% as.character(1:22)
M <- M[keep, , drop = FALSE]; ann <- ann[keep, , drop = FALSE]
o <- order(as.integer(ann$chr_s), ann$pos)
M <- M[o, , drop = FALSE]; ann <- ann[o, , drop = FALSE]
log("baseline universe: ", nrow(M), " autosomal probes x ", ncol(M), " samples")

expo <- as.numeric(pheno[[opt$exposure]]) / opt$`exposure-scale`
cv <- strsplit(opt$covars, ",")[[1]]
cv <- cv[cv %in% names(pheno)]
form <- paste("~ expo", if (length(cv)) paste("+", paste(cv, collapse = " + ")) else "",
              if (opt$`adjust-array`) "+ Array_Type" else "")
design <- model.matrix(as.formula(form), data = cbind(pheno, expo = expo))
log("design: ", form, "  ->  ", ncol(design), " columns, rank ", qr(design)$rank)
if (qr(design)$rank < ncol(design))
  log("!! design is rank deficient -- this is the array/cohort/chip collinearity")

cl <- clusterMaker(chr = ann$chr_s, pos = ann$pos, maxGap = opt$`max-gap`)
log("clusterMaker(maxGap=", opt$`max-gap`, "): ", length(unique(cl)), " clusters; ",
    "median probes/cluster = ", median(table(cl)))

log("bumphunter, B=", opt$`n-perm`, " (free permutation of the exposure column)")
bh <- bumphunter(M, design = design, chr = ann$chr_s, pos = ann$pos, cluster = cl,
                 coef = 2, cutoff = opt$cutoff, B = opt$`n-perm`,
                 smooth = TRUE, smoothFunction = loessByCluster,
                 nullMethod = "permutation", verbose = FALSE)
tab <- as.data.table(bh$table)
log("bumphunter regions: ", nrow(tab), "; at fwer<=0.05: ",
    sum(tab$fwer <= 0.05, na.rm = TRUE))
fwrite(tab, file.path(opt$`out-dir`, "bumphunter_regions.csv"))

# per-array refit, so the same cross-array generalisation metric can be
# computed for the baseline as for the replacement
per_arr <- list()
for (a in unique(pheno$Array_Type)) {
  s <- pheno$Array_Type == a
  if (sum(s) < 6) next
  cvs <- cv[vapply(cv, function(x) length(unique(pheno[[x]][s])) > 1, logical(1))]
  f2 <- paste("~ expo", if (length(cvs)) paste("+", paste(cvs, collapse = " + ")) else "")
  d2 <- model.matrix(as.formula(f2), data = cbind(pheno[s, , drop = FALSE], expo = expo[s]))
  fit <- limma::lmFit(M[, s, drop = FALSE], d2)
  per_arr[[a]] <- fit$coefficients[, 2]
  log("per-array fit ", a, ": n=", sum(s))
}
if (length(per_arr) == 2 && nrow(tab)) {
  ea <- eb <- numeric(nrow(tab))
  for (i in seq_len(nrow(tab))) {
    idx <- which(ann$chr_s == sub("^chr", "", tab$chr[i]) &
                 ann$pos >= tab$start[i] & ann$pos <= tab$end[i])
    ea[i] <- mean(per_arr[[1]][idx]); eb[i] <- mean(per_arr[[2]][idx])
  }
  tab[[paste0("effect_", names(per_arr)[1])]] <- ea
  tab[[paste0("effect_", names(per_arr)[2])]] <- eb
  ok <- is.finite(ea) & is.finite(eb)
  r_pear <- if (sum(ok) > 2) cor(ea[ok], eb[ok]) else NA_real_
  sconc <- if (any(ok)) mean(sign(ea[ok]) == sign(eb[ok])) else NA_real_
  log("cross-array region effects: r=", round(r_pear, 3),
      ", sign concordance=", round(100 * sconc, 1), "%")
  fwrite(tab, file.path(opt$`out-dir`, "bumphunter_regions.csv"))
} else { r_pear <- NA_real_; sconc <- NA_real_ }

blk_n <- NA_integer_
if (!opt$`skip-blocks`) {
  log("cpgCollapse (open sea, 500bp gap / 1500bp width) + blockFinder")
  gr <- try({
    gset <- h$gset
    if (is.null(gset)) stop("harmonised.rds carries no GenomicRatioSet; rerun 01 with --keep-gset")
    coll <- cpgCollapse(gset, what = "Beta", maxGap = 500, blockMaxGap = 250000,
                        verbose = FALSE)
    dsn <- design[, seq_len(min(2, ncol(design))), drop = FALSE]
    bf <- blockFinder(coll$object, design = dsn, coef = 2,
                      cutoff = opt$`block-cutoff`, B = min(opt$`n-perm`, 100L),
                      what = "Beta", smooth = TRUE, verbose = FALSE)
    # blockFinder returns a non-data.frame $table (NA) when no block clears the
    # cutoff, so this cannot be de-referenced blindly.
    bt <- if (is.data.frame(bf$table)) as.data.table(bf$table) else
      data.table(chr = character(), start = integer(), end = integer(),
                 value = numeric(), area = numeric(), L = integer(),
                 fwer = numeric())
    fwrite(bt, file.path(opt$`out-dir`, "blockfinder_blocks.csv"))
    log("blockFinder blocks: ", nrow(bt), "; collapsed units: ",
        nrow(coll$object), " (median width ",
        median(width(granges(coll$object))), " bp)")
    nrow(bt)
  }, silent = TRUE)
  if (inherits(gr, "try-error")) log("!! block stage failed: ", as.character(gr))
  else blk_n <- gr
}

write_json(list(
  n_probes = nrow(M), n_samples = ncol(M),
  design_formula = form, design_rank = qr(design)$rank,
  design_cols = ncol(design),
  n_clusters = length(unique(cl)),
  n_regions = nrow(tab),
  n_regions_fwer05 = sum(tab$fwer <= 0.05, na.rm = TRUE),
  cross_array_r = r_pear, cross_array_sign_concordance = sconc,
  n_blocks = blk_n, n_perm = opt$`n-perm`,
  runtime_s = as.numeric(difftime(Sys.time(), t0, units = "secs"))),
  file.path(opt$`out-dir`, "baseline_summary.json"), auto_unbox = TRUE, digits = 6)
log("done")
