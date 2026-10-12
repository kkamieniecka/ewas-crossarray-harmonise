# Class entry points and the analysis functions users call.
#
# This file is NOT generated from bin/ewasml.R, and it is the only part of the
# package that is not. The core in bin/ewasml.R is base R plus limma so that
# the pipeline's stage drivers can source it directly, matrix in and
# data.frame out; everything that knows about Bioconductor classes lives here.
# tests/test_pkg_identity.R skips this file by name.
#
# SummarizedExperiment, GenomicRanges, S4Vectors and IRanges are Imports: the
# object interface is the package's interface, and a GenomicRatioSet is what a
# methylation user arrives with. The pipeline keeps its dependency-free path
# by sourcing bin/ewasml.R rather than loading the package.
#
# ewas_regions() and ewas_blocks() are the two high-level functions. They do
# what bin/04_dmr_ml.R and bin/05_blocks_hsmm.R do: assay and design
# extraction, within-subject residuals, co-methylation clustering, smoothness
# selection by held-out subjects, region or block calling, the permutation
# null and subject-level stability selection. The matrix functions remain
# exported as building blocks.
#
# Edit this file in pkg/R-src/; tools/build_pkg.py copies it into the
# generated tree, where it will be overwritten.


#' Column of `colData` as a model vector
#'
#' Resolves an argument that may either name a `colData` column or be the
#' vector itself, and fails with a message naming the available columns.
#'
#' @param x A character scalar naming a column, or a vector of length
#' `nrow(cd)`.
#' @param cd The `colData` of the object, as a data frame.
#' @param what Name of the argument, used in error messages.
#' @return The resolved vector.
#' @keywords internal
#' @noRd
se_column <- function(x, cd, what) {
  if (is.character(x) && length(x) == 1L && !(x %in% names(cd)) &&
      nrow(cd) != 1L)
    stop(what, " = \"", x, "\" is not a colData column; available: ",
         paste(names(cd), collapse = ", "), call. = FALSE)
  if (is.character(x) && length(x) == 1L && x %in% names(cd)) x <- cd[[x]]
  if (length(x) != nrow(cd))
    stop(what, " has length ", length(x), " but the object has ", nrow(cd),
         " samples", call. = FALSE)
  x
}


#' Coerce an exposure to a numeric contrast
#'
#' Numeric, integer and logical exposures pass through. A two-level factor or
#' character vector is coded 0/1, and the reference level is reported, because
#' the sign of `effect_M` depends on it and a silent choice here is a silent
#' sign error downstream.
#'
#' @param x The exposure as resolved by `se_column()`.
#' @return A double vector.
#' @keywords internal
#' @noRd
se_exposure <- function(x) {
  if (is.numeric(x) || is.logical(x)) return(as.double(x))
  f <- factor(x)
  if (nlevels(f) != 2L)
    stop("exposure must be numeric, logical, or have exactly two levels; ",
         "got ", nlevels(f), " levels", call. = FALSE)
  message("exposure coded 1 for \"", levels(f)[2L], "\", 0 for \"",
          levels(f)[1L], "\"")
  as.double(f == levels(f)[2L])
}


#' Assay matrix on the M-value scale
#'
#' Picks the assay to analyse and guarantees it is on the M-value scale, which
#' is what every estimator in the package assumes.
#'
#' @details
#' `minfi::GenomicRatioSet()` stores `Beta` before `M` when both are present,
#' so taking the first assay silently analyses beta values. When `assay` is
#' left as `NULL` an assay named `M` is used if the object has one; otherwise,
#' if the object has exactly one assay, that one is used. Values that all lie
#' in `[0, 1]` are beta values and are converted with
#' `log2(beta / (1 - beta))`, clipped away from 0 and 1 first so that a probe
#' at either boundary gives a finite M-value rather than an infinity. The
#' conversion is reported with a message, never done silently.
#'
#' @param object A `SummarizedExperiment`.
#' @param assay Assay index or name, or `NULL` to choose as described above.
#' @return A numeric matrix on the M-value scale.
#' @keywords internal
#' @noRd
se_assay <- function(object, assay = NULL) {
  nms <- SummarizedExperiment::assayNames(object)
  n_assay <- length(SummarizedExperiment::assays(object))
  if (is.null(assay)) {
    if (!is.null(nms) && "M" %in% nms) {
      assay <- "M"
    } else if (n_assay == 1L) {
      assay <- 1L
    } else {
      stop("the object has ", n_assay, " assays (",
           paste(nms, collapse = ", "), ") and none is named \"M\"; name the ",
           "M-value assay \"M\" or pass assay = ", call. = FALSE)
    }
  }
  M <- SummarizedExperiment::assay(object, assay)
  if (!is.matrix(M)) M <- as.matrix(M)
  v <- M[is.finite(M)]
  if (length(v) && min(v) >= 0 && max(v) <= 1) {
    eps <- 1e-6
    message("assay \"", if (is.character(assay)) assay else nms[assay],
            "\" holds beta values; converting to M-values with ",
            "log2(beta / (1 - beta))")
    M[] <- log2(pmin(pmax(M, eps), 1 - eps) / (1 - pmin(pmax(M, eps), 1 - eps)))
  }
  M
}


#' Design columns resolved from colData
#'
#' @param object A `SummarizedExperiment`.
#' @param exposure,subject,covars As in [ewas_regions()].
#' @param n_samples Number of assay columns, for the length checks.
#' @return A list with `exposure` (double), `subject` (character) and
#' `covars` (matrix or `NULL`).
#' @keywords internal
#' @noRd
se_design <- function(object, exposure, subject, covars, n_samples) {
  cd <- as.data.frame(SummarizedExperiment::colData(object))
  if (nrow(cd) != n_samples)
    stop("colData has ", nrow(cd), " rows but the assay has ", n_samples,
         " columns", call. = FALSE)
  cmat <- NULL
  if (!is.null(covars) && length(covars)) {
    if (is.character(covars)) {
      miss <- setdiff(covars, names(cd))
      if (length(miss))
        stop("covars not in colData: ", paste(miss, collapse = ", "),
             call. = FALSE)
      sub <- cd[, covars, drop = FALSE]
      if (anyNA(sub))
        stop("covars contain missing values; model.matrix() would drop those ",
             "samples and the design would no longer match the assay",
             call. = FALSE)
      cmat <- stats::model.matrix(stats::reformulate(covars), data = sub)
      cmat <- cmat[, colnames(cmat) != "(Intercept)", drop = FALSE]
    } else {
      cmat <- as.matrix(covars)
      if (nrow(cmat) != n_samples)
        stop("covars has ", nrow(cmat), " rows but the assay has ", n_samples,
             " columns", call. = FALSE)
    }
  }
  list(exposure = se_exposure(se_column(exposure, cd, "exposure")),
       subject = as.character(se_column(subject, cd, "subject")),
       covars = cmat)
}


#' Fit the within-subject exposure effect on a SummarizedExperiment
#'
#' Extracts an assay and its design columns from a
#' [SummarizedExperiment::SummarizedExperiment] (or anything extending it,
#' including minfi's `GenomicRatioSet`) and passes them to [fit_within()].
#' The fit itself is the same code the pipeline runs; this function only
#' resolves the object.
#'
#' @details
#' `exposure`, `subject` and `covars` may name `colData` columns or be given
#' as vectors and a matrix. Named covariates are expanded with
#' [stats::model.matrix()] and the intercept is dropped, since the
#' within-subject transform removes it; a covariate with missing values is an
#' error rather than a silently shortened design.
#'
#' The estimators assume M-values. With `assay = NULL` the assay named `M` is
#' used when the object has one, which is what makes a `GenomicRatioSet`
#' behave as expected: `minfi` stores `Beta` first, so taking the first assay
#' would analyse beta values. A beta-valued assay is converted to M-values
#' with a message.
#'
#' @param object A `SummarizedExperiment`, or any object with `assay()` and
#' `colData()` methods.
#' @param exposure Character scalar naming a `colData` column, or a vector of
#' length `ncol(object)`. Must vary within subject.
#' @param subject Character scalar naming a `colData` column, or a vector of
#' length `ncol(object)`, giving the subject each sample belongs to.
#' @param covars Optional character vector naming `colData` columns, or a
#' matrix with `ncol(object)` rows.
#' @param assay Assay index or name. `NULL`, the default, takes the assay
#' named `M`, or the only assay when the object has one.
#' @param shrink_var,var_method Passed to [fit_within()].
#' @return The list returned by [fit_within()]: `beta`, `se`, `z`, `df`,
#' `n_subjects`, `n_samples`. `beta`, `se` and `z` carry the assay's row
#' names when it has them.
#' @seealso [fit_within()] for the matrix interface and for what the
#' within-subject transform does to the design.
#' @examples
#' set.seed(1)
#' M <- matrix(rnorm(40), nrow = 10,
#'             dimnames = list(paste0("cg", 1:10), NULL))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(M = M),
#'   colData = data.frame(dose = c(0, 1, 0, 1),
#'                        subject = c("s1", "s1", "s2", "s2")))
#' fit <- fit_within_se(se, "dose", "subject", shrink_var = FALSE)
#' head(fit$z)
#' @export
fit_within_se <- function(object, exposure, subject, covars = NULL,
                          assay = NULL, shrink_var = TRUE,
                          var_method = "limma") {
  M <- se_assay(object, assay)
  des <- se_design(object, exposure, subject, covars, ncol(M))
  fit_within(M, des$exposure, des$subject, covars = des$covars,
             shrink_var = shrink_var, var_method = var_method)
}


#' Probe chromosome and position from an annotated object
#'
#' Reads probe coordinates from `rowRanges()` when the object has them and
#' from `rowData()` otherwise, and returns them in the form
#' [comethylation_clusters()] and [call_regions()] take.
#'
#' @details
#' Coordinates come back in the object's own row order, unsorted: the region
#' caller treats adjacent rows as adjacent probes, so ordering is the caller's
#' decision and [sort_probes()] is how to make it. `rowData()` is searched for
#' a chromosome column named `chr`, `chrom` or `seqnames` and a position
#' column named `pos`, `position` or `start`, case-insensitively.
#'
#' @param object A `SummarizedExperiment` or `RangedSummarizedExperiment`.
#' @return A list with `chrom` (character) and `pos` (integer).
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(M = matrix(0, nrow = 3, ncol = 2)),
#'   rowData = data.frame(chr = "chr1", pos = c(10L, 20L, 30L)))
#' probe_coords(se)
#' @export
probe_coords <- function(object) {
  if (methods::is(object, "RangedSummarizedExperiment")) {
    # as.data.frame() rather than seqnames()/start(): those generics live in
    # GenomeInfoDb and BiocGenerics, and this keeps the layer's dependencies
    # to the two packages it already declares.
    d <- as.data.frame(SummarizedExperiment::rowRanges(object))
    return(list(chrom = as.character(d$seqnames), pos = as.integer(d$start)))
  }
  rd <- as.data.frame(SummarizedExperiment::rowData(object))
  find <- function(cands, what) {
    hit <- match(cands, tolower(names(rd)))
    hit <- hit[!is.na(hit)]
    if (!length(hit))
      stop("no ", what, " column in rowData(); looked for ",
           paste(cands, collapse = ", "), "; found ",
           paste(names(rd), collapse = ", "), call. = FALSE)
    rd[[hit[1L]]]
  }
  list(chrom = as.character(find(c("chr", "chrom", "seqnames"), "chromosome")),
       pos = as.integer(find(c("pos", "position", "start"), "position")))
}


#' Order an object by probe position
#'
#' @details
#' [call_regions()] and [call_blocks()] read runs of adjacent rows as runs of
#' adjacent probes, so an object whose rows are in array order — which is not
#' genomic order — produces regions that are not contiguous on the genome.
#' This applies the ordering once, before the fit, so every downstream stage
#' sees the same one.
#'
#' @param object A `SummarizedExperiment` or `RangedSummarizedExperiment`.
#' @return `object` with rows ordered by chromosome then position. The
#' chromosome order is [order()]'s on the chromosome names as characters,
#' which is stable but not karyotypic.
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(M = matrix(0, nrow = 3, ncol = 2)),
#'   rowData = data.frame(chr = "chr1", pos = c(30L, 10L, 20L)))
#' probe_coords(sort_probes(se))$pos
#' @export
sort_probes <- function(object) {
  co <- probe_coords(object)
  object[order(co$chrom, co$pos), ]
}


#' Region or block calls as a GRanges
#'
#' Converts the table returned by [call_regions()] or [call_blocks()] into a
#' [GenomicRanges::GRanges], with every remaining column as metadata.
#'
#' @details
#' Given the whole list from [call_regions()], the elements beside `table` —
#' `max_abs_z` — are kept in the object's `metadata()` rather than dropped,
#' since that value is the one a permutation null is compared against. The
#' table's `width` column is dropped: `GRanges` computes it, and `width` is
#' not a permitted metadata name.
#'
#' Coordinates pass through unchanged. Both callers report the first and last
#' probe or cluster position, 1-based and inclusive, which is what `GRanges`
#' means by `start` and `end`.
#'
#' @param x The list or table from [call_regions()], or the table from
#' [call_blocks()].
#' @return A `GRanges` with one range per called region or block.
#' @examples
#' tab <- data.frame(chr = "chr1", start = 100L, end = 400L, width = 301L,
#'                   n_probes = 4L, effect_M = 0.3, se = 0.02, z = 15,
#'                   p_pointwise = 1e-50, theta = 0.3,
#'                   probe_start = 0L, probe_end = 3L)
#' as_granges(list(table = tab, max_abs_z = 15))
#' @export
as_granges <- function(x) {
  meta <- list()
  if (!is.data.frame(x) && is.list(x)) {
    if (!"table" %in% names(x))
      stop("x is a list without a \"table\" element; pass the result of ",
           "call_regions() or call_blocks()", call. = FALSE)
    meta <- x[setdiff(names(x), "table")]
    x <- x$table
  }
  miss <- setdiff(c("chr", "start", "end"), names(x))
  if (length(miss))
    stop("table is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  x <- x[, setdiff(names(x), "width"), drop = FALSE]
  gr <- GenomicRanges::makeGRangesFromDataFrame(
    x, keep.extra.columns = TRUE, seqnames.field = "chr",
    start.field = "start", end.field = "end", ignore.strand = TRUE)
  if (length(meta)) S4Vectors::metadata(gr) <- meta
  gr
}


# ---------------------------------------------------------------------------
# High-level analysis: what the pipeline's stage drivers do, as two functions
# ---------------------------------------------------------------------------

#' Within-subject residuals for clustering
#'
#' Removes subject means from both the data and the design, then projects the
#' design out of the residuals, which is the quantity
#' [comethylation_clusters()] correlates. The same three lines are what stages
#' 04 and 05 of the pipeline run.
#'
#' @param M Assay matrix, probes by samples.
#' @param exposure,subject,covars Design, as resolved by `se_design()`.
#' @return A residual matrix with the dimensions of `M`.
#' @keywords internal
#' @noRd
within_residuals <- function(M, exposure, subject, covars) {
  resid <- demean_rows_by_group(M, subject)
  Xw <- demean_by_group(cbind(exposure, covars), subject)
  Xw <- Xw[, apply(abs(Xw), 2, max) > 1e-9, drop = FALSE]
  if (!ncol(Xw)) return(resid)
  resid - (resid %*% Xw %*% pinv(crossprod(Xw))) %*% t(Xw)
}


#' Chromosome-wise total-variation track
#'
#' @param chrom,pos Probe coordinates.
#' @param beta,se Probe-level effect and standard error.
#' @param lam0,decay_bp Penalty level and decay, passed to
#'   [distance_penalty()].
#' @param n_iter ADMM iterations, passed to [tv_denoise()].
#' @return A list with the denoised track `theta` and the precision weights
#'   `w`.
#' @keywords internal
#' @noRd
tv_track <- function(chrom, pos, beta, se, lam0, decay_bp, n_iter) {
  w <- ifelse(se > 0, 1 / se^2, 0)
  theta <- numeric(length(beta))
  for (cc in unique(chrom)) {
    idx <- which(chrom == cc)
    lam <- distance_penalty(as.numeric(pos[idx]), lam0 = lam0,
                            decay_bp = decay_bp, w = w[idx])
    theta[idx] <- tv_denoise(beta[idx], w[idx], lam, n_iter = n_iter)
  }
  list(theta = theta, w = w)
}


#' Choose the smoothness by leaving out whole subjects
#'
#' Held-out subjects rather than held-out samples: a subject's visits are the
#' unit the within-subject effect is estimated from, so splitting them across
#' train and test would leak the quantity being validated. The score is the
#' precision-weighted squared error of the training track against the held-out
#' subjects' own within-subject effects.
#'
#' @param M,chrom,pos,exposure,subject,covars Data and design.
#' @param grid Candidate `lam0` values.
#' @param cv_folds Number of subject folds.
#' @param decay_bp,n_iter,var_method,seed As in [ewas_regions()].
#' @return A list with `best` and the per-value table `cv`.
#' @keywords internal
#' @noRd
cv_lambda <- function(M, chrom, pos, exposure, subject, covars, grid,
                      cv_folds, decay_bp, n_iter, var_method, seed) {
  folds <- fold_assign(subject, cv_folds, seed = seed)
  rows <- list()
  skipped <- 0L
  for (lam0 in grid) {
    tot <- 0
    ntot <- 0L
    for (f in folds) {
      te <- subject %in% f
      tr <- !te
      if (length(unique(subject[te])) < 2 || length(unique(subject[tr])) < 4) {
        skipped <- skipped + 1L
        next
      }
      ctr <- if (!is.null(covars)) covars[tr, , drop = FALSE] else NULL
      cte <- if (!is.null(covars)) covars[te, , drop = FALSE] else NULL
      if (within_df(exposure[te], subject[te], cte) <= 0 ||
          within_df(exposure[tr], subject[tr], ctr) <= 0) {
        skipped <- skipped + 1L
        next
      }
      ftr <- fit_within(M[, tr, drop = FALSE], exposure[tr], subject[tr], ctr,
                        var_method = var_method)
      fte <- fit_within(M[, te, drop = FALSE], exposure[te], subject[te], cte,
                        var_method = var_method)
      th <- tv_track(chrom, pos, ftr$beta, ftr$se, lam0, decay_bp, n_iter)$theta
      wte <- ifelse(fte$se > 0, 1 / fte$se^2, 0)
      tot <- tot + sum(wte * (fte$beta - th)^2)
      ntot <- ntot + sum(wte > 0)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      lam0 = lam0, wsse = tot, n = ntot,
      score = if (ntot > 0) tot / ntot else NA_real_)
  }
  cv <- do.call(rbind, rows)
  if (!any(is.finite(cv$score)))
    stop("no cross-validation fold was estimable: with ",
         length(unique(subject)), " subjects and ", cv_folds, " folds every ",
         "split left too little within-subject variation. Use fewer folds, ",
         "drop single-visit subjects, or pass a single lambda.", call. = FALSE)
  list(best = cv$lam0[which.min(cv$score)], cv = cv)
}


#' Collapse co-methylated probes to cluster profiles
#'
#' The `cpgCollapse()` step of `minfi`'s block finder, but over clusters
#' defined by within-subject co-methylation rather than by fixed widths.
#'
#' @param M Assay matrix.
#' @param chrom,pos Probe coordinates.
#' @param cl Cluster index per probe, as returned by
#'   [comethylation_clusters()]; negative means unclustered.
#' @return A list with the cluster-mean matrix `M`, coordinates `chrom`,
#'   `start`, `end`, `mid`, and probe counts `n_probes`.
#' @keywords internal
#' @noRd
collapse_clusters <- function(M, chrom, pos, cl) {
  n_cl <- if (max(cl) >= 0) max(cl) + 1L else 0L
  if (!n_cl)
    stop("no co-methylation cluster was formed; raise max_gap or lower ",
         "rho_min", call. = FALSE)
  keep <- which(cl >= 0)
  cid <- cl[keep]
  o2 <- order(cid)
  keep <- keep[o2]
  cid <- cid[o2]
  runs <- split(keep, cid)
  Mc <- matrix(0, n_cl, ncol(M))
  cchr <- character(n_cl)
  cstart <- numeric(n_cl)
  cend <- numeric(n_cl)
  cnp <- integer(n_cl)
  for (k in seq_len(n_cl)) {
    rows <- runs[[as.character(k - 1L)]]
    Mc[k, ] <- colMeans(M[rows, , drop = FALSE])
    cchr[k] <- chrom[rows[1L]]
    cstart[k] <- min(pos[rows])
    cend[k] <- max(pos[rows])
    cnp[k] <- length(rows)
  }
  cmid <- (cstart + cend) / 2
  o <- order(cchr, cmid)
  list(M = Mc[o, , drop = FALSE], chrom = cchr[o], start = cstart[o],
       end = cend[o], mid = cmid[o], n_probes = cnp[o])
}


#' Permutation distribution of the genome-wide maximum |z|
#'
#' Does not call `set.seed()`: the replicates are drawn from the caller's RNG
#' stream, so reproducibility is the caller's `set.seed()` and this function
#' never resets a stream it does not own.
#'
#' @param M,co,expo,subj,cov Data, coordinates and design.
#' @param lam0,decay_bp,perm_iter,min_probes,min_effect,var_method Run
#'   settings, as in [ewas_regions()].
#' @param n_perm,perm_scheme Replicate count and permutation scheme.
#' @return A numeric vector of `n_perm` null maxima.
#' @keywords internal
#' @noRd
region_null <- function(M, co, expo, subj, cov, lam0, decay_bp, perm_iter,
                        min_probes, min_effect, var_method, n_perm,
                        perm_scheme) {
  fn <- if (perm_scheme == "within") within_subject_permutation else
    naive_permutation
  vapply(seq_len(n_perm), function(b) {
    pe <- fn(expo, subj)
    f <- fit_within(M, pe, subj, cov, var_method = var_method)
    th <- tv_track(co$chrom, co$pos, f$beta, f$se, lam0, decay_bp,
                   perm_iter)$theta
    call_regions(co$chrom, co$pos, f$beta, f$se, th, min_probes = min_probes,
                 min_effect = min_effect)$max_abs_z
  }, numeric(1))
}


#' Family-wise p-values and their Monte Carlo error
#'
#' @param reg Region table.
#' @param null_max Null maxima from [region_null()].
#' @param n_perm Replicate count.
#' @param say Progress reporter.
#' @return `reg` with `p_fwer` and `p_fwer_mcse` added.
#' @keywords internal
#' @noRd
add_fwer <- function(reg, null_max, n_perm, say) {
  reg$p_fwer <- vapply(reg$z, function(z)
    (1 + sum(null_max >= abs(z))) / (1 + length(null_max)), numeric(1))
  reg$p_fwer_mcse <- sqrt(reg$p_fwer * (1 - reg$p_fwer) / n_perm)
  borderline <- sum(reg$p_fwer - 2 * reg$p_fwer_mcse < 0.05 &
                    reg$p_fwer + 2 * reg$p_fwer_mcse > 0.05)
  if (borderline)
    say("  ", borderline, " region(s) have a p-value within two Monte Carlo ",
        "standard errors of 0.05; raise n_perm to resolve them")
  reg
}


#' Selection stability over subject half-samples
#'
#' A region counts as re-discovered when an interval called in the subsample
#' overlaps it, because boundaries move between subsamples. Like
#' [region_null()], this draws from the caller's RNG stream.
#'
#' @inheritParams region_null
#' @param reg Region table.
#' @param n_boot Number of half-samples.
#' @return `reg` with a `stability` column added.
#' @keywords internal
#' @noRd
add_stability <- function(reg, M, co, expo, subj, cov, lam0, decay_bp,
                          perm_iter, min_probes, min_effect, var_method,
                          n_boot) {
  keys <- sprintf("%s:%d-%d", reg$chr, reg$start, reg$end)
  fit_subset <- function(sub) {
    sel <- subj %in% sub
    cs <- if (!is.null(cov)) cov[sel, , drop = FALSE] else NULL
    if (length(unique(subj[sel])) < 3 ||
        within_df(expo[sel], subj[sel], cs) <= 0)
      return(character(0))
    f <- fit_within(M[, sel, drop = FALSE], expo[sel], subj[sel], cs,
                    var_method = var_method)
    th <- tv_track(co$chrom, co$pos, f$beta, f$se, lam0, decay_bp,
                   perm_iter)$theta
    r <- call_regions(co$chrom, co$pos, f$beta, f$se, th,
                      min_probes = min_probes, min_effect = min_effect)
    if (!nrow(r$table)) return(character(0))
    hit <- vapply(seq_len(nrow(reg)), function(i)
      any(r$table$chr == reg$chr[i] & r$table$start <= reg$end[i] &
          reg$start[i] <= r$table$end), logical(1))
    keys[hit]
  }
  ss <- stability_selection(fit_subset, subj, n_boot = n_boot, frac = 0.5)
  reg$stability <- ifelse(keys %in% names(ss$freq), ss$freq[keys], 0)
  reg$stability[is.na(reg$stability)] <- 0
  reg
}


#' Call differentially methylated regions from a methylation object

#'
#' The whole region analysis in one call: within-subject probe effects,
#' co-methylation clusters, smoothness chosen by held-out subjects,
#' total-variation region calling, and, when asked for, a permutation null and
#' subject-level selection stability. This is what `04_dmr_ml.R` does in the
#' pipeline.
#'
#' @details
#' The input is a [SummarizedExperiment::SummarizedExperiment] or anything
#' extending it, which includes `minfi`'s `GenomicRatioSet`. Rows are ordered
#' by chromosome and position internally, so an object in array order is
#' handled correctly; the returned ranges are in genomic order.
#'
#' Effects are identified within subject: subject means are differenced out of
#' both the data and the design, so subject identity, array generation and
#' chip are absorbed rather than modelled, and only covariates that vary
#' between a subject's visits contribute. The exposure must vary within
#' subject; a subject with one visit contributes nothing, which is reported
#' rather than silently dropped.
#'
#' `lambda` of length one is used as the penalty. A vector, the default, is
#' cross-validated over `cv_folds` folds of whole subjects, and the selected
#' value and the score table are returned in the object's `metadata()`.
#'
#' With `n_perm > 0` the exposure is permuted `n_perm` times, the whole region
#' call is repeated, and each region's family-wise p-value is the fraction of
#' null replicates whose maximum |z| reaches the region's own, with the
#' Monte Carlo standard error of that fraction beside it. The default scheme
#' permutes the exposure among the visits of the same subject, which leaves
#' subject, array, chip and the within-subject covariance untouched;
#' `"naive"` shuffles across subjects, the scheme `bumphunter` uses, and is
#' provided to quantify its mis-calibration rather than to be used.
#'
#' With `n_boot > 0`, regions are re-called on repeated half-samples of
#' subjects and `stability` reports the fraction of subsamples in which an
#' overlapping region was called again.
#'
#' @param object A `SummarizedExperiment`, `RangedSummarizedExperiment` or
#'   `GenomicRatioSet` with probe coordinates in `rowRanges()` or `rowData()`.
#' @param exposure Character scalar naming a `colData` column, or a vector of
#'   length `ncol(object)`. Must vary within subject.
#' @param subject Character scalar naming a `colData` column, or a vector of
#'   length `ncol(object)`.
#' @param covars Optional character vector naming `colData` columns, or a
#'   matrix with `ncol(object)` rows.
#' @param assay Assay index or name; `NULL` takes the assay named `M`, or the
#'   only assay. Beta values are converted with a message.
#' @param max_gap,rho_min Clustering: maximum distance in base pairs between
#'   adjacent probes, and the minimum within-subject co-methylation needed to
#'   join them. Passed to [comethylation_clusters()].
#' @param lambda Total-variation penalty. A single value is used as given; a
#'   vector is cross-validated over held-out subjects.
#' @param cv_folds Number of subject folds used when `lambda` is a vector.
#' @param decay_bp Distance scale over which the penalty decays, passed to
#'   [distance_penalty()].
#' @param min_probes,min_effect Region floors, passed to [call_regions()].
#' @param n_perm Permutation replicates for the family-wise p-value. `0`, the
#'   default, skips the null.
#' @param perm_scheme `"within"` permutes the exposure within subject;
#'   `"naive"` permutes across all samples.
#' @param n_boot Subject half-sample replicates for selection stability. `0`,
#'   the default, skips it.
#' @param n_iter,perm_iter ADMM iterations for the observed track and for the
#'   permutation and stability replicates.
#' @param var_method Variance moderation, `"limma"` or `"mom"`, passed to
#'   [fit_within()].
#' @param seed Seed for the deterministic subject-fold split used when
#'   `lambda` is cross-validated. The permutation and subsampling replicates
#'   are drawn from the session's RNG stream, so call `set.seed()` before
#'   `ewas_regions()` to make those reproducible; the function does not reset a
#'   stream it does not own.
#' @param verbose Report each stage as it runs.
#' @return A [GenomicRanges::GRanges] of called regions, in genomic order,
#'   with per-region metadata columns `n_probes`, `effect_M`, `se`, `z`,
#'   `p_pointwise`, `theta`, and `p_fwer`, `p_fwer_mcse` and `stability` when
#'   those were requested. The object's `metadata()` carries `max_abs_z`,
#'   `lambda`, `lambda_cv`, `n_clusters`, `fit` (degrees of freedom and
#'   subject and sample counts), `perm_scheme` and `null_max_abs_z`.
#' @seealso [ewas_blocks()] for long-range blocks, and [fit_within()],
#'   [comethylation_clusters()], [tv_denoise()] and [call_regions()] for the
#'   building blocks this assembles.
#' @examples
#' set.seed(1)
#' subject <- rep(sprintf("s%02d", 1:8), each = 3)
#' dose <- rep(c(0, 1, 2), times = 8)
#'
#' # probes in tight clusters separated by wide gaps, as on an array
#' pos <- as.integer(outer(c(0L, 150L, 300L), seq(1e5, by = 8000,
#'                                                length.out = 30L), "+"))
#' pos <- sort(pos)
#' n_probes <- length(pos)
#'
#' M <- matrix(rnorm(n_probes * length(subject), sd = 0.2), nrow = n_probes)
#' M <- M + rnorm(length(subject), sd = 0.5)[col(M)]
#' M <- M + rnorm(n_probes, sd = 1)[row(M)]
#' M[25:30, ] <- M[25:30, ] + 0.6 * rep(dose, each = 6)   # planted region
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(M = M),
#'   rowData = data.frame(chr = "chr1", pos = pos),
#'   colData = data.frame(dose = dose, subject = subject))
#' reg <- ewas_regions(se, "dose", "subject", lambda = 0.25, min_probes = 3L)
#' reg
#' @export
ewas_regions <- function(object, exposure, subject, covars = NULL,
                         assay = NULL, max_gap = 1000L, rho_min = 0.30,
                         lambda = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4),
                         cv_folds = 5L, decay_bp = 1000, min_probes = 3L,
                         min_effect = 0.05, n_perm = 0L,
                         perm_scheme = c("within", "naive"), n_boot = 0L,
                         n_iter = 3000L, perm_iter = 800L,
                         var_method = "limma", seed = 1L, verbose = FALSE) {
  perm_scheme <- match.arg(perm_scheme)
  object <- sort_probes(object)
  M <- se_assay(object, assay)
  co <- probe_coords(object)
  des <- se_design(object, exposure, subject, covars, ncol(M))
  expo <- des$exposure
  subj <- des$subject
  cov <- des$covars
  say <- function(...) if (verbose) message(...)

  n_vis <- table(subj)
  if (any(n_vis < 2))
    say(sum(n_vis < 2), " subject(s) have a single visit and contribute ",
        "nothing to a within-subject contrast")
  if (within_df(expo, subj, cov) <= 0)
    stop("no within-subject degrees of freedom: after differencing out ",
         "subject means the exposure and covariates leave nothing to ",
         "estimate from. The exposure must vary between visits of the same ",
         "subject.", call. = FALSE)

  say("probe-level within-subject fit")
  fit <- fit_within(M, expo, subj, cov, var_method = var_method)

  say("co-methylation clusters")
  cl <- comethylation_clusters(co$chrom, as.numeric(co$pos),
                               within_residuals(M, expo, subj, cov),
                               max_gap = max_gap, rho_min = rho_min)
  n_cl <- if (max(cl) >= 0) max(cl) + 1L else 0L

  cv <- NULL
  if (length(lambda) > 1L) {
    say("cross-validating the smoothness over held-out subjects")
    cvr <- cv_lambda(M, co$chrom, co$pos, expo, subj, cov, lambda, cv_folds,
                     decay_bp, perm_iter, var_method, seed)
    lam0 <- cvr$best
    cv <- cvr$cv
    say("  selected lambda = ", lam0)
  } else {
    lam0 <- as.numeric(lambda)
  }

  theta <- tv_track(co$chrom, co$pos, fit$beta, fit$se, lam0, decay_bp,
                    n_iter)$theta
  cl_arg <- ifelse(cl >= 0, cl, -seq_along(cl))
  obs <- call_regions(co$chrom, co$pos, fit$beta, fit$se, theta,
                      min_probes = min_probes, min_effect = min_effect,
                      cluster = cl_arg)
  reg <- obs$table
  say("regions called: ", nrow(reg))

  null_max <- numeric(0)
  if (nrow(reg) && n_perm > 0L) {
    say(perm_scheme, "-subject permutation null, ", n_perm, " replicates")
    null_max <- region_null(M, co, expo, subj, cov, lam0, decay_bp, perm_iter,
                            min_probes, min_effect, var_method, n_perm,
                            perm_scheme)
    reg <- add_fwer(reg, null_max, n_perm, say)
  }

  if (nrow(reg) && n_boot > 0L) {
    say("stability selection, ", n_boot, " subject half-samples")
    reg <- add_stability(reg, M, co, expo, subj, cov, lam0, decay_bp,
                         perm_iter, min_probes, min_effect, var_method, n_boot)
  }

  gr <- as_granges(list(
    table = reg, max_abs_z = obs$max_abs_z, lambda = lam0, lambda_cv = cv,
    n_clusters = n_cl,
    fit = list(df = fit$df, n_subjects = fit$n_subjects,
               n_samples = fit$n_samples),
    perm_scheme = if (length(null_max)) perm_scheme else NA_character_,
    null_max_abs_z = null_max))
  gr
}


#' Call long-range methylation blocks from a methylation object
#'
#' Collapses co-methylated probes into clusters, fits the within-subject
#' exposure effect per cluster, and calls blocks from a three-state hidden
#' semi-Markov model whose transition probabilities decay with genomic
#' distance. This is what `05_blocks_hsmm.R` does in the pipeline.
#'
#' @details
#' Blocks are a cluster-level statement, so the probe matrix is reduced first:
#' probes that are close and co-methylated within subject are averaged into
#' one profile, which is the `cpgCollapse()` idea with clusters defined by the
#' data rather than by fixed widths. The three states are hypomethylated,
#' neutral and hypermethylated with respect to the exposure; which column is
#' neutral is read from the fit rather than assumed, and when every cluster
#' effect falls on one side of zero only one direction can be called, which is
#' reported in `metadata()`.
#'
#' The model needs enough clusters to be identifiable. `min_fit_clusters`
#' refuses the fit below a floor rather than returning numbers from an
#' unidentified model; relaxing the clustering merges probes further and yields
#' fewer clusters, so the way out is to split clusters, not to relax them.
#'
#' @inheritParams ewas_regions
#' @param max_gap,rho_min Clustering, as in [ewas_regions()], with the looser
#'   defaults blocks are called at.
#' @param length_scale Distance scale of the state-transition decay, in base
#'   pairs. The 250 kb default matches the loess window it replaces.
#' @param min_post,min_clusters Block floors, passed to [call_blocks()].
#' @param min_fit_clusters Refuse to fit the block model on fewer clusters
#'   than this.
#' @return A [GenomicRanges::GRanges] of called blocks with metadata columns
#'   from [call_blocks()]. The object's `metadata()` carries `clusters`, a
#'   `GRanges` of the collapsed clusters with their effects and state
#'   posteriors, and `hsmm`, the fitted state means, standard deviations,
#'   stationary distribution, neutral column, log-likelihood and iteration
#'   count.
#' @seealso [ewas_regions()] for regions, and [fit_block_hsmm()] and
#'   [call_blocks()] for the building blocks.
#' @examples
#' set.seed(1)
#' subject <- rep(sprintf("s%02d", 1:10), each = 3)
#' dose <- rep(c(0, 1, 2), times = 10)
#'
#' # 60 clusters of three probes, so the block model has units to run along
#' pos <- sort(as.integer(outer(c(0L, 150L, 300L),
#'                              seq(1e5, by = 8000, length.out = 60L), "+")))
#' n_probes <- length(pos)
#'
#' M <- matrix(rnorm(n_probes * length(subject), sd = 0.2), nrow = n_probes)
#' M <- M + rnorm(length(subject), sd = 0.4)[col(M)]
#' M <- M + rnorm(n_probes, sd = 1)[row(M)]
#' shallow <- 60:120                                   # planted block
#' M[shallow, ] <- M[shallow, ] + 0.3 * rep(dose, each = length(shallow))
#'
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(M = M),
#'   rowData = data.frame(chr = "chr1", pos = pos),
#'   colData = data.frame(dose = dose, subject = subject))
#' blocks <- ewas_blocks(se, "dose", "subject", max_gap = 1000L,
#'                       rho_min = 0.2, min_fit_clusters = 20L)
#' blocks
#' @export
ewas_blocks <- function(object, exposure, subject, covars = NULL,
                        assay = NULL, max_gap = 1500L, rho_min = 0.20,
                        length_scale = 250000, min_post = 0.80,
                        min_clusters = 3L, min_fit_clusters = 200L,
                        var_method = "limma", seed = 1L, verbose = FALSE) {
  object <- sort_probes(object)
  M <- se_assay(object, assay)
  co <- probe_coords(object)
  des <- se_design(object, exposure, subject, covars, ncol(M))
  expo <- des$exposure
  subj <- des$subject
  cov <- des$covars
  say <- function(...) if (verbose) message(...)

  if (within_df(expo, subj, cov) <= 0)
    stop("no within-subject degrees of freedom: the exposure must vary ",
         "between visits of the same subject", call. = FALSE)

  say("co-methylation clusters")
  cl <- comethylation_clusters(co$chrom, as.numeric(co$pos),
                               within_residuals(M, expo, subj, cov),
                               max_gap = max_gap, rho_min = rho_min)
  n_cl <- if (max(cl) >= 0) max(cl) + 1L else 0L
  if (n_cl < min_fit_clusters)
    stop("too few co-methylation clusters (", n_cl, " < ", min_fit_clusters,
         ") to fit a block model; raise rho_min or lower max_gap to split ",
         "clusters, supply a larger panel, or lower min_fit_clusters if a ",
         "small test panel is intended", call. = FALSE)
  say("  ", n_cl, " clusters covering ", sum(cl >= 0), " probes")

  cc <- collapse_clusters(M, co$chrom, co$pos, cl)
  say("cluster-level within-subject fit")
  fit <- fit_within(cc$M, expo, subj, cov, var_method = var_method)

  say("distance-aware block model")
  hs <- fit_block_hsmm(fit$beta, fit$se, cc$mid, cc$chrom,
                       length_scale = length_scale, seed = seed)
  labs <- state_labels(hs$neutral)
  bl <- call_blocks(cc$chrom, cc$start, cc$end, hs$posterior,
                    min_post = min_post, min_clusters = min_clusters,
                    neutral = hs$neutral)
  say("blocks called: ", nrow(bl))

  clus <- data.frame(chr = cc$chrom, start = as.integer(cc$start),
                     end = as.integer(cc$end), n_probes = cc$n_probes,
                     effect_M = fit$beta, se = fit$se, z = fit$z,
                     stringsAsFactors = FALSE)
  for (k in seq_len(3L))
    clus[[paste0("post_", labs[k])]] <- hs$posterior[, k]

  as_granges(list(
    table = bl,
    clusters = as_granges(clus),
    hsmm = list(mu = hs$mu, sigma = hs$sigma, pi = hs$pi,
                neutral = hs$neutral, states = labs, loglik = hs$loglik,
                n_iter = hs$n_iter,
                one_sided = hs$neutral != 2L),
    fit = list(df = fit$df, n_subjects = fit$n_subjects,
               n_samples = fit$n_samples)))
}