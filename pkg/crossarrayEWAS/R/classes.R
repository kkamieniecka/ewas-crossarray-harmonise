# Class entry points: SummarizedExperiment in, GRanges out.
#
# This file is NOT generated from bin/ewasml.R, and it is the only part of the
# package that is not. The core is deliberately base R plus limma so that the
# pipeline's stage drivers can source it in an environment with no
# Bioconductor stack; putting assay extraction or GRanges construction in
# there would make that impossible. So the classes live here, the core stays
# matrix-in / data.frame-out, and tests/test_pkg_identity.R skips this file by
# name.
#
# SummarizedExperiment and GenomicRanges are Suggests rather than Imports for
# the same reason: the conda package (r-crossarrayewas) is installed beside
# the Galaxy wrappers, where the Bioconductor stack is not wanted. That choice
# also rules out S4 methods on those classes — setMethod() needs the generic's
# package at install time — which is why these are plain functions with
# requireNamespace() guards rather than methods on fit_within().
#
# Edit this file in pkg/R-src/; tools/build_pkg.py copies it into the
# generated tree, where it will be overwritten.


#' Column of `colData` as a model vector
#'
#' Resolves an argument that may either name a `colData` column or be the
#' vector itself, and fails with a message naming the available columns.
#'
#' @param x A character scalar naming a column, or a vector of length
#'   `nrow(cd)`.
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
#' The estimators assume M-values. An assay whose values all lie in `[0, 1]`
#' is almost certainly beta values, and gets a warning rather than a refusal,
#' because a legitimately M-value assay can land in that range when the probe
#' set is small.
#'
#' @param object A `SummarizedExperiment`, or any object with `assay()` and
#'   `colData()` methods.
#' @param exposure Character scalar naming a `colData` column, or a vector of
#'   length `ncol(object)`. Must vary within subject.
#' @param subject Character scalar naming a `colData` column, or a vector of
#'   length `ncol(object)`, giving the subject each sample belongs to.
#' @param covars Optional character vector naming `colData` columns, or a
#'   matrix with `ncol(object)` rows.
#' @param assay Assay index or name, passed to
#'   [SummarizedExperiment::assay()].
#' @param shrink_var,var_method Passed to [fit_within()].
#' @return The list returned by [fit_within()]: `beta`, `se`, `z`, `df`,
#'   `n_subjects`, `n_samples`. `beta`, `se` and `z` carry the assay's row
#'   names when it has them.
#' @seealso [fit_within()] for the matrix interface and for what the
#'   within-subject transform does to the design.
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   set.seed(1)
#'   M <- matrix(rnorm(40), nrow = 10,
#'               dimnames = list(paste0("cg", 1:10), NULL))
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(M = M),
#'     colData = data.frame(dose = c(0, 1, 0, 1),
#'                          subject = c("s1", "s1", "s2", "s2")))
#'   fit <- fit_within_se(se, "dose", "subject", shrink_var = FALSE)
#'   head(fit$z)
#' }
#' @export
fit_within_se <- function(object, exposure, subject, covars = NULL,
                          assay = 1L, shrink_var = TRUE,
                          var_method = "limma") {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("fit_within_se() needs the SummarizedExperiment package; install ",
         "it, or call fit_within() on a matrix instead.", call. = FALSE)
  M <- SummarizedExperiment::assay(object, assay)
  if (!is.matrix(M)) M <- as.matrix(M)
  cd <- as.data.frame(SummarizedExperiment::colData(object))
  if (nrow(cd) != ncol(M))
    stop("colData has ", nrow(cd), " rows but the assay has ", ncol(M),
         " columns", call. = FALSE)

  expo <- se_exposure(se_column(exposure, cd, "exposure"))
  subj <- as.character(se_column(subject, cd, "subject"))

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
      if (nrow(cmat) != ncol(M))
        stop("covars has ", nrow(cmat), " rows but the assay has ", ncol(M),
             " columns", call. = FALSE)
    }
  }

  # finite subset rather than range(na.rm = TRUE): an all-NA assay would warn,
  # and BiocCheck asks that the warning not simply be suppressed
  v <- M[is.finite(M)]
  if (length(v) && min(v) >= 0 && max(v) <= 1)
    warning("assay values all lie in [0, 1]: these look like beta values, ",
            "and the region and block callers assume M-values", call. = FALSE)

  fit_within(M, expo, subj, covars = cmat, shrink_var = shrink_var,
             var_method = var_method)
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
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(M = matrix(0, nrow = 3, ncol = 2)),
#'     rowData = data.frame(chr = "chr1", pos = c(10L, 20L, 30L)))
#'   probe_coords(se)
#' }
#' @export
probe_coords <- function(object) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("probe_coords() needs the SummarizedExperiment package",
         call. = FALSE)
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
#'   chromosome order is [order()]'s on the chromosome names as characters,
#'   which is stable but not karyotypic.
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(M = matrix(0, nrow = 3, ncol = 2)),
#'     rowData = data.frame(chr = "chr1", pos = c(30L, 10L, 20L)))
#'   probe_coords(sort_probes(se))$pos
#' }
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
#'   [call_blocks()].
#' @return A `GRanges` with one range per called region or block.
#' @examples
#' if (requireNamespace("GenomicRanges", quietly = TRUE)) {
#'   tab <- data.frame(chr = "chr1", start = 100L, end = 400L, width = 301L,
#'                     n_probes = 4L, effect_M = 0.3, se = 0.02, z = 15,
#'                     p_pointwise = 1e-50, theta = 0.3,
#'                     probe_start = 0L, probe_end = 3L)
#'   as_granges(list(table = tab, max_abs_z = 15))
#' }
#' @export
as_granges <- function(x) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("as_granges() needs the GenomicRanges package", call. = FALSE)
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