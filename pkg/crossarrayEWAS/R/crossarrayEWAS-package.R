#' crossarrayEWAS: design-aware EWAS region and block detection
#'
#' Numerical core of the cross-array harmonisation pipeline for longitudinal
#' epigenome-wide association studies. The exported functions cover the two
#' stages the pipeline replaces: probe-level within-subject estimation with
#' region calling (\code{\link{fit_within}},
#' \code{\link{comethylation_clusters}}, \code{\link{distance_penalty}},
#' \code{\link{tv_denoise}}, \code{\link{call_regions}}), and long-range
#' block detection (\code{\link{fit_block_hsmm}},
#' \code{\link{call_blocks}}).
#'
#' @section Provenance:
#' These routines are an R port of the pipeline's Python reference module, and
#' parts of both were written with AI assistance (Assisted-by: Claude). The
#' port's deliberate divergences from the reference are recorded at each
#' function under \sQuote{Notes carried over from the pipeline source}.
#'
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges IRanges
#' @importFrom Matrix sparseMatrix bandSparse Cholesky solve
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay assayNames assays colData rowData
#'   rowRanges
#' @importFrom jsonlite fromJSON
#' @importFrom limma squeezeVar
#' @importFrom methods is
#' @importFrom stats median model.matrix pnorm quantile reformulate var
"_PACKAGE"
