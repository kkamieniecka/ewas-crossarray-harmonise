# Tests for the class entry points. Hand-written, not generated: the layer
# they cover is hand-written too (pkg/R-src/classes.R). Copied into the
# package tree by tools/build_pkg.py.
#
# Every block skips when the suggested package is absent, because the layer is
# optional by design — the conda build installs the package without the
# Bioconductor stack.

se_fixture <- function(n_probes = 12L, ranged = FALSE, beta_scale = FALSE) {
  set.seed(7)
  n_samp <- 8L
  M <- matrix(stats::rnorm(n_probes * n_samp), nrow = n_probes,
              dimnames = list(paste0("cg", seq_len(n_probes)),
                              paste0("s", seq_len(n_samp))))
  if (beta_scale) M <- 1 / (1 + 2^(-M))
  cd <- data.frame(
    dose = c(0, 1, 0, 1, 0, 1, 0, 1),
    subject = rep(paste0("p", 1:4), each = 2L),
    arm = factor(rep(c("450k", "epic"), each = 4L)),
    age = c(40, 41, 55, 56, 61, 62, 33, 34))
  rd <- data.frame(chr = rep(c("chr1", "chr2"), each = n_probes / 2L),
                   pos = rep(seq.int(100L, by = 100L,
                                     length.out = n_probes / 2L), 2L))
  if (ranged) {
    rr <- GenomicRanges::GRanges(rd$chr,
            IRanges::IRanges(start = rd$pos, width = 2L))
    names(rr) <- rownames(M)
    SummarizedExperiment::SummarizedExperiment(
      assays = list(M = M), rowRanges = rr, colData = cd)
  } else {
    SummarizedExperiment::SummarizedExperiment(
      assays = list(M = M), rowData = rd, colData = cd)
  }
}

test_that("fit_within_se reproduces fit_within on the same data", {
  skip_if_not_installed("SummarizedExperiment")
  se <- se_fixture()
  direct <- fit_within(SummarizedExperiment::assay(se, "M"),
                       se$dose, se$subject)
  expect_identical(fit_within_se(se, "dose", "subject"), direct)
  # vectors instead of colData names must give the same answer
  expect_identical(fit_within_se(se, se$dose, se$subject), direct)
  expect_identical(names(direct$beta), rownames(se))
})

test_that("fit_within_se resolves and validates its design arguments", {
  skip_if_not_installed("SummarizedExperiment")
  se <- se_fixture()
  expect_error(fit_within_se(se, "treatment", "subject"),
               "not a colData column")
  expect_error(fit_within_se(se, "treatment", "subject"), "dose")
  expect_error(fit_within_se(se, c(0, 1), "subject"), "length 2")
  # arm is a factor, so the 0/1 coding message fires before the fit refuses
  expect_error(suppressMessages(fit_within_se(se, "arm", "subject")),
               "no within-subject variation")
})

test_that("a two-level exposure is coded against its first level", {
  skip_if_not_installed("SummarizedExperiment")
  se <- se_fixture()
  se$state <- ifelse(se$dose > 0, "on", "off")
  expect_message(fit <- fit_within_se(se, "state", "subject"),
                 "coded 1 for \"on\"")
  expect_equal(fit$beta, fit_within_se(se, "dose", "subject")$beta)
})

test_that("named covariates are expanded and the intercept dropped", {
  skip_if_not_installed("SummarizedExperiment")
  se <- se_fixture()
  by_name <- fit_within_se(se, "dose", "subject", covars = c("age", "arm"))
  X <- stats::model.matrix(~ age + arm,
                           data = as.data.frame(
                             SummarizedExperiment::colData(se)))
  by_matrix <- fit_within_se(se, "dose", "subject",
                             covars = X[, colnames(X) != "(Intercept)",
                                        drop = FALSE])
  expect_identical(by_name, by_matrix)
  # arm is constant within subject, so the within transform drops it and the
  # residual degrees of freedom only pay for age
  expect_identical(by_name$df, fit_within_se(se, "dose", "subject",
                                             covars = "age")$df)
})

test_that("missing covariate values are an error, not a shortened design", {
  skip_if_not_installed("SummarizedExperiment")
  se <- se_fixture()
  se$age[3] <- NA
  expect_error(fit_within_se(se, "dose", "subject", covars = "age"),
               "missing values")
  expect_error(fit_within_se(se, "dose", "subject", covars = "height"),
               "not in colData")
})

test_that("an assay on the beta scale is warned about", {
  skip_if_not_installed("SummarizedExperiment")
  expect_warning(fit_within_se(se_fixture(beta_scale = TRUE),
                               "dose", "subject"),
                 "look like beta values")
  expect_silent(fit_within_se(se_fixture(), "dose", "subject"))
})

test_that("probe_coords reads rowRanges and rowData alike", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("GenomicRanges")
  flat <- probe_coords(se_fixture())
  ranged <- probe_coords(se_fixture(ranged = TRUE))
  expect_identical(flat, ranged)
  expect_type(flat$chrom, "character")
  expect_type(flat$pos, "integer")
})

test_that("probe_coords accepts the other column spellings", {
  skip_if_not_installed("SummarizedExperiment")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(M = matrix(0, nrow = 2L, ncol = 2L)),
    rowData = data.frame(Seqnames = c("chr1", "chr1"),
                         Position = c(5L, 9L)))
  expect_identical(probe_coords(se),
                   list(chrom = c("chr1", "chr1"), pos = c(5L, 9L)))
  bare <- SummarizedExperiment::SummarizedExperiment(
    assays = list(M = matrix(0, nrow = 2L, ncol = 2L)),
    rowData = data.frame(gene = c("A", "B")))
  expect_error(probe_coords(bare), "no chromosome column")
})

test_that("sort_probes orders rows and carries the assay with them", {
  skip_if_not_installed("SummarizedExperiment")
  se <- se_fixture()
  shuffled <- se[c(4L, 1L, 7L, 2L, 3L, 5L, 6L, 8L:12L), ]
  sorted <- sort_probes(shuffled)
  co <- probe_coords(sorted)
  expect_false(is.unsorted(order(co$chrom, co$pos)))
  expect_identical(rownames(sorted), rownames(se))
  expect_identical(SummarizedExperiment::assay(sorted, "M"),
                   SummarizedExperiment::assay(se, "M"))
})

test_that("as_granges converts a region call and keeps max_abs_z", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("SummarizedExperiment")
  se <- sort_probes(se_fixture())
  co <- probe_coords(se)
  fit <- fit_within_se(se, "dose", "subject")
  theta <- rep(c(0.4, 0), times = c(4L, 8L))
  reg <- call_regions(co$chrom, co$pos, fit$beta, fit$se, theta,
                      min_probes = 3L, min_effect = 0)
  gr <- as_granges(reg)
  expect_s4_class(gr, "GRanges")
  expect_identical(length(gr), nrow(reg$table))
  expect_identical(as.integer(GenomicRanges::start(gr)), reg$table$start)
  expect_identical(GenomicRanges::width(gr), reg$table$width)
  expect_false("width" %in% names(GenomicRanges::mcols(gr)))
  expect_identical(GenomicRanges::mcols(gr)$n_probes, reg$table$n_probes)
  expect_identical(S4Vectors::metadata(gr)$max_abs_z, reg$max_abs_z)
  # the table on its own is accepted, and then there is no metadata
  expect_identical(S4Vectors::metadata(as_granges(reg$table)), list())
})

test_that("as_granges converts a block call and an empty table", {
  skip_if_not_installed("GenomicRanges")
  blocks <- data.frame(chr = c("chr1", "chr2"), start = c(100L, 50L),
                       end = c(900L, 400L), width = c(801L, 351L),
                       n_clusters = c(9L, 4L),
                       direction = c("hypo", "hyper"),
                       posterior = c(0.97, 0.88),
                       stringsAsFactors = FALSE)
  gr <- as_granges(blocks)
  expect_identical(GenomicRanges::mcols(gr)$direction, blocks$direction)
  expect_identical(as.character(GenomicRanges::seqnames(gr)), blocks$chr)
  empty <- as_granges(blocks[0, ])
  expect_s4_class(empty, "GRanges")
  expect_identical(length(empty), 0L)
})

test_that("as_granges rejects tables it cannot place", {
  skip_if_not_installed("GenomicRanges")
  expect_error(as_granges(data.frame(chr = "chr1", start = 1L)), "end")
  expect_error(as_granges(list(z = 1)), "without a \"table\" element")
})