# Tests for the class entry points. Hand-written, not generated: the layer
# they cover is hand-written too (pkg/R-src/classes.R). Copied into the
# package tree by tools/build_pkg.py.
#
# SummarizedExperiment and GenomicRanges are Imports, so nothing here skips:
# if they are missing the package did not install.

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
  se <- se_fixture()
  direct <- fit_within(SummarizedExperiment::assay(se, "M"),
                       se$dose, se$subject)
  expect_identical(fit_within_se(se, "dose", "subject"), direct)
  # vectors instead of colData names must give the same answer
  expect_identical(fit_within_se(se, se$dose, se$subject), direct)
  expect_identical(names(direct$beta), rownames(se))
})

test_that("fit_within_se resolves and validates its design arguments", {
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
  se <- se_fixture()
  se$state <- ifelse(se$dose > 0, "on", "off")
  expect_message(fit <- fit_within_se(se, "state", "subject"),
                 "coded 1 for \"on\"")
  expect_equal(fit$beta, fit_within_se(se, "dose", "subject")$beta)
})

test_that("named covariates are expanded and the intercept dropped", {
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
  se <- se_fixture()
  se$age[3] <- NA
  expect_error(fit_within_se(se, "dose", "subject", covars = "age"),
               "missing values")
  expect_error(fit_within_se(se, "dose", "subject", covars = "height"),
               "not in colData")
})

test_that("an assay on the beta scale is converted, not analysed as is", {
  se <- se_fixture(beta_scale = TRUE)
  expect_message(fit_within_se(se, "dose", "subject"), "beta values")
  b <- SummarizedExperiment::assay(se, "M")
  direct <- fit_within(log2(b / (1 - b)), se$dose, se$subject)
  expect_equal(suppressMessages(fit_within_se(se, "dose", "subject")), direct,
               tolerance = 1e-8)
})

test_that("the M assay is preferred over Beta when both are present", {
  se <- se_fixture()
  M <- SummarizedExperiment::assay(se, "M")
  both <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Beta = 1 / (1 + 2^(-M)), M = M),
    rowData = as.data.frame(SummarizedExperiment::rowData(se)),
    colData = SummarizedExperiment::colData(se))
  # no conversion message: the M assay was taken, not the first one
  expect_silent(fit <- fit_within_se(both, "dose", "subject"))
  expect_identical(fit, fit_within(M, se$dose, se$subject))
  # and an ambiguous object says so rather than guessing
  two <- SummarizedExperiment::SummarizedExperiment(
    assays = list(first = M, second = M),
    colData = SummarizedExperiment::colData(se))
  expect_error(fit_within_se(two, "dose", "subject"), "none is named")
})
test_that("probe_coords reads rowRanges and rowData alike", {
  flat <- probe_coords(se_fixture())
  ranged <- probe_coords(se_fixture(ranged = TRUE))
  expect_identical(flat, ranged)
  expect_type(flat$chrom, "character")
  expect_type(flat$pos, "integer")
})

test_that("probe_coords accepts the other column spellings", {
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
  expect_error(as_granges(data.frame(chr = "chr1", start = 1L)), "end")
  expect_error(as_granges(list(z = 1)), "without a \"table\" element")
})

# --- the high-level functions ----------------------------------------------

ewas_fixture <- function(n_sub = 8L, n_clust = 60L, effect = 0.6) {
  set.seed(11)
  visits <- rep(3L, n_sub)
  subject <- rep(sprintf("p%02d", seq_len(n_sub)), visits)
  n <- length(subject)
  dose <- rep(c(0, 1, 2), times = n_sub)
  sizes <- rep(3L, n_clust)
  pos <- integer(0)
  at <- 100000L
  for (k in seq_len(n_clust)) {
    pos <- c(pos, at + c(0L, 150L, 300L))
    at <- at + 300L + 8000L
  }
  n_probes <- length(pos)
  M <- matrix(stats::rnorm(n_probes * n, sd = 0.2), nrow = n_probes)
  M <- M + stats::rnorm(n, sd = 0.5)[col(M)]            # subject-visit level
  M <- M + stats::rnorm(n_probes, sd = 1)[row(M)]       # probe baseline
  planted <- which(rep(seq_len(n_clust), sizes) %in% 20:21)
  M[planted, ] <- M[planted, ] + effect * rep(dose, each = length(planted))
  rr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(pos, width = 2L))
  names(rr) <- sprintf("cg%05d", seq_len(n_probes))
  rownames(M) <- names(rr)
  list(se = SummarizedExperiment::SummarizedExperiment(
         assays = list(M = M), rowRanges = rr,
         colData = data.frame(dose = dose, subject = subject,
                              age = stats::rnorm(n, 40, 5))),
       planted = range(pos[planted]))
}

test_that("ewas_regions recovers a planted region and reports its run", {
  fx <- ewas_fixture()
  gr <- ewas_regions(fx$se, "dose", "subject", lambda = 0.25, min_probes = 3L,
                     min_effect = 0.05)
  expect_s4_class(gr, "GRanges")
  expect_true(length(gr) >= 1L)
  hit <- GenomicRanges::start(gr) <= fx$planted[2] &
    fx$planted[1] <= GenomicRanges::end(gr)
  expect_true(any(hit))
  # the strongest region is the planted one, not something else
  expect_true(hit[which.max(abs(gr$z))])
  md <- S4Vectors::metadata(gr)
  expect_identical(md$lambda, 0.25)
  expect_true(md$n_clusters > 0L)
  expect_true(md$fit$df > 0L)
  expect_identical(md$fit$n_subjects, 8L)
})

test_that("ewas_regions assembles the same calls as the building blocks", {
  fx <- ewas_fixture()
  gr <- ewas_regions(fx$se, "dose", "subject", lambda = 0.25, min_probes = 3L,
                     min_effect = 0.05)
  se <- sort_probes(fx$se)
  co <- probe_coords(se)
  M <- SummarizedExperiment::assay(se, "M")
  fit <- fit_within(M, se$dose, se$subject)
  design <- stats::model.matrix(~ subject + dose,
                                data = as.data.frame(
                                  SummarizedExperiment::colData(se)))
  resid <- t(stats::lm.fit(design, t(M))$residuals)
  cl <- comethylation_clusters(co$chrom, as.numeric(co$pos), resid,
                               max_gap = 1000L, rho_min = 0.30)
  w <- ifelse(fit$se > 0, 1 / fit$se^2, 0)
  lam <- distance_penalty(as.numeric(co$pos), lam0 = 0.25, decay_bp = 1000,
                          w = w)
  theta <- tv_denoise(fit$beta, w, lam, n_iter = 3000L)
  obs <- call_regions(co$chrom, co$pos, fit$beta, fit$se, theta,
                      min_probes = 3L, min_effect = 0.05,
                      cluster = ifelse(cl >= 0, cl, -seq_along(cl)))
  expect_identical(as.integer(GenomicRanges::start(gr)), obs$table$start)
  expect_equal(gr$z, obs$table$z)
})

test_that("ewas_regions adds permutation and stability columns on request", {
  fx <- ewas_fixture()
  gr <- ewas_regions(fx$se, "dose", "subject", lambda = 0.25, n_perm = 5L,
                     n_boot = 4L, perm_iter = 100L)
  expect_true(all(c("p_fwer", "p_fwer_mcse", "stability") %in%
                    names(GenomicRanges::mcols(gr))))
  expect_true(all(gr$p_fwer > 0 & gr$p_fwer <= 1))
  expect_true(all(gr$stability >= 0 & gr$stability <= 1))
  md <- S4Vectors::metadata(gr)
  expect_length(md$null_max_abs_z, 5L)
  expect_identical(md$perm_scheme, "within")
  # replicates come from the session RNG, so the same seed set by the caller
  # gives the same p-values
  set.seed(99)
  a <- ewas_regions(fx$se, "dose", "subject", lambda = 0.25, n_perm = 5L,
                    perm_iter = 100L)
  set.seed(99)
  b <- ewas_regions(fx$se, "dose", "subject", lambda = 0.25, n_perm = 5L,
                    perm_iter = 100L)
  expect_equal(a$p_fwer, b$p_fwer)
})

test_that("ewas_regions cross-validates lambda when given a grid", {
  fx <- ewas_fixture()
  gr <- ewas_regions(fx$se, "dose", "subject", lambda = c(0.25, 1),
                     cv_folds = 4L, perm_iter = 100L)
  md <- S4Vectors::metadata(gr)
  expect_true(md$lambda %in% c(0.25, 1))
  expect_identical(nrow(md$lambda_cv), 2L)
  expect_true(all(is.finite(md$lambda_cv$score)))
})

test_that("ewas_regions refuses a design with no within-subject information", {
  fx <- ewas_fixture()
  se <- fx$se
  # exposure constant within subject: nothing left after differencing
  SummarizedExperiment::colData(se)$dose <-
    as.numeric(factor(se$subject))
  expect_error(ewas_regions(se, "dose", "subject", lambda = 0.25),
               "no within-subject degrees of freedom")
})

test_that("ewas_blocks returns blocks with clusters and states attached", {
  fx <- ewas_fixture(n_clust = 120L, effect = 0.3)
  gr <- ewas_blocks(fx$se, "dose", "subject", max_gap = 1000L,
                    rho_min = 0.20, min_fit_clusters = 50L)
  expect_s4_class(gr, "GRanges")
  md <- S4Vectors::metadata(gr)
  expect_s4_class(md$clusters, "GRanges")
  expect_true(length(md$clusters) >= 50L)
  expect_length(md$hsmm$mu, 3L)
  expect_true(md$hsmm$neutral %in% 1:3)
  # when every cluster effect lies on one side of zero the neutral state sorts
  # to an end and the labels say so, which is what one_sided reports
  expect_length(md$hsmm$states, 3L)
  expect_true("neutral" %in% md$hsmm$states)
  expect_identical(md$hsmm$one_sided, md$hsmm$neutral != 2L)
  expect_true(all(c("effect_M", "se", "z") %in%
                    names(GenomicRanges::mcols(md$clusters))))
})

test_that("ewas_blocks declines an unidentifiable block model", {
  fx <- ewas_fixture(n_clust = 20L)
  expect_error(ewas_blocks(fx$se, "dose", "subject", min_fit_clusters = 200L),
               "too few co-methylation clusters")
})
