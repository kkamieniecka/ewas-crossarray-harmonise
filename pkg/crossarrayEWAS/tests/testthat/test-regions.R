track <- function() {
  chrom <- rep("chr1", 12)
  pos <- as.integer(seq(1000, 3200, length.out = 12))
  beta <- c(rep(0.4, 6), rep(-0.1, 6))
  se <- c(rep(0.1, 6), rep(0.05, 6))
  lam <- distance_penalty(pos, lam0 = 0.5, w = 1 / se^2)
  theta <- tv_denoise(beta, 1 / se^2, lam, n_iter = 3000L)
  list(chrom = chrom, pos = pos, beta = beta, se = se, theta = theta)
}

test_that("region boundaries fall where the denoised track jumps", {
  d <- track()
  res <- call_regions(d$chrom, d$pos, d$beta, d$se, d$theta,
                      min_probes = 3L, min_effect = 0.01)
  expect_true(nrow(res$table) >= 1)
  expect_true(all(res$table$start %in% d$pos))
  expect_true(all(res$table$end %in% d$pos))
  expect_true(all(res$table$end >= res$table$start))
})

test_that("the region effect is the precision-weighted mean of its probes", {
  d <- track()
  res <- call_regions(d$chrom, d$pos, d$beta, d$se, d$theta,
                      min_probes = 3L, min_effect = 0.01)
  eff <- grep("^(beta|effect)", names(res$table), value = TRUE)[1]
  for (i in seq_len(nrow(res$table))) {
    k <- (res$table$probe_start[i] + 1L):(res$table$probe_end[i] + 1L)
    w <- 1 / d$se[k]^2
    expect_equal(res$table[[eff]][i], sum(w * d$beta[k]) / sum(w),
                 tolerance = 1e-6)
  }
})

test_that("min_probes and min_effect filter regions out", {
  d <- track()
  expect_identical(
    nrow(call_regions(d$chrom, d$pos, d$beta, d$se, d$theta,
                      min_probes = 99L)$table), 0L)
  expect_identical(
    nrow(call_regions(d$chrom, d$pos, d$beta, d$se, d$theta,
                      min_effect = 10)$table), 0L)
})

test_that("a cluster change forces a boundary", {
  d <- track()
  one <- call_regions(d$chrom, d$pos, d$beta, d$se, d$theta, min_probes = 3L,
                      min_effect = 0.01, cluster = rep(1L, 12))
  two <- call_regions(d$chrom, d$pos, d$beta, d$se, d$theta, min_probes = 3L,
                      min_effect = 0.01, cluster = rep(c(1L, 2L), each = 6))
  expect_gte(nrow(two$table), nrow(one$table))
})
