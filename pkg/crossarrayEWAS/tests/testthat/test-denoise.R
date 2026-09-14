test_that("the penalty decays with distance and respects its floor", {
  pos <- as.integer(c(1000, 1100, 2000, 20000))
  lam <- distance_penalty(pos, lam0 = 1, decay_bp = 1000)
  expect_length(lam, length(pos) - 1L)
  expect_true(all(diff(lam) < 0))                    # monotone in distance
  expect_true(all(distance_penalty(pos, 1, lam_floor = 0.5) >= 0.5))
})

test_that("precision weights make lam0 dimensionless", {
  pos <- as.integer(seq(1000, 1900, by = 100))
  w <- rep(4, length(pos))
  expect_equal(distance_penalty(pos, 1, w = w),
               4 * distance_penalty(pos, 1), tolerance = 1e-12)
})

test_that("tv_denoise interpolates between the data and a constant", {
  y <- c(rep(1, 6), rep(-1, 6))
  w <- rep(1, 12)
  loose <- tv_denoise(y, w, lam = rep(1e-8, 11), n_iter = 2000L)
  tight <- tv_denoise(y, w, lam = rep(1e6, 11), n_iter = 2000L)
  expect_equal(loose, y, tolerance = 1e-3)
  expect_equal(tight, rep(mean(y), 12), tolerance = 1e-3)
})

test_that("the fit is piecewise constant and segments shrink as lam grows", {
  y <- c(rep(1, 5), rep(0.2, 5), rep(-0.8, 5))
  w <- rep(1, 15)
  n_seg <- function(lam)
    length(unique(round(tv_denoise(y, w, rep(lam, 14), n_iter = 2000L), 6)))
  expect_lte(n_seg(0.5), n_seg(0.01))
  expect_identical(n_seg(1e6), 1L)
})

test_that("each segment is the precision-weighted mean of its probes", {
  y <- c(rep(1, 5), rep(-1, 5))
  w <- c(rep(1, 5), rep(9, 5))
  theta <- tv_denoise(y, w, rep(0.05, 9), n_iter = 3000L)
  for (v in unique(round(theta, 6))) {
    k <- which(round(theta, 6) == v)
    expect_equal(v, sum(w[k] * y[k]) / sum(w[k]), tolerance = 1e-4)
  }
})
