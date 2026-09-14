block_fit <- function(n = 40) {
  set.seed(7)
  chrom <- rep("chr1", n)
  pos <- as.integer(seq(1e5, 5e6, length.out = n))
  y <- c(rnorm(n %/% 2, 0.6, 0.1), rnorm(n - n %/% 2, 0, 0.1))
  se <- rep(0.1, n)
  list(chrom = chrom, pos = pos, y = y, se = se,
       fit = fit_block_hsmm(y, se, pos, chrom, n_iter = 30L))
}

test_that("the fitted posterior is a proper distribution over three states", {
  d <- block_fit()
  post <- d$fit$posterior
  expect_identical(dim(post), c(length(d$y), 3L))
  expect_equal(unname(rowSums(post)), rep(1, length(d$y)), tolerance = 1e-8)
  expect_true(all(post >= -1e-12))
})

test_that("states are ordered by mu and the neutral state sits at zero", {
  d <- block_fit()
  expect_false(is.unsorted(d$fit$mu))
  expect_equal(d$fit$mu[d$fit$neutral], 0, tolerance = 1e-12)
  expect_true(d$fit$neutral %in% 1:3)
})

test_that("the fit is deterministic and seed is ignored", {
  d <- block_fit()
  again <- fit_block_hsmm(d$y, d$se, d$pos, d$chrom, n_iter = 30L, seed = 999L)
  expect_equal(again$posterior, d$fit$posterior, tolerance = 1e-12)
  expect_equal(again$loglik, d$fit$loglik, tolerance = 1e-12)
})

test_that("blocks respect min_clusters and min_post", {
  d <- block_fit()
  blocks <- call_blocks(d$chrom, d$pos, d$pos + 500L, d$fit$posterior,
                        neutral = d$fit$neutral, min_clusters = 3L)
  expect_true(all(c("chr", "start", "end", "width", "n_clusters",
                    "direction", "posterior") %in% names(blocks)))
  if (nrow(blocks)) {
    expect_true(all(blocks$n_clusters >= 3L))
    expect_true(all(blocks$posterior >= 0.80))
    expect_true(all(blocks$width == blocks$end - blocks$start + 1L))
  }
  expect_identical(nrow(call_blocks(d$chrom, d$pos, d$pos + 500L,
                                    d$fit$posterior,
                                    neutral = d$fit$neutral,
                                    min_clusters = 999L)), 0L)
})

test_that("mislabelling the neutral column changes block directions", {
  d <- block_fit()
  right <- call_blocks(d$chrom, d$pos, d$pos + 500L, d$fit$posterior,
                       neutral = d$fit$neutral, min_clusters = 3L)
  wrong_col <- if (d$fit$neutral == 1L) 3L else 1L
  wrong <- call_blocks(d$chrom, d$pos, d$pos + 500L, d$fit$posterior,
                       neutral = wrong_col, min_clusters = 3L)
  if (nrow(right) && nrow(wrong))
    expect_false(identical(sort(right$direction), sort(wrong$direction)))
})

test_that("state labels stay unambiguous when the pinned state is at an end", {
  expect_identical(state_labels(2), c("hypo", "neutral", "hyper"))
  lab <- state_labels(1)
  expect_length(lab, 3L)
  expect_identical(lab[1], "neutral")
  expect_identical(length(unique(lab)), 3L)
})
