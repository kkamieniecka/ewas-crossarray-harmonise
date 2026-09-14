test_that("moderate_variance shrinks towards the prior", {
  set.seed(1)
  s2 <- rchisq(500, df = 6) / 6
  for (m in c("limma", "mom")) {
    out <- moderate_variance(s2, df = 6, method = m)
    expect_length(out, length(s2))
    expect_true(all(out > 0))
    expect_lt(var(out), var(s2))                     # shrinkage
    expect_lt(max(out), max(s2))
    expect_gt(min(out), min(s2))
  }
})

test_that("the two variance fits agree on the same model", {
  set.seed(2)
  s2 <- rchisq(500, df = 10) / 10
  a <- moderate_variance(s2, df = 10, method = "limma")
  b <- moderate_variance(s2, df = 10, method = "mom")
  expect_gt(cor(a, b), 0.99)
  expect_lt(abs(median(a) - median(b)) / median(a), 0.05)
})

test_that("too few probes leaves the variances untouched", {
  s2 <- rchisq(20, df = 4) / 4
  expect_identical(moderate_variance(s2, df = 4, method = "mom"), s2)
})
