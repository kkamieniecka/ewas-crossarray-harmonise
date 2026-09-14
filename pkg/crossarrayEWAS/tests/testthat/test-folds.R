test_that("fold assignment is deterministic and order-invariant", {
  subj <- sprintf("S%02d", 1:11)
  a <- fold_assign(subj, 3)
  expect_identical(a, fold_assign(subj, 3))
  b <- fold_assign(rev(subj), 3)
  expect_identical(lapply(a, sort), lapply(b, sort))
  expect_setequal(unlist(a), subj)
  expect_identical(length(unlist(a)), length(subj))   # each subject once
})

test_that("a different seed gives a different split", {
  subj <- sprintf("S%02d", 1:11)
  expect_false(identical(fold_assign(subj, 3, seed = 0),
                         fold_assign(subj, 3, seed = 7)))
})

test_that("visits of one subject never straddle folds", {
  subj <- rep(sprintf("S%02d", 1:8), each = 3)
  folds <- fold_assign(subj, 4)
  expect_setequal(unlist(folds), unique(subj))
  expect_identical(length(unlist(folds)), length(unique(subj)))
})

test_that("within_df reports non-identifiability as zero", {
  # exposure constant within subject: the within transform annihilates it
  expect_identical(within_df(c(0, 0, 1, 1), c("A", "A", "B", "B")), 0L)
  expect_gt(within_df(c(0, 1, 0, 1), c("A", "A", "B", "B")), 0L)
})

test_that("a time-invariant covariate is removed exactly", {
  ex <- rep(c(0, 1), 4)
  su <- rep(c("A", "B", "C", "D"), each = 2)
  age <- rep(c(40, 55, 61, 33), each = 2)
  expect_identical(within_df(ex, su, covars = cbind(age)),
                   within_df(ex, su))
})
