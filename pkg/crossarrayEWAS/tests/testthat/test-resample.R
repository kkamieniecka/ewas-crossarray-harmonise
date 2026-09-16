test_that("within-subject permutation preserves each subject's exposures", {
  set.seed(1)
  ex <- c(0, 1, 2, 10, 11, 12)
  su <- c("A", "A", "A", "B", "B", "B")
  for (i in 1:20) {
    p <- within_subject_permutation(ex, su)
    expect_setequal(p[su == "A"], ex[su == "A"])
    expect_setequal(p[su == "B"], ex[su == "B"])
  }
})

test_that("a single-visit subject contributes nothing under the within null", {
  set.seed(1)
  ex <- c(0, 1, 99)
  su <- c("A", "A", "B")
  expect_true(all(replicate(20, within_subject_permutation(ex, su)[3]) == 99))
})

test_that("the naive permutation breaks subject pairing", {
  set.seed(1)
  ex <- c(0, 1, 2, 10, 11, 12)
  su <- c("A", "A", "A", "B", "B", "B")
  p <- replicate(50, naive_permutation(ex, su))
  expect_setequal(as.vector(p[, 1]), ex)
  # at least one draw moves a value across subjects
  expect_true(any(apply(p, 2, function(v) !setequal(v[1:3], ex[1:3]))))
})

test_that("stability selection returns frequencies in [0, 1]", {
  subjects <- sprintf("S%02d", 1:12)
  fit_fn <- function(sub) {
    out <- "always"
    if ("S01" %in% sub) out <- c(out, "sometimes")
    out
  }
  ss <- stability_selection(fit_fn, subjects, n_boot = 40L, frac = 0.5)
  expect_identical(ss$n_boot, 40L)
  expect_true(all(ss$freq >= 0 & ss$freq <= 1))
  expect_equal(unname(ss$freq[["always"]]), 1)
  expect_lt(ss$freq[["sometimes"]], 1)
})

test_that("stability selection is reproducible from the caller's seed", {
  subjects <- sprintf("S%02d", 1:12)
  fit_fn <- function(sub) sub[1]
  set.seed(3L)
  a <- stability_selection(fit_fn, subjects, 20L)$freq
  set.seed(3L)
  b <- stability_selection(fit_fn, subjects, 20L)$freq
  expect_identical(a, b)
})

test_that("stability selection leaves the caller's stream where it found it", {
  subjects <- sprintf("S%02d", 1:12)
  fit_fn <- function(sub) sub[1]
  set.seed(11L)
  before <- runif(1)
  set.seed(11L)
  invisible(stability_selection(fit_fn, subjects, 5L))
  set.seed(11L)
  expect_identical(runif(1), before)
})
