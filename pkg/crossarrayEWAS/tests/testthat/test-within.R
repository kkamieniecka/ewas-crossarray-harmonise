planted_fit <- function(effect = 1.5, n_probes = 30, n_sub = 8) {
  set.seed(42)
  subject <- rep(sprintf("S%02d", seq_len(n_sub)), each = 2)
  exposure <- rep(c(0, 1), times = n_sub)
  M <- matrix(rnorm(n_probes * length(subject)), nrow = n_probes)
  M <- M + rep(rnorm(n_sub, sd = 3), each = 2)[col(M)]   # subject intercepts
  M[1:5, exposure == 1] <- M[1:5, exposure == 1] + effect
  list(M = M, exposure = exposure, subject = subject)
}

test_that("fit_within recovers a planted effect and absorbs subject means", {
  d <- planted_fit()
  fit <- fit_within(d$M, d$exposure, d$subject)
  expect_length(fit$beta, nrow(d$M))
  expect_identical(fit$n_subjects, 8L)
  expect_identical(fit$n_samples, ncol(d$M))
  expect_gt(mean(fit$beta[1:5]), 1)                  # planted probes
  expect_lt(abs(mean(fit$beta[6:30])), 0.5)          # null probes
  expect_true(all(abs(fit$z[1:5]) > abs(median(fit$z[6:30]))))
})

test_that("large subject intercepts do not move the estimate", {
  d <- planted_fit()
  shift <- rep(rnorm(8, sd = 50), each = 2)
  fit_a <- fit_within(d$M, d$exposure, d$subject)
  fit_b <- fit_within(d$M + shift[col(d$M)], d$exposure, d$subject)
  expect_equal(fit_a$beta, fit_b$beta, tolerance = 1e-8)
})

test_that("variance moderation only changes the standard errors", {
  d <- planted_fit(n_probes = 300)   # moderation no-ops below 100 probes
  a <- fit_within(d$M, d$exposure, d$subject, shrink_var = TRUE)
  b <- fit_within(d$M, d$exposure, d$subject, shrink_var = FALSE)
  expect_equal(a$beta, b$beta, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(a$se, b$se)))
})
