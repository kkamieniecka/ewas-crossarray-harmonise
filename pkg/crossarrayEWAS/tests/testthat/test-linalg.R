test_that("pinv satisfies the pseudo-inverse identity when the design is singular", {
  # After the within transform, a covariate constant within every subject
  # becomes exactly zero, so this is the shape the stage-04 driver hands it.
  X <- cbind(1, c(0, 1, 0, 1), c(2, 2, 2, 2))
  A <- crossprod(X)
  expect_lt(qr(A)$rank, ncol(A))
  P <- pinv(A)
  expect_equal(A %*% P %*% A, A, tolerance = 1e-8)
  expect_equal(P %*% A %*% P, P, tolerance = 1e-8)
})

test_that("pinv agrees with solve when the matrix is well conditioned", {
  set.seed(4)
  A <- crossprod(matrix(rnorm(25), 5))
  expect_equal(pinv(A), solve(A), tolerance = 1e-8)
})

test_that("pinv drops directions below the rcond cutoff", {
  # One singular value 1e-12 relative to the largest: kept at the numpy
  # default, discarded once rcond rises above it.
  A <- diag(c(1, 1e-12))
  expect_gt(max(abs(pinv(A))), 1e11)
  expect_equal(pinv(A, rcond = 1e-6), diag(c(1, 0)), tolerance = 1e-8)
})

test_that("demean_by_group leaves every group mean at zero", {
  set.seed(5)
  X <- matrix(rnorm(24), nrow = 8)
  g <- rep(c("a", "b", "c", "d"), each = 2)
  Z <- demean_by_group(X, g)
  means <- rowsum(Z, g) / as.vector(table(g))
  expect_equal(unname(as.matrix(means)), matrix(0, 4, ncol(X)),
               tolerance = 1e-12)
  # The transform is a projection: applying it twice changes nothing.
  expect_equal(demean_by_group(Z, g), Z, tolerance = 1e-12)
})

test_that("demean_rows_by_group centres each subject's samples in place", {
  set.seed(6)
  M <- matrix(rnorm(30), nrow = 5)            # probes x samples
  subj <- c("S1", "S1", "S1", "S2", "S2", "S2")
  R <- demean_rows_by_group(M, subj)
  for (s in unique(subj))
    expect_equal(rowMeans(R[, subj == s, drop = FALSE]),
                 rep(0, nrow(M)), tolerance = 1e-12)
  expect_equal(dim(R), dim(M))
  # A probe-specific constant is a subject-invariant shift, so it survives
  # the transform unchanged; a subject-specific one is annihilated.
  shift <- rnorm(nrow(M))
  expect_equal(demean_rows_by_group(M + shift, subj), R, tolerance = 1e-12)
  by_subject <- ifelse(subj == "S1", 3, -7)
  expect_equal(demean_rows_by_group(M + rep(by_subject, each = nrow(M)), subj),
               R, tolerance = 1e-12)
})

test_that("the two demean directions are transposes of each other", {
  set.seed(7)
  X <- matrix(rnorm(18), nrow = 6)
  g <- rep(c("p", "q"), each = 3)
  expect_equal(unname(demean_rows_by_group(t(X), g)),
               unname(t(demean_by_group(X, g))), tolerance = 1e-12)
})

test_that("demean_by_group keeps input row names but invents them when absent", {
  # Documented quirk rather than intent: with no row names on the input, the
  # group codes leak out of rowsum() into the result. Numerically harmless,
  # and the Python twin cannot show it, but a caller that relies on the names
  # would be misled. Fixing it means editing bin/ewasml.R, so the behaviour
  # is pinned here to make that a deliberate change.
  X <- matrix(as.double(1:6), ncol = 1)
  g <- rep(c("p", "q"), each = 3)
  expect_identical(rownames(demean_by_group(X, g)),
                   c("1", "1", "1", "2", "2", "2"))
  named <- X
  rownames(named) <- paste0("s", 1:6)
  expect_identical(rownames(demean_by_group(named, g)), paste0("s", 1:6))
})
