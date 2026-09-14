test_that("read_f64 round-trips a matrix with dimnames", {
  stem <- tempfile()
  m <- matrix(as.double(1:6), nrow = 2,
              dimnames = list(c("cg1", "cg2"), c("s1", "s2", "s3")))
  writeBin(as.vector(m), paste0(stem, ".f64"))
  writeLines(jsonlite::toJSON(list(nrow = 2L, ncol = 3L,
                                   rownames = rownames(m),
                                   colnames = colnames(m)),
                              auto_unbox = TRUE),
             paste0(stem, "_dims.json"))
  got <- read_f64(stem)
  expect_equal(got$mat, matrix(as.double(1:6), nrow = 2))
  expect_identical(got$rownames, c("cg1", "cg2"))
  expect_identical(got$colnames, c("s1", "s2", "s3"))
})
