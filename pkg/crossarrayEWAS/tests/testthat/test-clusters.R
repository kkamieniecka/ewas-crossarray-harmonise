test_that("only genomically adjacent, co-varying probes cluster", {
  set.seed(3)
  shared <- rnorm(10)
  resid <- rbind(shared + rnorm(10, sd = 0.05),     # probes 1-3 co-vary
                 shared + rnorm(10, sd = 0.05),
                 shared + rnorm(10, sd = 0.05),
                 rnorm(10), rnorm(10), rnorm(10))
  chrom <- rep("chr1", 6)
  pos <- as.integer(c(100, 200, 300, 400, 500, 600))
  cid <- comethylation_clusters(chrom, pos, resid, max_gap = 150L)
  expect_length(cid, 6)
  expect_identical(length(unique(cid[1:3])), 1L)     # one cluster
  expect_false(cid[1] == cid[6])
})

test_that("a gap wider than max_gap breaks the cluster", {
  set.seed(4)
  shared <- rnorm(10)
  resid <- rbind(shared + rnorm(10, sd = 0.05), shared + rnorm(10, sd = 0.05))
  far <- comethylation_clusters(rep("chr1", 2), as.integer(c(100, 50000)),
                                resid, max_gap = 1000L)
  near <- comethylation_clusters(rep("chr1", 2), as.integer(c(100, 300)),
                                 resid, max_gap = 1000L)
  expect_false(far[1] == far[2] && far[1] > 0)
  expect_identical(near[1], near[2])
})

test_that("a chromosome change breaks the cluster", {
  set.seed(5)
  shared <- rnorm(10)
  resid <- rbind(shared + rnorm(10, sd = 0.05), shared + rnorm(10, sd = 0.05))
  cid <- comethylation_clusters(c("chr1", "chr2"), as.integer(c(100, 200)),
                                resid, max_gap = 1000L)
  expect_false(cid[1] == cid[2] && cid[1] > 0)
})
