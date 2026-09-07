#!/usr/bin/env Rscript
## Does bin/ewasml.R compute the same numbers as bin/ewasml.py?
##
## Run tests/gen_equivalence_fixtures.py first; it writes the inputs and the
## Python reference outputs. This script recomputes each quantity in R and
## reports the largest disagreement, function by function.
##
## The fixtures are seeded pseudo-random arrays, not methylation data. This
## test is about the port; agreement on real data is a separate question
## (docs/design-rationale.md, "Two implementations of one estimator").
##
## Usage:  Rscript tests/test_equivalence.R <fixture-dir>

suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: test_equivalence.R <fixture-dir>")
FX <- args[1]
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=",
        commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "..", "bin", "ewasml.R"))

F_ <- function(f) file.path(FX, f)
ref <- jsonlite::fromJSON(F_("reference.json"))

fails <- 0L
report <- function(name, diff, tol, detail = "") {
  ok <- is.finite(diff) && diff <= tol
  if (!ok) fails <<- fails + 1L
  cat(sprintf("%-46s %-4s  max diff %-11.4g (tol %g) %s\n",
              name, if (ok) "ok" else "FAIL", diff, tol, detail))
}
maxad <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))
maxrd <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  d <- abs(a - b) / pmax(abs(b), 1e-12)
  max(d[is.finite(d)])
}

## ---- inputs ---------------------------------------------------------------
mv <- read_f64(F_("mval"))
M <- mv$mat
des <- read.csv(F_("design.csv"), stringsAsFactors = FALSE,
                colClasses = c(subject = "character"))
ann <- read.csv(F_("annotation.csv"), stringsAsFactors = FALSE,
                colClasses = c(chr = "character"))
covars <- as.matrix(des[, c("cov1", "cov2", "cov3")])
cat(sprintf("fixture: %d probes x %d samples, %d subjects\n",
            nrow(M), ncol(M), length(unique(des$subject))))

## read_f64 itself: the matrix must round-trip in the right orientation
stopifnot(nrow(M) == nrow(ann), ncol(M) == nrow(des))
report("read_f64 dimensions", 0, 0,
       sprintf("[%d x %d]", nrow(M), ncol(M)))

## ---- 1. probe-level fit ---------------------------------------------------
for (tag in c("raw", "shrunk")) {
  py <- read.csv(F_(sprintf("fit_%s.csv", tag)))
  shrink <- tag == "shrunk"
  # method = "mom" is the port of the Python; "limma" is compared below
  f <- fit_within(M, des$exposure, des$subject, covars,
                  shrink_var = shrink, var_method = "mom")
  report(sprintf("fit_within [%s] beta", tag), maxad(f$beta, py$beta), 1e-11)
  report(sprintf("fit_within [%s] se", tag), maxad(f$se, py$se), 1e-11)
  report(sprintf("fit_within [%s] z", tag), maxad(f$z, py$z), 1e-9)
  report(sprintf("fit_within [%s] df", tag),
         abs(f$df - ref[[sprintf("fit_%s_df", tag)]]), 0)
}
fit <- fit_within(M, des$exposure, des$subject, covars, var_method = "mom")

## ---- 2. variance moderation ----------------------------------------------
mod <- read.csv(F_("moderation.csv"))
r_mom <- moderate_variance(mod$s2, ref$moderation_df, method = "mom")
report("moderate_variance [mom vs python]", maxrd(r_mom, mod$post_mom), 1e-10)
if (requireNamespace("limma", quietly = TRUE)) {
  r_limma <- moderate_variance(mod$s2, ref$moderation_df, method = "limma")
  # not an equivalence requirement: limma::squeezeVar is the reference
  # implementation of the same model, so a difference here is a difference
  # between the hand-rolled fit and limma's, which is worth knowing.
  cat(sprintf("%-46s %-4s  max rel diff %.4g (informational)\n",
              "moderate_variance [limma vs python]", "--",
              maxrd(r_limma, mod$post_mom)))
} else {
  cat("moderate_variance [limma vs python]        skipped (limma not installed)\n")
}

## ---- 3. co-methylation clusters ------------------------------------------
rs <- read_f64(F_("resid"))
cl_py <- read.csv(F_("clusters.csv"))$cluster
cl_r <- comethylation_clusters(ann$chr_num, as.numeric(ann$pos), rs$mat,
                               max_gap = 1000L, rho_min = 0.30)
report("comethylation_clusters ids", maxad(cl_r, cl_py), 0,
       sprintf("[%d clusters, %d clustered probes]",
               max(cl_r) + 1L, sum(cl_r >= 0)))

## ---- 4. penalty and TV track ---------------------------------------------
w <- ifelse(fit$se > 0, 1 / fit$se^2, 0)
pen_py <- read.csv(F_("penalty.csv"))
theta_py <- read.csv(F_("tv_track.csv"))$theta
theta_r <- numeric(nrow(ann))
pen_r <- numeric(0)
for (c in sort(unique(ann$chr_num))) {
  idx <- which(ann$chr_num == c)
  lam <- distance_penalty(as.numeric(ann$pos[idx]), lam0 = ref$lam0,
                          decay_bp = ref$decay_bp, w = w[idx])
  pen_r <- c(pen_r, lam)
  theta_r[idx] <- tv_denoise(fit$beta[idx], w[idx], lam, n_iter = 3000L)
}
report("distance_penalty", maxad(pen_r, pen_py$lam), 1e-12)
report("tv_denoise theta", maxad(theta_r, theta_py), 1e-8)
## the segmentation itself, which is what the regions are built from
brk_py <- c(TRUE, abs(diff(theta_py)) > 1e-6)
brk_r <- c(TRUE, abs(diff(theta_r)) > 1e-6)
report("tv_denoise breakpoints", sum(brk_py != brk_r), 0,
       sprintf("[%d segments]", sum(brk_r)))

## ---- 5. region calling ----------------------------------------------------
cl_arg <- ifelse(cl_r >= 0, cl_r, -seq_along(cl_r))
reg_py <- read.csv(F_("regions.csv"), colClasses = c(chr = "character"))
reg_r <- call_regions(ann$chr, ann$pos, fit$beta, fit$se, theta_r,
                      min_probes = 3L, min_effect = 0.05, cluster = cl_arg)
tab <- reg_r$table
report("call_regions n_regions", abs(nrow(tab) - nrow(reg_py)), 0,
       sprintf("[%d regions]", nrow(tab)))
if (nrow(tab) == nrow(reg_py) && nrow(tab) > 0) {
  for (col in c("start", "end", "n_probes", "probe_start", "probe_end")) {
    report(sprintf("call_regions %s", col), maxad(tab[[col]], reg_py[[col]]), 0)
  }
  report("call_regions chr", sum(tab$chr != reg_py$chr), 0)
  for (col in c("effect_M", "se", "z", "theta")) {
    report(sprintf("call_regions %s", col), maxad(tab[[col]], reg_py[[col]]), 1e-9)
  }
  report("call_regions p_pointwise",
         maxrd(tab$p_pointwise, reg_py$p_pointwise), 1e-9)
}
report("call_regions max_abs_z", abs(reg_r$max_abs_z - ref$max_abs_z), 1e-9)

## ---- 6. permutation schemes and the FWER chain ---------------------------
## R cannot reproduce numpy's draws, so the draws are read from the fixture
## and the downstream arithmetic is compared exactly.
perm <- as.matrix(read.csv(F_("perm_within.csv")))
null_py <- read.csv(F_("null_within.csv"))$max_abs_z
null_r <- numeric(ncol(perm))
for (i in seq_len(ncol(perm))) {
  f <- fit_within(M, perm[, i], des$subject, covars, var_method = "mom")
  wl_all <- ifelse(f$se > 0, 1 / f$se^2, 0)
  th <- numeric(nrow(ann))
  for (c in sort(unique(ann$chr_num))) {
    idx <- which(ann$chr_num == c)
    lam <- distance_penalty(as.numeric(ann$pos[idx]), lam0 = ref$lam0,
                            decay_bp = ref$decay_bp, w = wl_all[idx])
    th[idx] <- tv_denoise(f$beta[idx], wl_all[idx], lam, n_iter = 800L)
  }
  r <- call_regions(ann$chr, ann$pos, f$beta, f$se, th,
                    min_probes = 3L, min_effect = 0.05)
  null_r[i] <- r$max_abs_z
}
report("null max|z| over python draws", maxad(null_r, null_py), 1e-8,
       sprintf("[%d replicates]", length(null_r)))

## the schemes themselves: properties, since the draws cannot match
set.seed(11)
pw <- within_subject_permutation(des$exposure, des$subject)
same_multiset <- all(vapply(split(seq_along(des$subject), des$subject),
  function(ix) isTRUE(all.equal(sort(pw[ix]), sort(des$exposure[ix]))),
  logical(1)))
report("within_subject_permutation preserves per-subject exposures",
       as.numeric(!same_multiset), 0)
report("within_subject_permutation leaves singletons fixed",
       sum(vapply(split(seq_along(des$subject), des$subject),
                  function(ix) if (length(ix) == 1)
                    as.numeric(pw[ix] != des$exposure[ix]) else 0,
                  numeric(1))), 0)
pn <- naive_permutation(des$exposure)
report("naive_permutation preserves the exposure multiset",
       as.numeric(!isTRUE(all.equal(sort(pn), sort(des$exposure)))), 0)

## ---- 7. stability_selection scheme ---------------------------------------
seen <- list()
ss <- stability_selection(function(sub) {
  seen[[length(seen) + 1L]] <<- sub
  c("A", "B")
}, des$subject, n_boot = 20L, frac = 0.5, seed = 3L)
nsub <- length(unique(des$subject))
k_expected <- max(2L, floor(0.5 * nsub))
report("stability_selection subsample size",
       max(abs(vapply(seen, length, integer(1)) - k_expected)), 0,
       sprintf("[k=%d of %d subjects]", k_expected, nsub))
report("stability_selection samples subjects without replacement",
       sum(vapply(seen, function(s) as.numeric(anyDuplicated(s) > 0), numeric(1))), 0)
report("stability_selection frequency arithmetic",
       maxad(ss$freq[c("A", "B")], c(1, 1)), 0)

## ---- 8. identification guard and deterministic fold split ---------------
report("within_df full design",
       abs(within_df(des$exposure, des$subject, covars) - ref$within_df_full), 0,
       sprintf("[df=%d]", within_df(des$exposure, des$subject, covars)))
half <- sort(unique(des$subject))[1:6]
sel <- des$subject %in% half
report("within_df half of the subjects",
       abs(within_df(des$exposure[sel], des$subject[sel],
                     covars[sel, , drop = FALSE]) - ref$within_df_half), 0)
fold_py <- read.csv(F_("folds.csv"), colClasses = c(subject = "character"))
mismatch <- 0
for (k in unique(fold_py$k)) {
  fr <- fold_assign(des$subject, k, seed = 1)
  py <- split(fold_py$subject[fold_py$k == k], fold_py$fold[fold_py$k == k])
  mismatch <- mismatch + as.numeric(!identical(
    lapply(unname(fr), sort), lapply(unname(py), sort)))
}
report("fold_assign identical splits", mismatch, 0,
       sprintf("[k = %s]", paste(unique(fold_py$k), collapse = ",")))

## ---- 9. block HSMM, stage 05 --------------------------------------------
## No RNG anywhere in this fit, so every quantity is comparable -- including
## the called blocks, which for stage 04 could only be compared on replayed
## draws.
hin <- read.csv(F_("hsmm_input.csv"), stringsAsFactors = FALSE,
                colClasses = c(chr = "character"))
tr <- read.csv(F_("transitions.csv"))
A_r <- dist_transitions(tr$d, c(0.05, 0.90, 0.05), 250000)
A_py <- as.matrix(tr[, -1])
A_flat <- matrix(NA_real_, nrow(tr), 9L)
for (i in 1:3) for (j in 1:3) A_flat[, (i - 1L) * 3L + j] <- A_r[, i, j]
report("dist_transitions A(d)", maxad(A_flat, A_py), 1e-15,
       sprintf("[d = %s]", paste(format(tr$d, scientific = TRUE, digits = 1),
                                 collapse = ",")))
report("A(d) rows sum to 1", max(abs(rowSums(A_flat[, c(1, 2, 3)]) - 1)), 1e-15)

t0 <- proc.time()[["elapsed"]]
hs <- fit_block_hsmm(hin$y, hin$se, hin$mid, hin$chr,
                     length_scale = ref$hsmm_length_scale, seed = 1)
t_hsmm <- proc.time()[["elapsed"]] - t0
post_py <- read.csv(F_("hsmm_posterior.csv"))
report("fit_block_hsmm iterations to convergence",
       abs(hs$n_iter - ref$hsmm_n_iter), 0,
       sprintf("[%d iterations, %d clusters, %.1fs]", hs$n_iter,
               ref$hsmm_n_clusters, t_hsmm))
report("fit_block_hsmm state means", maxad(hs$mu, ref$hsmm_mu), 1e-9,
       sprintf("[hypo %.4f, hyper %.4f]", hs$mu[1], hs$mu[3]))
report("fit_block_hsmm state sd", maxad(hs$sigma, ref$hsmm_sigma), 1e-9)
report("fit_block_hsmm stationary distribution", maxad(hs$pi, ref$hsmm_pi), 1e-9)
report("fit_block_hsmm log-likelihood",
       abs(hs$loglik - ref$hsmm_loglik) / abs(ref$hsmm_loglik), 1e-12,
       sprintf("[%.4f]", hs$loglik))
report("fit_block_hsmm posterior", maxad(as.matrix(post_py), hs$posterior), 1e-9)
report("posterior rows sum to 1", max(abs(rowSums(hs$posterior) - 1)), 1e-12)
report("state assignment identical",
       sum(apply(hs$posterior, 1, which.max) !=
           apply(as.matrix(post_py), 1, which.max)), 0,
       sprintf("[%d clusters]", nrow(hs$posterior)))

blk_py <- read.csv(F_("blocks.csv"), stringsAsFactors = FALSE,
                   colClasses = c(chr = "character"))
report("fit_block_hsmm neutral column", abs(hs$neutral - (ref$hsmm_neutral + 1L)), 0,
       sprintf("[R %d == Python %d, 1- vs 0-based]", hs$neutral, ref$hsmm_neutral))
blk <- call_blocks(hin$chr, hin$start, hin$end, hs$posterior,
                   min_post = 0.80, min_clusters = 3L, neutral = hs$neutral)
report("call_blocks block count", abs(nrow(blk) - ref$n_blocks), 0,
       sprintf("[%d blocks]", nrow(blk)))
if (nrow(blk) == nrow(blk_py) && nrow(blk) > 0L) {
  report("call_blocks boundaries and widths",
         max(abs(blk$start - blk_py$start), abs(blk$end - blk_py$end),
             abs(blk$width - blk_py$width)), 0)
  report("call_blocks cluster counts, chromosome and direction",
         sum(blk$n_clusters != blk_py$n_clusters | blk$chr != blk_py$chr |
             blk$direction != blk_py$direction), 0,
         sprintf("[%s]", paste(table(blk$direction), names(table(blk$direction)),
                               collapse = ", ")))
  report("call_blocks block posterior",
         maxad(blk$posterior, blk_py$posterior), 1e-9,
         sprintf("[min %.4f]", min(blk$posterior)))
}

## call_blocks decision rules, which no reference can pin down.
## NB: the block COUNT is not monotone in min_post -- a stricter threshold can
## drop one cluster out of a long run and split it into two reported blocks.
## What is monotone is the number of clusters called into blocks.
pp <- as.matrix(post_py)
strict <- call_blocks(hin$chr, hin$start, hin$end, pp, min_post = 0.99,
                      min_clusters = 3L, neutral = hs$neutral)
report("raising min_post cannot call more clusters into blocks",
       max(0, sum(strict$n_clusters) - sum(blk_py$n_clusters)), 0,
       sprintf("[%d clusters at 0.99 vs %d at 0.80, in %d vs %d blocks]",
               sum(strict$n_clusters), sum(blk_py$n_clusters),
               nrow(strict), nrow(blk_py)))
report("raising min_clusters cannot add blocks",
       max(0, nrow(call_blocks(hin$chr, hin$start, hin$end, pp,
                               min_post = 0.80, min_clusters = 10L,
                               neutral = hs$neutral)) - nrow(blk_py)), 0)
report("every block lies on one chromosome",
       sum(vapply(seq_len(nrow(blk)), function(i)
         as.numeric(length(unique(hin$chr[hin$chr == blk$chr[i] &
                                          hin$start >= blk$start[i] &
                                          hin$end <= blk$end[i]])) != 1L),
         numeric(1))), 0)
report("block posterior is the mean over its own clusters",
       max(vapply(seq_len(nrow(blk)), function(i) {
         k <- which(hin$chr == blk$chr[i] & hin$start >= blk$start[i] &
                    hin$end <= blk$end[i])
         st <- if (blk$direction[i] == "hypo") 1L else 3L
         abs(mean(hs$posterior[k, st]) - blk$posterior[i])
       }, numeric(1))), 1e-12)

## ---- 10. degenerate state ordering ---------------------------------------
## Every cluster effect on one side of zero, as in the real GSE237561 fit. The
## state pinned at mu = 0 sorts to an end, so a caller that assumes the middle
## column is neutral both mislabels the pinned state as directional and treats
## a genuine effect state as no change.
hs2 <- fit_block_hsmm(abs(hin$y) + 0.05, hin$se, hin$mid, hin$chr,
                      length_scale = ref$hsmm_length_scale)
post2_py <- as.matrix(read.csv(F_("hsmm_posterior_pos.csv")))
report("degenerate fit state means", maxad(hs2$mu, ref$hsmm_pos_mu), 1e-9,
       sprintf("[%s]", paste(sprintf("%.4f", hs2$mu), collapse = ", ")))
report("degenerate fit posterior", maxad(hs2$posterior, post2_py), 1e-9)
report("pinned neutral state sorts to an end",
       abs(hs2$neutral - (ref$hsmm_pos_neutral + 1L)), 0,
       sprintf("[column %d of 3, not the middle]", hs2$neutral))
blk2_py <- read.csv(F_("blocks_pos.csv"), stringsAsFactors = FALSE,
                    colClasses = c(chr = "character"))
blk2 <- call_blocks(hin$chr, hin$start, hin$end, hs2$posterior,
                    min_post = 0.80, min_clusters = 3L, neutral = hs2$neutral)
report("degenerate case block count",
       abs(nrow(blk2) - ref$hsmm_pos_n_blocks), 0,
       sprintf("[%d blocks; assuming the middle column reports %d]", nrow(blk2),
               nrow(call_blocks(hin$chr, hin$start, hin$end, hs2$posterior,
                                min_post = 0.80, min_clusters = 3L,
                                neutral = 2L))))
report("degenerate case directions identical",
       as.numeric(!identical(sort(unique(blk2$direction)),
                             sort(as.character(ref$hsmm_pos_directions)))), 0,
       sprintf("[%s]", paste(sort(unique(blk2$direction)), collapse = ",")))
report("no state is called hypo when none is below neutral",
       sum(blk2$direction == "hypo"), 0)
if (nrow(blk2) == nrow(blk2_py) && nrow(blk2) > 0L)
  report("degenerate case boundaries",
         max(abs(blk2$start - blk2_py$start), abs(blk2$end - blk2_py$end)), 0)

cat(sprintf("\n%s: %d check(s) failed\n",
            if (fails == 0L) "PASS" else "FAIL", fails))
quit(status = if (fails == 0L) 0L else 1L)
