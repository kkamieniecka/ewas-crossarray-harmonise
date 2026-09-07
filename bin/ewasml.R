## ewasml.R -- numerical core for the section 2.6 replacement (R port).
##
## This is a line-for-line port of bin/ewasml.py. The Python module remains in
## the tree until the two agree on GSE237561; tests/test_equivalence.R measures
## the agreement function by function. Where the two must differ, the reason is
## stated at the function.
##
## Dependencies: Matrix (banded solve), jsonlite (matrix header), base R
## (digamma/trigamma/pnorm). limma is optional and only for moderate_variance().
##
## Deliberate divergences from the Python, all verifiable in the test:
##  1. moderate_variance(method = "limma") calls limma::squeezeVar instead of
##     the hand-rolled method-of-moments fit. "mom" reproduces the Python.
##  2. tv_denoise() factorises the tridiagonal system once (CHOLMOD) instead of
##     calling LAPACK dgtsv every iteration. Same system, same fixed point,
##     ~1e-12 arithmetic difference; see the test for the observed effect on
##     segmentation.
##  3. Resampling functions use R's RNG, so permutation- and bootstrap-derived
##     numbers cannot be reproduced across languages. Only the schemes are
##     testable; the test drives them from Python-generated draws.

suppressPackageStartupMessages({
  library(Matrix)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# io
# ---------------------------------------------------------------------------

read_f64 <- function(stem) {
  d <- jsonlite::fromJSON(paste0(stem, "_dims.json"))
  n <- as.integer(d$nrow) * as.integer(d$ncol)
  x <- readBin(paste0(stem, ".f64"), "double", n = n, size = 8,
               endian = "little")
  if (length(x) != n)
    stop(sprintf("%s.f64 has %d values, header says %dx%d",
                 stem, length(x), d$nrow, d$ncol))
  # the dump is column-major, which is R's native layout
  mat <- matrix(x, nrow = as.integer(d$nrow), ncol = as.integer(d$ncol))
  rn <- if (!is.null(d$rownames)) as.character(d$rownames) else NULL
  cn <- if (!is.null(d$colnames)) as.character(d$colnames) else NULL
  list(mat = mat, rownames = rn, colnames = cn)
}

# ---------------------------------------------------------------------------
# probe-level within-subject estimator
# ---------------------------------------------------------------------------

#' Moore-Penrose pseudo-inverse with numpy.linalg.pinv's cutoff rule.
pinv <- function(A, rcond = 1e-15) {
  s <- svd(A)
  cutoff <- rcond * max(s$d)
  inv_d <- ifelse(s$d > cutoff, 1 / s$d, 0)
  s$v %*% (inv_d * t(s$u))
}

#' Subtract the group mean from every column (the within transform).
demean_by_group <- function(X, groups) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  g <- match(groups, sort(unique(groups)))
  cnt <- tabulate(g, nbins = max(g))
  means <- rowsum(X, g, reorder = TRUE) / cnt
  X - means[g, , drop = FALSE]
}

#' Row-wise (per probe) subtraction of the group mean across samples.
demean_rows_by_group <- function(M, groups) {
  g <- match(groups, sort(unique(groups)))
  ng <- max(g)
  cnt <- tabulate(g, nbins = ng)
  # samples x groups averaging operator; sparse so this stays O(n_samples)
  G <- sparseMatrix(i = seq_along(g), j = g, x = 1 / cnt[g],
                    dims = c(length(g), ng))
  means <- as.matrix(M %*% G)          # probes x groups
  M - means[, g, drop = FALSE]
}

#' Solve trigamma(y) = x (Newton, as in limma's trigammaInverse).
trigamma_inv <- function(x) {
  y <- 0.5 + 1 / x
  for (i in seq_len(60)) {
    f <- psigamma(y, 1) - x
    d <- psigamma(y, 2)
    if (d == 0) break
    y_new <- y - f / d
    if (y_new <= 0) y_new <- y / 2
    if (abs(y_new - y) < 1e-10 * max(1, y)) {
      y <- y_new
      break
    }
    y <- y_new
  }
  y
}

#' Empirical-Bayes shrinkage of per-probe residual variance.
#'
#' method = "mom"   : the method-of-moments fit ported from ewasml.py
#' method = "limma" : limma::squeezeVar, i.e. the reference implementation of
#'                    the same model. Preferred in R -- it is maintained,
#'                    handles the degenerate cases, and this repository already
#'                    depends on limma at stage 02.
moderate_variance <- function(s2, df, method = c("limma", "mom")) {
  method <- match.arg(method)
  ok <- is.finite(s2) & s2 > 0
  if (sum(ok) < 100) return(s2)

  if (method == "limma") {
    if (!requireNamespace("limma", quietly = TRUE))
      stop("moderate_variance(method='limma') needs limma; ",
           "install it or pass method='mom'")
    out <- s2
    out[ok] <- limma::squeezeVar(s2[ok], df = df)$var.post
    return(out)
  }

  e <- log(s2[ok]) - digamma(df / 2) + log(df / 2)
  ebar <- mean(e)
  target <- var(e) - psigamma(df / 2, 1)      # var() is the ddof=1 estimator
  d0 <- if (target <= 0) Inf else 2 * trigamma_inv(target)
  if (!is.finite(d0)) {
    s0sq <- exp(ebar + digamma(df / 2) - log(df / 2))
    return(rep(s0sq, length(s2)))
  }
  s0sq <- exp(ebar + digamma(d0 / 2) - log(d0 / 2))
  out <- s2
  out[ok] <- (d0 * s0sq + df * s2[ok]) / (d0 + df)
  out
}

FNV_OFFSET <- 2166136261
FNV_PRIME <- 16777619

#' 32-bit FNV-1a. Small, exactly specified, and reproducible in any language
#' -- which is the only reason to hand-roll a hash here. The multiply is split
#' into 16-bit halves because h * FNV_PRIME exceeds the 53-bit exact range of
#' a double.
fnv1a <- function(s) {
  h <- FNV_OFFSET
  for (b in utf8ToInt(enc2utf8(s))) {
    # h can exceed 2^31, which bitwXor cannot represent; the byte only
    # touches the low 8 bits, so XOR there and leave the rest alone
    low8 <- h %% 256
    h <- h - low8 + bitwXor(as.integer(low8), as.integer(b))
    lo <- h %% 65536; hi <- h %/% 65536
    h <- (lo * FNV_PRIME + ((hi * FNV_PRIME) %% 65536) * 65536) %% 4294967296
  }
  h
}

#' Assign subjects to cross-validation folds deterministically.
#'
#' Subjects are ordered by a hash of (seed, subject id) and dealt into
#' contiguous folds. This replaces a seeded RNG shuffle so that the fold split
#' -- and therefore the selected smoothness -- is identical in the R and
#' Python implementations and stable across library versions. Changing `seed`
#' gives an independent split.
fold_assign <- function(subject, n_folds, seed = 0) {
  uniq <- sort(unique(as.character(subject)))
  hv <- vapply(uniq, function(s) fnv1a(paste0(seed, ":", s)), numeric(1))
  keys <- uniq[order(hv, uniq)]
  n <- length(keys); k <- max(1L, as.integer(n_folds))
  sizes <- rep(n %/% k, k)
  if (n %% k) sizes[seq_len(n %% k)] <- sizes[seq_len(n %% k)] + 1L
  out <- split(keys, rep(seq_len(k), times = sizes))
  out[lengths(out) > 0]
}

#' Residual degrees of freedom the within-subject fit would have.
#'
#' n_samples - n_subjects - rank(within design). Returns 0 when the exposure
#' itself is annihilated by the within transform, i.e. when the longitudinal
#' contrast is not identified at all. Lets a caller decide whether a
#' cross-validation fold or a subsample is estimable before paying for the
#' fit -- a fold of mostly single-visit subjects is not, and this is the cheap
#' test for it.
within_df <- function(exposure, subject, covars = NULL) {
  X <- matrix(as.double(exposure), ncol = 1)
  if (!is.null(covars) && length(covars)) X <- cbind(X, as.matrix(covars))
  Xw <- demean_by_group(X, subject)
  keep <- apply(abs(Xw), 2, max) > 1e-9
  if (!keep[1]) return(0L)
  as.integer(length(subject) - length(unique(subject)) - sum(keep))
}

#' Within-subject (fixed-effects) estimate of the exposure effect per probe.
#'
#' M         probes x samples matrix of M-values
#' exposure  length-samples exposure
#' subject   length-samples subject identifier
#' covars    samples x k TIME-VARYING covariates only; the within transform
#'           removes time-invariant ones exactly and including them makes the
#'           design rank-deficient.
#'
#' Subject intercepts are eliminated by demeaning rather than estimated, which
#' is what makes the permutation and subsampling loops affordable.
fit_within <- function(M, exposure, subject, covars = NULL,
                       shrink_var = TRUE, var_method = "limma") {
  X <- matrix(as.double(exposure), ncol = 1)
  if (!is.null(covars) && length(covars))
    X <- cbind(X, as.matrix(covars))
  Xw <- demean_by_group(X, subject)
  keep <- apply(abs(Xw), 2, max) > 1e-9
  if (!keep[1])
    stop("exposure has no within-subject variation; a longitudinal contrast ",
         "is not identified")
  Xw <- Xw[, keep, drop = FALSE]
  Mw <- demean_rows_by_group(M, subject)

  XtX <- crossprod(Xw)
  XtXi <- pinv(XtX)
  B <- Mw %*% Xw %*% XtXi                    # probes x p
  resid <- Mw - B %*% t(Xw)
  n_sub <- length(unique(subject))
  df <- ncol(M) - n_sub - ncol(Xw)
  if (df <= 0) stop("no residual degrees of freedom")
  s2 <- rowSums(resid^2) / df
  if (shrink_var) s2 <- moderate_variance(s2, df, method = var_method)
  v_beta <- XtXi[1, 1]
  se <- sqrt(s2 * v_beta)
  beta <- B[, 1]
  z <- ifelse(se > 0, beta / se, 0)
  list(beta = beta, se = se, z = z, df = df,
       n_subjects = n_sub, n_samples = ncol(M))
}

# ---------------------------------------------------------------------------
# replaces clusterMaker: data-driven co-methylation clustering
# ---------------------------------------------------------------------------

#' Cluster probes that are both physically close and empirically co-methylated.
#'
#' Returns an integer cluster id per probe (-1 = singleton). `resid` is the
#' within-subject residual matrix, so the correlation is between-visit
#' co-variation -- not between-subject or between-array variation, which is
#' what a region is supposed to mean.
#'
#' Only genomically adjacent pairs are ever linked, so the union-find in the
#' Python reduces to run-length grouping over the sorted order. Cluster ids are
#' assigned in genomic order; the Python numbers them by the smallest input
#' index in each group, which is the same thing whenever the input is sorted
#' (stage 04 always sorts).
comethylation_clusters <- function(chrom, pos, resid, max_gap = 1000L,
                                   rho_min = 0.30, chunk = 50000L) {
  n <- length(pos)
  if (n == 0) return(integer(0))
  ord <- order(chrom, pos)

  # normalise residual rows once so correlation is a dot product
  R <- resid - rowMeans(resid)
  nrm <- sqrt(rowSums(R^2))
  nrm[nrm == 0] <- Inf
  R <- R / nrm

  i <- ord[-n]; j <- ord[-1]
  same_chr <- chrom[i] == chrom[j]
  near <- (pos[j] - pos[i]) <= max_gap
  cand <- which(same_chr & near)

  rho <- rep(-Inf, n - 1)
  if (length(cand)) {
    # chunked so the two gathered copies never dominate memory
    for (st in seq(1L, length(cand), by = chunk)) {
      sl <- cand[st:min(st + chunk - 1L, length(cand))]
      rho[sl] <- rowSums(R[i[sl], , drop = FALSE] * R[j[sl], , drop = FALSE])
    }
  }
  link <- rep(FALSE, n - 1)
  link[cand] <- rho[cand] >= rho_min

  grp <- cumsum(c(TRUE, !link))          # group id along the sorted order
  sz <- tabulate(grp)
  cid_sorted <- ifelse(sz[grp] >= 2, grp, -1L)
  keep <- sort(unique(cid_sorted[cid_sorted > 0]))
  remap <- integer(0)
  if (length(keep)) {
    remap <- seq_along(keep) - 1L
    names(remap) <- as.character(keep)
    cid_sorted <- ifelse(cid_sorted > 0,
                         remap[as.character(cid_sorted)], -1L)
  }
  out <- integer(n)
  out[ord] <- as.integer(cid_sorted)
  out
}

# ---------------------------------------------------------------------------
# replaces loess smoothing: weighted total-variation denoising
# ---------------------------------------------------------------------------

soft_threshold <- function(x, t) sign(x) * pmax(abs(x) - t, 0)

#' Weighted total-variation denoising by ADMM.
#'
#' Solves  min_theta 0.5*sum_i w_i (y_i - theta_i)^2 + sum_i lam_i |d theta_i|
#' with y in genomic order within one contiguous run.
#'
#' lam is length n-1. A distance-decayed lam is what makes this a *genomic*
#' smoother: a large gap gets lam ~ 0 and is free to break, which is precisely
#' the behaviour loess with a fixed span cannot express.
#'
#' The returned track is DEBIASED: the l1 penalty chooses where the breakpoints
#' are, and each segment is then reported as the precision-weighted mean of y
#' over that segment (the relaxed-lasso convention). The reported region effect
#' is therefore an unbiased weighted mean in M-value units, directly comparable
#' with a per-probe coefficient.
#'
#' Unlike the Python, the tridiagonal system is factorised once: (W + rho D'D)
#' does not change across iterations, so re-solving it 3000 times with a fresh
#' LAPACK call is wasted work.
tv_denoise <- function(y, w, lam, n_iter = 3000L, rho = NULL,
                       tol = 1e-6, over_relax = 1.7) {
  n <- length(y)
  if (n == 1) return(as.double(y))
  y <- as.double(y); w <- as.double(w); lam <- as.double(lam)
  if (is.null(rho)) {
    pos_w <- w[w > 0]
    rho <- if (length(pos_w)) median(pos_w) else 1
    rho <- max(rho, 1e-8)
  }

  z <- numeric(n - 1)
  u <- numeric(n - 1)

  diag_DtD <- c(1, rep(2, max(n - 2, 0)), 1)
  A <- bandSparse(n, n, k = c(0, 1),
                  diagonals = list(w + rho * diag_DtD, rep(-rho, n - 1)),
                  symmetric = TRUE)
  ch <- Cholesky(A, perm = FALSE, LDL = FALSE)

  Wy <- w * y
  eps_abs <- tol * sqrt(n)
  theta <- y
  for (it in seq_len(n_iter)) {
    v <- z - u
    rhs <- Wy
    rhs[-n] <- rhs[-n] - rho * v
    rhs[-1] <- rhs[-1] + rho * v
    theta <- as.vector(solve(ch, rhs))
    dtheta <- diff(theta)
    dhat <- over_relax * dtheta + (1 - over_relax) * z
    z_old <- z
    z <- soft_threshold(dhat + u, lam / rho)
    u <- u + dhat - z
    r_prim <- sqrt(sum((dtheta - z)^2))
    r_dual <- rho * sqrt(sum((z - z_old)^2))
    if (r_prim < eps_abs + tol * sqrt(sum(z^2)) &&
        r_dual < eps_abs + tol * rho * sqrt(sum(u^2)))
      break
  }
  # The primal iterate theta is only approximately piecewise constant, but the
  # split variable z IS exactly sparse (it comes straight out of a soft
  # threshold), so z defines the breakpoints. Thresholding |diff(theta)| would
  # make the segmentation depend on the scale of y.
  seg <- cumsum(c(0, z != 0)) + 1
  sums <- rowsum(cbind(w * y, w), seg, reorder = TRUE)
  wsum <- sums[, 2]
  val <- ifelse(wsum > 0, sums[, 1] / wsum, NA_real_)
  if (anyNA(val)) {
    # zero-precision segment: fall back to the primal iterate's mean
    fb <- rowsum(theta, seg, reorder = TRUE)[, 1] / tabulate(seg)
    val[is.na(val)] <- fb[is.na(val)]
  }
  as.vector(val[match(seg, sort(unique(seg)))])
}

#' lam_i = lam0 * scale * exp(-d_i / decay_bp) over successive distances.
#'
#' `scale` is median(w) when precision weights are supplied, which makes lam0
#' DIMENSIONLESS: lam0 ~ 1 means "one probe's worth of precision resists one
#' unit of jump", so the same grid is usable on M-values, beta values or a
#' different cohort. Without it lam0 must be retuned whenever the noise level
#' changes -- the practical reason a fixed loess span does not transfer between
#' studies.
distance_penalty <- function(pos, lam0, decay_bp = 1000, w = NULL,
                             lam_floor = 0) {
  d <- diff(as.double(pos))
  scale <- 1
  if (!is.null(w)) {
    pos_w <- w[w > 0]
    if (length(pos_w)) scale <- median(pos_w)
  }
  pmax(lam0 * scale * exp(-d / decay_bp), lam_floor)
}

# ---------------------------------------------------------------------------
# region calling from a piecewise-constant track
# ---------------------------------------------------------------------------

#' Turn the denoised track into regions.
#'
#' Segment boundaries are where the piecewise-constant solution jumps (or the
#' chromosome/cluster changes), so boundaries are estimated, not thresholded
#' off a smooth curve. The region statistic is the precision-weighted mean
#' effect over its own standard error, which accounts for probes in a region
#' carrying very different amounts of information.
call_regions <- function(chrom, pos, beta, se, theta, min_probes = 3L,
                         min_effect = 0.05, tol = 1e-6, cluster = NULL) {
  n <- length(beta)
  w <- ifelse(se > 0, 1 / se^2, 0)
  same_chr <- chrom[-1] == chrom[-n]
  jump <- abs(diff(theta)) > tol
  same_cl <- if (is.null(cluster)) rep(TRUE, n - 1) else
    cluster[-1] == cluster[-n]
  brk <- c(TRUE, (!same_chr) | jump | (!same_cl))
  seg <- cumsum(brk)

  np <- tabulate(seg)
  agg <- rowsum(cbind(w * beta, w), seg, reorder = TRUE)
  wsum <- agg[, 2]
  bbar <- agg[, 1] / wsum
  se_bar <- sqrt(1 / wsum)
  zbar <- ifelse(se_bar > 0, bbar / se_bar, 0)

  first <- which(brk)
  last <- c(first[-1] - 1L, n)
  ok <- np >= min_probes & wsum > 0 & abs(bbar) >= min_effect
  ok[is.na(ok)] <- FALSE

  if (!any(ok))
    return(list(table = data.frame(), max_abs_z = 0))

  k <- which(ok)
  tab <- data.frame(
    chr = as.character(chrom[first[k]]),
    start = as.integer(pos[first[k]]),
    end = as.integer(pos[last[k]]),
    width = as.integer(pos[last[k]] - pos[first[k]] + 1),
    n_probes = as.integer(np[k]),
    effect_M = bbar[k],
    se = se_bar[k],
    z = zbar[k],
    p_pointwise = 2 * pnorm(-abs(zbar[k])),
    theta = theta[first[k]],
    probe_start = as.integer(first[k] - 1L),   # 0-based, as in the Python
    probe_end = as.integer(last[k] - 1L),
    stringsAsFactors = FALSE
  )
  list(table = tab, max_abs_z = max(abs(zbar[k])))
}

# ---------------------------------------------------------------------------
# inference: the two permutation schemes
# ---------------------------------------------------------------------------

#' Shuffle the exposure among the visits of the same subject.
#'
#' This is the null that matches a longitudinal design: subject identity, array
#' type, cohort, chip, number of visits and the within-subject covariance are
#' all untouched, and only the pairing between visit and exposure is broken. A
#' subject with a single usable visit contributes nothing under this null,
#' exactly as it contributes nothing to the within-subject estimator.
within_subject_permutation <- function(exposure, subject) {
  out <- as.double(exposure)
  for (s in unique(subject)) {
    idx <- which(subject == s)
    if (length(idx) > 1) out[idx] <- sample(out[idx])
  }
  out
}

#' Free permutation of the exposure across all samples -- the scheme bumphunter
#' uses on the design column. Retained only so the pipeline can demonstrate
#' that it is mis-calibrated here (it breaks the subject pairing and reassigns
#' exposures across array types).
naive_permutation <- function(exposure, subject = NULL) {
  sample(as.double(exposure))
}

# ---------------------------------------------------------------------------
# stability selection
# ---------------------------------------------------------------------------

#' Meinshausen-Buhlmann stability selection at the SUBJECT level.
#'
#' fit_fn(subject_subset) must return a character vector of region keys
#' selected on that subsample. Subsampling subjects (not samples) is required:
#' resampling samples would split a subject's visit series across train and
#' test and leak the within-subject effect being tested.
stability_selection <- function(fit_fn, subjects, n_boot = 100L, frac = 0.5,
                                seed = 1L) {
  set.seed(seed)
  uniq <- sort(unique(subjects))
  k <- max(2L, floor(frac * length(uniq)))
  counts <- list()
  for (b in seq_len(n_boot)) {
    sub <- sample(uniq, size = k, replace = FALSE)
    for (key in unique(fit_fn(sub)))
      counts[[key]] <- (if (is.null(counts[[key]])) 0 else counts[[key]]) + 1
  }
  list(freq = vapply(counts, function(v) v / n_boot, numeric(1)),
       n_boot = n_boot)
}
