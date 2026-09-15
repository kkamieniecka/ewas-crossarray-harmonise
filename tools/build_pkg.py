#!/usr/bin/env python3
"""Split bin/ewasml.R into the crossarrayEWAS package skeleton.

Function bodies are copied verbatim; only the documentation block above each
function is rewritten into roxygen form, keeping the original prose as
@details.
"""
import os, re, shutil, textwrap

SRC = os.environ.get("EWASML_SRC", "bin/ewasml.R")
PKG = "pkg/crossarrayEWAS"
PKGNAME = "crossarrayEWAS"

# function -> R/ file
LAYOUT = {
    "read_f64": "io",
    "pinv": "linalg", "demean_by_group": "linalg",
    "demean_rows_by_group": "linalg", "trigamma_inv": "linalg",
    "moderate_variance": "variance",
    "fnv1a": "folds", "fold_assign": "folds", "within_df": "folds",
    "fit_within": "within",
    "comethylation_clusters": "clusters",
    "soft_threshold": "denoise", "tv_denoise": "denoise",
    "distance_penalty": "denoise",
    "call_regions": "regions",
    "within_subject_permutation": "resample",
    "naive_permutation": "resample", "stability_selection": "resample",
    "dist_transitions": "blocks", "fit_block_hsmm": "blocks",
    "state_labels": "blocks", "call_blocks": "blocks",
}
FILE_ORDER = ["io", "linalg", "variance", "folds", "within", "clusters",
              "denoise", "regions", "resample", "blocks"]

# Everything the stage drivers call has to be exported, or they cannot run
# against the installed package. pinv() is exported for exactly that reason:
# bin/04_dmr_ml.R uses it to residualise the within-transformed design.
INTERNAL = {"trigamma_inv", "fnv1a", "soft_threshold", "dist_transitions"}

# shared example fixtures, inlined per example so each is self-contained
FIT = """set.seed(1)
n_probes <- 20; n_sub <- 6
subject <- rep(sprintf("S%02d", seq_len(n_sub)), each = 2)
exposure <- rep(c(0, 1), times = n_sub)
M <- matrix(rnorm(n_probes * length(subject)), nrow = n_probes)
M[1:3, exposure == 1] <- M[1:3, exposure == 1] + 1.2"""

TRACK = """chrom <- rep("chr1", 12)
pos <- as.integer(seq(1000, 3200, length.out = 12))
beta <- c(rep(0.4, 6), rep(-0.1, 6))
se <- rep(0.1, 12)
lam <- distance_penalty(pos, lam0 = 1, w = 1 / se^2)
theta <- tv_denoise(beta, w = 1 / se^2, lam = lam, n_iter = 200L)"""

HSMM = """set.seed(1)
chrom <- rep("chr1", 40)
pos <- as.integer(seq(1e5, 5e6, length.out = 40))
y <- c(rnorm(15, 0.5, 0.1), rnorm(25, 0, 0.1))
se <- rep(0.1, 40)
fit <- fit_block_hsmm(y, se, pos, chrom, n_iter = 20L)"""

DOCS = {
 "read_f64": dict(
  title="Read a raw float64 matrix written by the pipeline",
  desc=("Reads the ``<stem>.f64`` / ``<stem>_dims.json`` pair that the "
        "Nextflow stages use to pass probe-level matrices between processes. "
        "This is a transport format, not an analysis format: it exists because "
        "a column-major double dump plus a JSON header round-trips between R "
        "and Python without either language's serialisation getting involved."),
  params=[("stem", "path prefix without extension; ``<stem>.f64`` holds the "
           "doubles and ``<stem>_dims.json`` the dimensions and dimnames.")],
  ret=("list with ``mat`` (numeric matrix, probes x samples), ``rownames`` and "
       "``colnames`` (character or ``NULL``)."),
  ex='tmp <- tempfile()\nm <- matrix(as.double(1:6), nrow = 2)\nwriteBin(as.vector(m), paste0(tmp, ".f64"))\nwriteLines(jsonlite::toJSON(list(nrow = 2L, ncol = 3L), auto_unbox = TRUE),\n           paste0(tmp, "_dims.json"))\nread_f64(tmp)$mat'),
 "pinv": dict(
   title="Moore-Penrose pseudo-inverse of a small dense matrix",
   desc=("Singular values at or below ``rcond * max(d)`` are treated as zero "
         "rather than inverted, so a rank-deficient cross-product matrix "
         "gives the least-squares solution instead of an error. This is what "
         "makes the within-subject design safe to residualise against: after "
         "the within transform, a covariate that is constant within every "
         "subject becomes exactly zero and leaves the design singular.\n"
         "\n"
         "Intended for the few-column design matrices of a within-subject "
         "fit, where a dense SVD costs nothing. It is not a substitute for a "
         "sparse solver on probe-scale matrices."),
   params=[("A", "Numeric matrix to invert. Dense; no symmetry assumed."),
           ("rcond", "Relative singular-value cutoff. Singular values at or "
                     "below ``rcond * max(d)`` are inverted to zero. The "
                     "default, ``1e-15``, matches ``numpy.linalg.pinv`` so "
                     "the R and Python cores agree.")],
   ret="Numeric matrix, the pseudo-inverse of ``A``.",
   ex=('X <- cbind(1, c(0, 1, 0, 1), c(2, 2, 2, 2))  # third column redundant\n'
       'A <- crossprod(X)\n'
       'qr(A)$rank            # 2, not 3\n'
       'P <- pinv(A)\n'
       'all.equal(A %*% P %*% A, A)  # the defining property still holds')),
 "trigamma_inv": dict(title="Invert the trigamma function", desc="", params=[], ret="", ex=""),
 "fnv1a": dict(title="32-bit FNV-1a hash of a string", desc="", params=[], ret="", ex=""),
 "soft_threshold": dict(title="Soft-thresholding operator", desc="", params=[], ret="", ex=""),
 "dist_transitions": dict(title="Distance-dependent transition array", desc="", params=[], ret="", ex=""),
 "demean_by_group": dict(
  title="Subtract group means from a design matrix (the within transform)",
  desc=("Applies the within-subject transform to a samples x terms design "
        "matrix, which is how subject intercepts are eliminated rather than "
        "estimated."),
  params=[("X", "numeric matrix, samples x terms."),
          ("groups", "grouping vector of length ``nrow(X)``, normally subject "
           "identifiers.")],
  ret=("numeric matrix of the same dimensions, each group's mean removed. "
       "Row names are carried over from ``X`` when it has them; when it does "
       "not, the returned matrix carries the integer group codes as row "
       "names, which is an artefact of the grouped-sum step rather than "
       "information about the samples."),
  ex='X <- cbind(exposure = c(0, 1, 0, 1), age = c(40, 41, 55, 56))\ndemean_by_group(X, groups = c("A", "A", "B", "B"))'),
 "demean_rows_by_group": dict(
  title="Subtract group means from every probe of a methylation matrix",
  desc=("Row-wise (per probe) within-subject transform of a probes x samples "
        "matrix. Uses a sparse averaging operator so cost stays linear in the "
        "number of samples."),
  params=[("M", "numeric matrix, probes x samples."),
          ("groups", "grouping vector of length ``ncol(M)``, normally subject "
           "identifiers.")],
  ret="numeric matrix of the same dimensions, each group's mean removed per row.",
  ex='M <- matrix(as.double(1:8), nrow = 2)\ndemean_rows_by_group(M, groups = c("A", "A", "B", "B"))'),
 "moderate_variance": dict(
  title="Empirical-Bayes shrinkage of per-probe residual variance",
  desc=("Shrinks per-probe residual variances towards a fitted prior, which is "
        "what keeps a probe with a fluke-small variance out of the top of the "
        "region ranking. With fewer than 100 finite positive variances the "
        "prior cannot be fitted and the input is returned unchanged."),
  params=[("s2", "numeric vector of per-probe residual variances."),
          ("df", "residual degrees of freedom of the fit that produced ``s2`` "
           "(scalar)."),
          ("method", "``\"limma\"`` for limma::squeezeVar (default, preferred "
           "in R), ``\"mom\"`` for the method-of-moments fit ported from the "
           "Python reference.")],
  ret="numeric vector of shrunken variances, same length as ``s2``.",
  ex='set.seed(1)\ns2 <- rchisq(300, df = 8) / 8\nsummary(moderate_variance(s2, df = 8, method = "mom"))'),
 "fold_assign": dict(
  title="Deterministic assignment of subjects to cross-validation folds",
  desc=("Hash-ordered fold assignment, so the split - and therefore the "
        "selected smoothness - is identical in R and Python and stable across "
        "library versions."),
  params=[("subject", "vector of subject identifiers; duplicates are collapsed "
           "so a subject's visits never straddle two folds."),
          ("n_folds", "number of folds."),
          ("seed", "integer mixed into the hash; a different value gives an "
           "independent split.")],
  ret="list of character vectors, one per non-empty fold, holding subject ids.",
  ex='fold_assign(sprintf("S%02d", 1:10), n_folds = 3)'),
 "within_df": dict(
  title="Residual degrees of freedom of the within-subject fit",
  desc=("Cheap identifiability test for a fold or a subsample: returns 0 when "
        "the within transform annihilates the exposure, i.e. when the "
        "longitudinal contrast is not identified at all."),
  params=[("exposure", "numeric exposure, one value per sample."),
          ("subject", "subject identifier, one per sample."),
          ("covars", "optional samples x k matrix of time-varying covariates.")],
  ret="integer scalar; 0 if the exposure is not identified within subject.",
  ex='within_df(exposure = c(0, 1, 0, 1), subject = c("A", "A", "B", "B"))\nwithin_df(exposure = c(0, 0, 1, 1), subject = c("A", "A", "B", "B"))'),
 "fit_within": dict(
  title="Within-subject (fixed-effects) probe-level fit",
  desc=("Fits the exposure effect within subject for every probe at once, "
        "eliminating subject intercepts by demeaning rather than estimating "
        "them, which is what makes the permutation and subsampling loops "
        "affordable."),
  params=[("M", "numeric matrix, probes x samples (M-values)."),
          ("exposure", "numeric exposure, one value per sample."),
          ("subject", "subject identifier, one per sample."),
          ("covars", "optional samples x k matrix of TIME-VARYING covariates "
           "only; time-invariant ones are removed exactly by the transform and "
           "including them makes the design rank-deficient."),
          ("shrink_var", "apply empirical-Bayes variance moderation."),
          ("var_method", "method passed to ``moderate_variance()``.")],
  ret=("list with ``beta``, ``se``, ``z`` (one per probe), ``df``, "
       "``n_subjects`` and ``n_samples``."),
  ex=FIT + '\nfit <- fit_within(M, exposure, subject)\nstr(fit[c("df", "n_subjects", "n_samples")])\nhead(fit$z, 3)'),
 "comethylation_clusters": dict(
  title="Group adjacent probes that co-vary within subject",
  desc=("Clusters genomically adjacent probes whose within-subject residuals "
        "correlate, so a region means co-variation of the longitudinal signal "
        "rather than between-subject or between-array variation."),
  params=[("chrom", "chromosome per probe."),
          ("pos", "genomic position per probe."),
          ("resid", "probes x samples matrix of within-subject residuals."),
          ("max_gap", "largest gap (bp) across which two probes may be linked."),
          ("rho_min", "minimum residual correlation for a link."),
          ("chunk", "probes per correlation chunk; memory/speed trade-off only.")],
  ret=("integer vector of cluster ids, one per input probe, in input order; "
       "probes in no multi-probe cluster are returned as -1."),
  ex=FIT + '\nresid <- demean_rows_by_group(M, subject)\nchrom <- rep("chr1", n_probes)\npos <- as.integer(seq(1000, 1000 + 200 * n_probes, length.out = n_probes))\ntable(comethylation_clusters(chrom, pos, resid, max_gap = 5000L))'),
 "tv_denoise": dict(
  title="Total-variation denoising of a probe-level effect track",
  desc=("Fits a piecewise-constant effect track by total-variation "
        "regularisation, then reports each segment as the precision-weighted "
        "mean of the input over that segment, so a region effect stays an "
        "unbiased weighted mean in M-value units."),
  params=[("y", "numeric effect estimate per probe."),
          ("w", "precision weight per probe, normally ``1 / se^2``."),
          ("lam", "penalty per gap; length ``length(y) - 1``, as returned by "
           "``distance_penalty()``."),
          ("n_iter", "maximum ADMM iterations."),
          ("rho", "ADMM step size; ``NULL`` picks a scale-free default."),
          ("tol", "convergence tolerance on the primal iterate."),
          ("over_relax", "ADMM over-relaxation factor.")],
  ret="numeric vector of fitted values, one per probe.",
  ex=TRACK + '\nround(theta, 3)'),
 "distance_penalty": dict(
  title="Distance-decaying penalty for the denoising step",
  desc=("Builds the per-gap penalty ``lam0 * scale * exp(-d / decay_bp)``. "
        "When precision weights are supplied ``scale`` is their median, which "
        "makes ``lam0`` dimensionless and the same grid usable on M-values, "
        "beta values or another cohort."),
  params=[("pos", "genomic position per probe, sorted."),
          ("lam0", "dimensionless penalty scale."),
          ("decay_bp", "distance (bp) over which the penalty decays by 1/e."),
          ("w", "optional precision weights, used only to set ``scale``."),
          ("lam_floor", "lower bound on the returned penalty.")],
  ret="numeric vector of penalties, one per gap (``length(pos) - 1``).",
  ex='pos <- as.integer(c(1000, 1100, 1400, 9000))\ndistance_penalty(pos, lam0 = 1, decay_bp = 1000)'),
 "call_regions": dict(
  title="Turn a denoised effect track into regions",
  desc=("Region boundaries are the jumps of the piecewise-constant solution "
        "(or a change of chromosome or cluster), so they are estimated rather "
        "than thresholded off a smooth curve. The region statistic is the "
        "precision-weighted mean effect over its own standard error."),
  params=[("chrom", "chromosome per probe."),
          ("pos", "genomic position per probe."),
          ("beta", "effect estimate per probe."),
          ("se", "standard error per probe."),
          ("theta", "denoised track from ``tv_denoise()``."),
          ("min_probes", "minimum probes per reported region."),
          ("min_effect", "minimum absolute weighted mean effect."),
          ("tol", "tolerance for treating two ``theta`` values as equal."),
          ("cluster", "optional cluster id per probe; a change of id forces a "
           "boundary.")],
  ret=("list with ``table`` (one row per region: ``chr``, ``start``, ``end``, "
       "effect, standard error, z, probe count and 0-based probe index range) "
       "and ``max_abs_z``."),
  ex=TRACK + '\nres <- call_regions(chrom, pos, beta, se, theta, min_probes = 3L)\nres$table[, c("chr", "start", "end")]'),
 "within_subject_permutation": dict(
  title="Permute the exposure within subject",
  desc=("The null that matches a longitudinal design: subject identity, array "
        "type, cohort, chip, visit count and the within-subject covariance are "
        "untouched, and only the pairing between visit and exposure is broken."),
  params=[("exposure", "numeric exposure, one value per sample."),
          ("subject", "subject identifier, one per sample.")],
  ret="numeric vector, the exposure permuted within each subject.",
  ex='set.seed(1)\nwithin_subject_permutation(c(0, 1, 0, 1), c("A", "A", "B", "B"))'),
 "naive_permutation": dict(
  title="Permute the exposure across all samples",
  desc=("Free permutation of the exposure, the scheme bumphunter uses on the "
        "design column. Provided so the mis-calibration of that null on a "
        "longitudinal, multi-array panel can be demonstrated rather than "
        "asserted."),
  params=[("exposure", "numeric exposure, one value per sample."),
          ("subject", "ignored; accepted so the two permutation schemes are "
           "interchangeable as arguments.")],
  ret="numeric vector, the exposure permuted across samples.",
  ex='set.seed(1)\nnaive_permutation(c(0, 1, 0, 1))'),
 "stability_selection": dict(
  title="Subject-level stability selection",
  desc=("Meinshausen-Buhlmann stability selection resampling SUBJECTS, not "
        "samples: resampling samples would split a subject's visit series "
        "across train and test and leak the within-subject effect under test."),
  params=[("fit_fn", "function of a subject subset returning a character "
           "vector of selected region keys."),
          ("subjects", "vector of subject identifiers."),
          ("n_boot", "number of subsamples."),
          ("frac", "fraction of subjects per subsample."),
          ("seed", "RNG seed.")],
  ret="list with ``freq`` (named selection frequencies) and ``n_boot``.",
  ex='subjects <- sprintf("S%02d", 1:12)\nfit_fn <- function(sub) if (length(sub) > 3) c("chr1:100-200") else character(0)\nstability_selection(fit_fn, subjects, n_boot = 20L)$freq'),
 "fit_block_hsmm": dict(
  title="Distance-aware three-state hidden semi-Markov block model",
  desc=("Fits hypo / neutral / hyper states over cluster-level effects with "
        "transition probabilities that decay with genomic distance, so block "
        "boundaries are not limited by a fixed smoothing window. The neutral "
        "state is pinned at mu = 0. Baum-Welch from these fixed starting "
        "values is deterministic."),
  params=[("y", "effect estimate per cluster."),
          ("se", "standard error per cluster."),
          ("pos", "representative genomic position per cluster."),
          ("chrom", "chromosome per cluster."),
          ("length_scale", "distance (bp) over which state persistence decays."),
          ("n_iter", "maximum EM iterations."),
          ("tol", "relative log-likelihood tolerance."),
          ("seed", "accepted and ignored; the fit is deterministic.")],
  ret=("list with ``mu``, ``sigma``, ``pi`` (states ordered by ``mu``), "
       "``length_scale``, ``posterior`` (clusters x 3), ``loglik``, "
       "``n_iter`` and ``neutral``, the 1-based column pinned at mu = 0."),
  ex=HSMM + '\nround(fit$mu, 3)\nfit$neutral'),
 "state_labels": dict(
  title="Direction labels for the three block states",
  desc=("Labels the state columns given which one is pinned at mu = 0. When "
        "all cluster effects fall on one side of zero the pinned state sorts "
        "to an end, two states sit on the same side, and those get a column "
        "suffix rather than being silently collapsed."),
  params=[("neutral", "1-based column of the pinned state, i.e. "
           "``fit_block_hsmm()$neutral``.")],
  ret="character vector of length 3.",
  ex='state_labels(2)\nstate_labels(1)'),
 "call_blocks": dict(
  title="Call blocks from the state posterior",
  desc=("Merges runs of clusters whose posterior for a non-neutral state "
        "exceeds ``min_post``, so each block carries a calibrated confidence "
        "rather than a binary permutation call. Pass ``neutral`` from the fit: "
        "assuming the middle column is neutral mislabels every block when all "
        "cluster effects fall on one side of zero."),
  params=[("chrom", "chromosome per cluster."),
          ("start", "cluster start position."),
          ("end", "cluster end position."),
          ("post", "clusters x 3 posterior matrix from ``fit_block_hsmm()``."),
          ("min_post", "minimum posterior for a cluster to join a block."),
          ("min_clusters", "minimum clusters per reported block."),
          ("neutral", "1-based neutral column, ``fit_block_hsmm()$neutral``.")],
  ret=("data frame with ``chr``, ``start``, ``end``, ``width`` "
       "(inclusive, ``end - start + 1``), "
       "``n_clusters``, ``direction`` and ``posterior``; zero rows if nothing "
       "passes."),
  ex=HSMM + '\ncall_blocks(chrom, pos, pos + 500L, fit$posterior,\n            neutral = fit$neutral, min_clusters = 3L)'),
}


def segments(lines):
    """Return [(name, doc_lines, code_lines)] in source order."""
    starts = [i for i, l in enumerate(lines)
              if re.match(r"^[A-Za-z_.][A-Za-z0-9_.]* *<- *function", l)]
    out = []
    for i in starts:
        name = lines[i].split("<-")[0].strip()
        # documentation block immediately above (roxygen or plain comments)
        j = i - 1
        doc = []
        while j >= 0 and lines[j].lstrip().startswith("#") \
                and not re.match(r"^# -{10,}", lines[j]):
            doc.append(lines[j]); j -= 1
        doc = list(reversed(doc))
        # body: close the signature first, then balance braces. Strings and
        # trailing comments are stripped before counting.
        def clean(s):
            s = re.sub(r'"(\\.|[^"\\])*"', '""', s)
            s = re.sub(r"'(\\.|[^'\\])*'", "''", s)
            return s.split("#")[0]

        par, sig_end = 0, i
        for k in range(i, len(lines)):
            par += clean(lines[k]).count("(") - clean(lines[k]).count(")")
            sig_end = k
            if par <= 0:
                break
        if "{" not in "".join(clean(l) for l in lines[i:sig_end + 1]):
            k = sig_end                        # one-line function body
        else:
            depth = 0
            for k in range(i, len(lines)):
                depth += clean(lines[k]).count("{") - clean(lines[k]).count("}")
                if k >= sig_end and depth <= 0:
                    break
        out.append((name, doc, lines[i:k + 1]))
    return out


def roxygen(name, doc, meta):
    """Build a roxygen block: our tags first, original prose as @details."""
    prose = [re.sub(r"^#'? ?", "", l).rstrip() for l in doc]
    while prose and not prose[0]:
        prose.pop(0)
    while prose and not prose[-1]:
        prose.pop()
    L = ["#' " + meta["title"], "#'"]
    if meta.get("desc"):
        for para in meta["desc"].split("\n\n"):
            L += ["#' " + x for x in textwrap.wrap(para, 76)] + ["#'"]
    if prose:
        L += ["#' @details", "#' Notes carried over from the pipeline source:", "#'"]
        L += ["#' " + p if p else "#'" for p in prose] + ["#'"]
    for pname, ptext in meta.get("params", []):
        w = textwrap.wrap(f"@param {pname} {ptext}", 76,
                          subsequent_indent="  ")
        L += ["#' " + x for x in w]
    if meta.get("ret"):
        L += ["#' " + x for x in textwrap.wrap("@return " + meta["ret"], 76,
                                               subsequent_indent="  ")]
    if name in INTERNAL:
        L += ["#' @keywords internal", "#' @noRd"]
    else:
        if meta.get("ex"):
            L += ["#' @examples"] + ["#' " + x if x else "#'"
                                     for x in meta["ex"].split("\n")]
        L += ["#' @export"]
    return L


DESCRIPTION = """Package: crossarrayEWAS
Type: Package
Title: Design-Aware Region and Block Detection for Cross-Array Longitudinal EWAS
Version: 0.99.0
Authors@R: person("Katarzyna", "Kamieniecka", role = c("aut", "cre"),
    email = "ENTER.YOUR.ADDRESS@example.org")
Description: Estimation and region-calling routines for epigenome-wide
    association studies that combine Illumina HumanMethylation450 and
    MethylationEPIC samples and follow subjects over more than one visit.
    Probe-level effects are fitted within subject, so subject identity, array
    type and chip are absorbed rather than modelled; adjacent probes are
    grouped by within-subject co-methylation; regions are obtained by
    total-variation denoising with a distance-decaying penalty, which places
    boundaries where the effect track jumps instead of thresholding a smooth
    curve; and long-range blocks are called from a three-state hidden
    semi-Markov model whose transition probabilities decay with genomic
    distance. Selection stability and permutation nulls resample subjects, not
    samples, so a subject's visit series is never split.
License: MIT + file LICENSE
Encoding: UTF-8
Depends: R (>= 4.5.0)
Imports: Matrix, jsonlite, limma, stats
Suggests: testthat (>= 3.0.0), knitr, rmarkdown, BiocStyle,
    SummarizedExperiment, GenomicRanges
biocViews: DNAMethylation, DifferentialMethylation, Epigenetics,
    MethylationArray, Regression, Software
BiocType: Software
URL: https://github.com/kkamieniecka/ewas-crossarray-harmonise
BugReports: https://github.com/kkamieniecka/ewas-crossarray-harmonise/issues
Config/testthat/edition: 3
RoxygenNote: 8.1.0
"""

PKGDOC = '''#' crossarrayEWAS: design-aware EWAS region and block detection
#'
#' Numerical core of the cross-array harmonisation pipeline for longitudinal
#' epigenome-wide association studies. The exported functions cover the two
#' stages the pipeline replaces: probe-level within-subject estimation with
#' region calling (\\code{\\link{fit_within}},
#' \\code{\\link{comethylation_clusters}}, \\code{\\link{distance_penalty}},
#' \\code{\\link{tv_denoise}}, \\code{\\link{call_regions}}), and long-range
#' block detection (\\code{\\link{fit_block_hsmm}},
#' \\code{\\link{call_blocks}}).
#'
#' @section Provenance:
#' These routines are an R port of the pipeline's Python reference module, and
#' parts of both were written with AI assistance (Assisted-by: Claude). The
#' port's deliberate divergences from the reference are recorded at each
#' function under \\sQuote{Notes carried over from the pipeline source}.
#'
#' @keywords internal
#' @importFrom Matrix sparseMatrix bandSparse Cholesky solve
#' @importFrom jsonlite fromJSON
#' @importFrom stats median var quantile pnorm
"_PACKAGE"
'''

NAMESPACE_TMPL = """# Generated by build_pkg.py; regenerate with roxygen2::roxygenise().
{exports}
importFrom(Matrix,Cholesky)
importFrom(Matrix,bandSparse)
importFrom(Matrix,solve)
importFrom(Matrix,sparseMatrix)
importFrom(jsonlite,fromJSON)
importFrom(stats,median)
importFrom(stats,pnorm)
importFrom(stats,quantile)
importFrom(stats,var)
"""

NEWS = """# crossarrayEWAS 0.99.0

* First version, extracted from the `ewas-crossarray-harmonise` pipeline
  (`bin/ewasml.R`). Function bodies are unchanged from the pipeline source;
  the package adds documentation, a namespace and unit tests.
* Not yet present, and required before submission: entry points accepting
  `SummarizedExperiment` / `GenomicRatioSet` and returning `GRanges`, a
  `BiocStyle` vignette, and the cross-array harmonisation stage (still
  straight-line script code in `scripts/01_harmonise.R`).
"""

README = """# crossarrayEWAS

Design-aware region and block detection for longitudinal EWAS that combines
Illumina 450K and EPIC samples. This is the numerical core of the
[ewas-crossarray-harmonise](https://github.com/kkamieniecka/ewas-crossarray-harmonise)
pipeline, packaged so the Galaxy wrappers and the Nextflow stages can depend on
a version rather than vendoring the script.

## Installation

Once accepted into Bioconductor:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("crossarrayEWAS")
```

Development version:

```r
BiocManager::install("kkamieniecka/ewas-crossarray-harmonise",
                     subdir = "pkg/crossarrayEWAS")
```

## Scope

`fit_within()` estimates the exposure effect within subject for every probe,
absorbing subject, array type and chip. `comethylation_clusters()` groups
adjacent probes that co-vary within subject. `distance_penalty()` and
`tv_denoise()` produce a piecewise-constant effect track whose jumps are the
region boundaries, which `call_regions()` reports with a precision-weighted
effect and standard error. `fit_block_hsmm()` and `call_blocks()` call
long-range blocks from a distance-aware three-state hidden semi-Markov model.
`stability_selection()` and `within_subject_permutation()` resample subjects,
not samples.

Status: pre-submission skeleton. See `NEWS.md` for what is still missing.
"""

TESTS = {
"linalg": '''test_that("pinv satisfies the pseudo-inverse identity when the design is singular", {
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
''',
"io": '''test_that("read_f64 round-trips a matrix with dimnames", {
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
''',
"variance": '''test_that("moderate_variance shrinks towards the prior", {
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
''',
"folds": '''test_that("fold assignment is deterministic and order-invariant", {
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
''',
"within": '''planted_fit <- function(effect = 1.5, n_probes = 30, n_sub = 8) {
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
''',
"clusters": '''test_that("only genomically adjacent, co-varying probes cluster", {
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
''',
"denoise": '''test_that("the penalty decays with distance and respects its floor", {
  pos <- as.integer(c(1000, 1100, 2000, 20000))
  lam <- distance_penalty(pos, lam0 = 1, decay_bp = 1000)
  expect_length(lam, length(pos) - 1L)
  expect_true(all(diff(lam) < 0))                    # monotone in distance
  expect_true(all(distance_penalty(pos, 1, lam_floor = 0.5) >= 0.5))
})

test_that("precision weights make lam0 dimensionless", {
  pos <- as.integer(seq(1000, 1900, by = 100))
  w <- rep(4, length(pos))
  expect_equal(distance_penalty(pos, 1, w = w),
               4 * distance_penalty(pos, 1), tolerance = 1e-12)
})

test_that("tv_denoise interpolates between the data and a constant", {
  y <- c(rep(1, 6), rep(-1, 6))
  w <- rep(1, 12)
  loose <- tv_denoise(y, w, lam = rep(1e-8, 11), n_iter = 2000L)
  tight <- tv_denoise(y, w, lam = rep(1e6, 11), n_iter = 2000L)
  expect_equal(loose, y, tolerance = 1e-3)
  expect_equal(tight, rep(mean(y), 12), tolerance = 1e-3)
})

test_that("the fit is piecewise constant and segments shrink as lam grows", {
  y <- c(rep(1, 5), rep(0.2, 5), rep(-0.8, 5))
  w <- rep(1, 15)
  n_seg <- function(lam)
    length(unique(round(tv_denoise(y, w, rep(lam, 14), n_iter = 2000L), 6)))
  expect_lte(n_seg(0.5), n_seg(0.01))
  expect_identical(n_seg(1e6), 1L)
})

test_that("each segment is the precision-weighted mean of its probes", {
  y <- c(rep(1, 5), rep(-1, 5))
  w <- c(rep(1, 5), rep(9, 5))
  theta <- tv_denoise(y, w, rep(0.05, 9), n_iter = 3000L)
  for (v in unique(round(theta, 6))) {
    k <- which(round(theta, 6) == v)
    expect_equal(v, sum(w[k] * y[k]) / sum(w[k]), tolerance = 1e-4)
  }
})
''',
"regions": '''track <- function() {
  chrom <- rep("chr1", 12)
  pos <- as.integer(seq(1000, 3200, length.out = 12))
  beta <- c(rep(0.4, 6), rep(-0.1, 6))
  se <- c(rep(0.1, 6), rep(0.05, 6))
  lam <- distance_penalty(pos, lam0 = 0.5, w = 1 / se^2)
  theta <- tv_denoise(beta, 1 / se^2, lam, n_iter = 3000L)
  list(chrom = chrom, pos = pos, beta = beta, se = se, theta = theta)
}

test_that("region boundaries fall where the denoised track jumps", {
  d <- track()
  res <- call_regions(d$chrom, d$pos, d$beta, d$se, d$theta,
                      min_probes = 3L, min_effect = 0.01)
  expect_true(nrow(res$table) >= 1)
  expect_true(all(res$table$start %in% d$pos))
  expect_true(all(res$table$end %in% d$pos))
  expect_true(all(res$table$end >= res$table$start))
})

test_that("the region effect is the precision-weighted mean of its probes", {
  d <- track()
  res <- call_regions(d$chrom, d$pos, d$beta, d$se, d$theta,
                      min_probes = 3L, min_effect = 0.01)
  eff <- grep("^(beta|effect)", names(res$table), value = TRUE)[1]
  for (i in seq_len(nrow(res$table))) {
    k <- (res$table$probe_start[i] + 1L):(res$table$probe_end[i] + 1L)
    w <- 1 / d$se[k]^2
    expect_equal(res$table[[eff]][i], sum(w * d$beta[k]) / sum(w),
                 tolerance = 1e-6)
  }
})

test_that("min_probes and min_effect filter regions out", {
  d <- track()
  expect_identical(
    nrow(call_regions(d$chrom, d$pos, d$beta, d$se, d$theta,
                      min_probes = 99L)$table), 0L)
  expect_identical(
    nrow(call_regions(d$chrom, d$pos, d$beta, d$se, d$theta,
                      min_effect = 10)$table), 0L)
})

test_that("a cluster change forces a boundary", {
  d <- track()
  one <- call_regions(d$chrom, d$pos, d$beta, d$se, d$theta, min_probes = 3L,
                      min_effect = 0.01, cluster = rep(1L, 12))
  two <- call_regions(d$chrom, d$pos, d$beta, d$se, d$theta, min_probes = 3L,
                      min_effect = 0.01, cluster = rep(c(1L, 2L), each = 6))
  expect_gte(nrow(two$table), nrow(one$table))
})
''',
"resample": '''test_that("within-subject permutation preserves each subject's exposures", {
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

test_that("stability selection is reproducible for a given seed", {
  subjects <- sprintf("S%02d", 1:12)
  fit_fn <- function(sub) sub[1]
  expect_identical(stability_selection(fit_fn, subjects, 20L, seed = 3L)$freq,
                   stability_selection(fit_fn, subjects, 20L, seed = 3L)$freq)
})
''',
"blocks": '''block_fit <- function(n = 40) {
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
''',
}


def write_metadata(pkg, exported):
    open(f"{pkg}/DESCRIPTION", "w").write(DESCRIPTION)
    open(f"{pkg}/NAMESPACE", "w").write(NAMESPACE_TMPL.format(
        exports="\n".join(f"export({n})" for n in exported)))
    open(f"{pkg}/NEWS.md", "w").write(NEWS)
    open(f"{pkg}/README.md", "w").write(README)
    # "MIT + file LICENSE" wants a two-field DCF stub; the full text sits
    # beside it as LICENSE.note.
    open(f"{pkg}/LICENSE", "w").write(
        "YEAR: 2026\nCOPYRIGHT HOLDER: Katarzyna Kamieniecka\n")
    shutil.copy("LICENSE", f"{pkg}/LICENSE.note")
    open(f"{pkg}/R/{PKGNAME}-package.R", "w").write(PKGDOC)
    open(f"{pkg}/tests/testthat.R", "w").write(
        f"library(testthat)\nlibrary({PKGNAME})\n\ntest_check(\"{PKGNAME}\")\n")
    for name, body in TESTS.items():
        open(f"{pkg}/tests/testthat/test-{name}.R", "w").write(body)
    open(f"{pkg}/.Rbuildignore", "w").write("^\\.github$\n^build_pkg\\.py$\n")
    open(f"{pkg}/inst/scripts/README.md", "w").write(
        "# Pipeline drivers\n\n"
        "The command-line stage drivers (`04_dmr_ml.R`, `05_blocks_hsmm.R`) "
        "live in the pipeline repository and are not yet installed here. They "
        "call the exported functions plus one internal helper (`pinv()`), "
        "which has to be replaced by `MASS::ginv()` or dropped before they can "
        "run against the installed package.\n")


def main():
    lines = open(SRC).read().split("\n")
    segs = segments(lines)
    assert len(segs) == 22, len(segs)
    assert set(n for n, _, _ in segs) == set(LAYOUT), \
        set(n for n, _, _ in segs) ^ set(LAYOUT)

    if os.path.isdir(PKG):
        shutil.rmtree(PKG)
    for d in ("R", "man", "tests/testthat", "inst/scripts", "vignettes"):
        os.makedirs(f"{PKG}/{d}", exist_ok=True)

    per_file = {f: [] for f in FILE_ORDER}
    for name, doc, code in segs:
        meta = DOCS[name]
        block = roxygen(name, doc, meta) + code + [""]
        per_file[LAYOUT[name]].append("\n".join(block))

    # top-level constants live outside any function body and must come along
    consts = [l for l in lines if re.match(r"^FNV_(OFFSET|PRIME) *<-", l)]
    assert len(consts) == 2, consts
    per_file["folds"].insert(0, "\n".join(consts) + "\n")

    header = ("# Generated by build_pkg.py from "
              "ewas-crossarray-harmonise/scripts/ewasml.R.\n"
              "# Function bodies are verbatim; only documentation was added.\n")
    for f in FILE_ORDER:
        with open(f"{PKG}/R/{f}.R", "w") as fh:
            fh.write(header + "\n" + "\n".join(per_file[f]))

    exported = sorted(n for n in LAYOUT if n not in INTERNAL)
    write_metadata(PKG, exported)
    return exported, {f: len(per_file[f]) for f in FILE_ORDER}


if __name__ == "__main__":
    exp, counts = main()
    print(len(exp), "exported:", " ".join(exp))
    print("files:", counts)
