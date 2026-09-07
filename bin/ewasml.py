"""
ewasml -- region and block detection for cross-array longitudinal EWAS.

Replaces the two region-level algorithms of the minfi/EWASGalaxy stack
(Aryee et al., Bioinformatics 2014;30:1363, sections 2.6 and 2.7).

What was wrong with them for THIS design
----------------------------------------
Section 2.6 (bump hunting) does three things: (i) `clusterMaker` groups probes
by a fixed maximum genomic gap, (ii) a probe-level regression coefficient is
smoothed with loess inside each cluster, (iii) significance comes from
permuting the design matrix column of interest.

  (i) fails under cross-array harmonisation. The 450K/EPIC intersection has a
      median inter-probe gap of 335 bp, so the 300 bp default splits ~52% of
      neighbouring probe pairs, and the clustering a study gets depends on
      which arrays happened to be in it. Clusters are a property of the array
      panel, not of the biology.
  (ii) loess with a fixed span treats a 20 bp step and a 20 kb step as equal
      neighbours, and gives no breakpoints -- the region boundary is whatever
      the arbitrary cutoff on the smoothed curve happens to cross.
  (iii) is invalid for repeated measures. Permuting the exposure column across
      all samples destroys the within-subject pairing, so the permutation null
      has more independent information than the data. In this design it is
      additionally broken because subject, cohort, chip and array type are
      nested: a free permutation mixes samples across arrays and manufactures
      between-array variance that no real reassignment of exposure could.

Section 2.7 (block finding) collapses open-sea probes with a fixed 500 bp gap /
1500 bp width rule and then loess-smooths with a >=250 kb window. Under
harmonisation the open-sea pool is exactly the content that differs most
between arrays (78% of EPIC-only probes are open sea), so fixed-width collapse
produces non-comparable units, and the paper itself notes that block
*boundaries* are resolution-limited.

What this module does instead
-----------------------------
`fit_within`          within-subject (fixed-effects) probe-level estimator.
                      Subject intercepts absorb array, cohort and chip, which
                      are perfectly confounded here, so the exposure effect is
                      identified from within-subject contrasts only.
`comethylation_clusters`
                      replaces clusterMaker: probes are joined when they are
                      both close AND empirically co-methylated in these data,
                      so clusters are data-defined and array-panel-independent.
`tv_denoise`          replaces loess: weighted total-variation (1D fused-lasso)
                      denoising, solved by ADMM. Convex, precision-weighted,
                      with a distance-decayed penalty so a long gap is free to
                      break. Produces explicit piecewise-constant segments, ie
                      real region boundaries, and one tunable smoothness
                      parameter chosen by held-out-subject cross-validation.
`within_subject_permutation`
                      replaces design-column permutation: the exposure is
                      shuffled among the visits OF THE SAME SUBJECT, which
                      preserves subject, array, cohort, chip and the
                      within-subject covariance and destroys only the
                      association with time.
`stability_selection` Meinshausen-Buhlmann subject-level subsampling, giving a
                      selection frequency per region as a robustness score
                      that a single permutation p-value cannot express.
`fit_block_hsmm`      replaces cpgCollapse + wide-window loess: a
                      distance-aware 3-state hidden Markov model over
                      open-sea cluster effects, where the transition matrix is
                      an explicit function of genomic distance and emissions
                      carry each cluster's own standard error. Posterior
                      decoding gives soft, probabilistic block boundaries.

All estimators are deterministic given `seed`.
"""

from __future__ import annotations

import json
import numpy as np
from dataclasses import dataclass, field
from scipy.linalg import solve_banded
from scipy import sparse
from scipy.stats import norm


# ----------------------------------------------------------------------------
# I/O for the matrix written by 01_harmonise.R
# ----------------------------------------------------------------------------
def read_f64(stem: str):
    """Read the column-major float64 dump + JSON header written by R."""
    with open(f"{stem}_dims.json") as fh:
        d = json.load(fh)
    x = np.fromfile(f"{stem}.f64", dtype="<f8")
    if x.size != d["nrow"] * d["ncol"]:
        raise ValueError(f"{stem}.f64 has {x.size} values, header says "
                         f"{d['nrow']}x{d['ncol']}")
    mat = x.reshape((d["nrow"], d["ncol"]), order="F")
    rn = d.get("rownames"); cn = d.get("colnames")
    rn = list(rn) if rn is not None else None
    cn = list(cn) if cn is not None else None
    return mat, rn, cn


# ----------------------------------------------------------------------------
# probe-level within-subject estimator
# ----------------------------------------------------------------------------
@dataclass
class WithinFit:
    beta: np.ndarray      # effect per probe, M-value units per exposure unit
    se: np.ndarray        # standard error
    z: np.ndarray         # beta / se
    df: int               # residual degrees of freedom
    n_subjects: int
    n_samples: int


def _demean_by_group(X: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Subtract the group mean from every column (within transform)."""
    out = np.array(X, dtype=float, copy=True)
    uniq, inv = np.unique(groups, return_inverse=True)
    counts = np.bincount(inv).astype(float)
    for j in range(out.shape[1]):
        sums = np.bincount(inv, weights=out[:, j])
        out[:, j] -= (sums / counts)[inv]
    return out


def fit_within(M: np.ndarray, exposure: np.ndarray, subject: np.ndarray,
               covars: np.ndarray | None = None,
               shrink_var: bool = True) -> WithinFit:
    """
    Within-subject (fixed-effects) estimate of the exposure effect per probe.

    M         probes x samples matrix of M-values
    exposure  length-samples exposure (eg days on clozapine / 100)
    subject   length-samples subject identifier
    covars    samples x k time-varying covariates (cell fractions, smoking).
              Time-INVARIANT covariates (age, sex, array type, cohort) must NOT
              be passed: the within transform removes them exactly, and
              including them makes the design rank-deficient.

    Subject intercepts are eliminated by demeaning rather than estimated, so
    cost is one matrix product regardless of the number of subjects -- which is
    what makes the permutation and subsampling loops below affordable.
    """
    M = np.asarray(M, dtype=float)
    X = exposure.reshape(-1, 1).astype(float)
    if covars is not None and covars.size:
        X = np.hstack([X, np.asarray(covars, dtype=float)])
    Xw = _demean_by_group(X, subject)
    # drop covariate columns that the within transform annihilated
    keep = np.abs(Xw).max(axis=0) > 1e-9
    if not keep[0]:
        raise ValueError("exposure has no within-subject variation; a "
                         "longitudinal contrast is not identified")
    Xw = Xw[:, keep]
    Mw = _demean_rows_by_group(M, subject)

    XtX = Xw.T @ Xw
    XtXi = np.linalg.pinv(XtX)
    B = Mw @ Xw @ XtXi                      # probes x p
    resid = Mw - B @ Xw.T
    n_sub = len(np.unique(subject))
    df = M.shape[1] - n_sub - Xw.shape[1]
    if df <= 0:
        raise ValueError("no residual degrees of freedom")
    s2 = (resid ** 2).sum(axis=1) / df
    if shrink_var:
        # limma-style variance moderation: shrink each probe's residual
        # variance towards the global trend. Without it, region statistics are
        # dominated by probes that happen to have a tiny sample variance.
        s2 = _moderate_variance(s2, df)
    v_beta = XtXi[0, 0]
    se = np.sqrt(s2 * v_beta)
    beta = B[:, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(se > 0, beta / se, 0.0)
    return WithinFit(beta=beta, se=se, z=z, df=df,
                     n_subjects=n_sub, n_samples=M.shape[1])


def _demean_rows_by_group(M: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Row-wise (per probe) subtraction of the group mean across samples."""
    uniq, inv = np.unique(groups, return_inverse=True)
    counts = np.bincount(inv).astype(float)
    S = np.zeros((M.shape[0], len(uniq)))
    np.add.at(S.T, inv, M.T)
    return M - (S / counts)[:, inv]


def _moderate_variance(s2: np.ndarray, df: int) -> np.ndarray:
    """
    Empirical-Bayes shrinkage of per-probe residual variance towards a prior
    (the limma model): s2_post = (d0*s0^2 + df*s2) / (d0 + df), with (d0, s0^2)
    from the method-of-moments fit on log s2.
    """
    ok = np.isfinite(s2) & (s2 > 0)
    if ok.sum() < 100:
        return s2
    from scipy.special import digamma, polygamma
    e = np.log(s2[ok]) - digamma(df / 2) + np.log(df / 2)
    ebar = e.mean()
    target = e.var(ddof=1) - polygamma(1, df / 2)
    d0 = np.inf if target <= 0 else 2.0 * _trigamma_inv(target)
    if not np.isfinite(d0):
        s0sq = np.exp(ebar + digamma(df / 2) - np.log(df / 2))
        return np.full_like(s2, s0sq)
    s0sq = np.exp(ebar + digamma(d0 / 2) - np.log(d0 / 2))
    out = np.array(s2, copy=True)
    out[ok] = (d0 * s0sq + df * s2[ok]) / (d0 + df)
    return out


def _trigamma_inv(x: float) -> float:
    """Solve polygamma(1, y) = x for y (Newton, as in limma's trigammaInverse)."""
    from scipy.special import polygamma
    y = 0.5 + 1.0 / x
    for _ in range(60):
        f = polygamma(1, y) - x
        d = polygamma(2, y)
        if d == 0:
            break
        step = f / d
        y_new = y - step
        if y_new <= 0:
            y_new = y / 2.0
        if abs(y_new - y) < 1e-10 * max(1.0, y):
            y = y_new
            break
        y = y_new
    return y


# ----------------------------------------------------------------------------
# replaces clusterMaker: data-driven co-methylation clustering
# ----------------------------------------------------------------------------
def comethylation_clusters(chrom: np.ndarray, pos: np.ndarray,
                           resid: np.ndarray, max_gap: int = 1000,
                           rho_min: float = 0.30) -> np.ndarray:
    """
    Cluster probes that are both physically close and empirically
    co-methylated, returning an integer cluster id per probe (-1 = singleton).

    `resid` is the within-subject residual matrix (probes x samples); the
    correlation is therefore between-visit co-variation, not between-subject
    or between-array variation, which is what a region is supposed to mean.

    Unlike a fixed maximum-gap rule, the same locus yields the same cluster
    whether it was measured on 450K, EPIC or the intersection, because the
    criterion is a property of the data at that locus.
    """
    order = np.lexsort((pos, chrom))
    n = len(pos)
    parent = np.arange(n)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    # normalise residual rows once so correlation is a dot product
    R = resid - resid.mean(axis=1, keepdims=True)
    nrm = np.sqrt((R ** 2).sum(axis=1))
    nrm[nrm == 0] = np.inf
    R = R / nrm[:, None]

    for k in range(len(order) - 1):
        i, j = order[k], order[k + 1]
        if chrom[i] != chrom[j]:
            continue
        if pos[j] - pos[i] > max_gap:
            continue
        if float(R[i] @ R[j]) >= rho_min:
            union(i, j)

    roots = np.array([find(i) for i in range(n)])
    uniq, inv, counts = np.unique(roots, return_inverse=True, return_counts=True)
    cid = np.where(counts[inv] >= 2, inv, -1)
    # renumber the surviving clusters consecutively
    keep = np.unique(cid[cid >= 0])
    remap = {old: new for new, old in enumerate(keep)}
    return np.array([remap.get(c, -1) for c in cid], dtype=int)


# ----------------------------------------------------------------------------
# replaces loess smoothing: weighted total-variation denoising
# ----------------------------------------------------------------------------
def tv_denoise(y: np.ndarray, w: np.ndarray, lam: np.ndarray,
               n_iter: int = 3000, rho: float | None = None,
               tol: float = 1e-6, over_relax: float = 1.7) -> np.ndarray:
    """
    Solve  min_theta  0.5*sum_i w_i (y_i - theta_i)^2 + sum_i lam_i |d theta_i|
    by ADMM. y must already be in genomic order within one contiguous run.

    w    per-probe precision weight (1/se^2)
    lam  length n-1 penalty on each successive difference. Passing a
         distance-decayed lam is what makes this a *genomic* smoother: a large
         gap gets lam ~ 0 and is free to break, which is precisely the
         behaviour loess with a fixed span cannot express.

    The returned track is DEBIASED: the l1 penalty is used to choose where the
    breakpoints are, and each resulting segment is then reported as the
    precision-weighted mean of y over that segment. This is the relaxed-lasso
    convention, and it matters for interpretation -- the reported region effect
    is an unbiased weighted mean in M-value units, not a shrunken one, so it
    can be compared directly with a per-probe coefficient.
    """
    n = len(y)
    if n == 1:
        return y.astype(float).copy()
    y = np.asarray(y, float); w = np.asarray(w, float)
    lam = np.asarray(lam, float)
    if rho is None:
        # scale the augmented-Lagrangian parameter to the data precision, so
        # convergence does not depend on the units of y
        pos_w = w[w > 0]
        rho = float(np.median(pos_w)) if pos_w.size else 1.0
        rho = max(rho, 1e-8)

    theta = y.astype(float).copy()
    z = np.zeros(n - 1)
    u = np.zeros(n - 1)

    # (W + rho*D'D) is symmetric tridiagonal; build once in banded form.
    diag_DtD = np.concatenate([[1.0], np.full(max(n - 2, 0), 2.0), [1.0]])
    ab = np.zeros((3, n))
    ab[0, 1:] = -rho
    ab[1, :] = w + rho * diag_DtD
    ab[2, :-1] = -rho

    Wy = w * y
    eps_abs = tol * np.sqrt(n)
    for it in range(n_iter):
        v = z - u
        rhs = Wy.copy()
        rhs[:-1] -= rho * v
        rhs[1:] += rho * v
        theta = solve_banded((1, 1), ab, rhs, check_finite=False)
        dtheta = np.diff(theta)
        # over-relaxation accelerates the fused-lasso ADMM appreciably
        dhat = over_relax * dtheta + (1.0 - over_relax) * z
        z_old = z
        z = _soft(dhat + u, lam / rho)
        u = u + dhat - z
        r_prim = np.linalg.norm(dtheta - z)
        r_dual = rho * np.linalg.norm(z - z_old)
        if r_prim < eps_abs + tol * np.linalg.norm(z) and \
           r_dual < eps_abs + tol * rho * np.linalg.norm(u):
            break
    # The ADMM primal iterate theta is only approximately piecewise constant,
    # but the split variable z IS exactly sparse (it comes straight out of a
    # soft threshold), so z is what defines the breakpoints. Using |diff(theta)|
    # against a tolerance instead would make the segmentation depend on the
    # scale of y.
    keep_jump = z != 0.0
    seg = np.concatenate([[0], np.cumsum(keep_jump)])
    out = np.empty_like(theta)
    for s in range(seg[-1] + 1):
        idx = seg == s
        ws = w[idx]
        out[idx] = (ws * y[idx]).sum() / ws.sum() if ws.sum() > 0 else theta[idx].mean()
    return out


def _soft(x: np.ndarray, t: np.ndarray) -> np.ndarray:
    return np.sign(x) * np.maximum(np.abs(x) - t, 0.0)


def distance_penalty(pos: np.ndarray, lam0: float, decay_bp: float = 1000.0,
                     w: np.ndarray | None = None,
                     lam_floor: float = 0.0) -> np.ndarray:
    """
    lam_i = lam0 * scale * exp(-d_i / decay_bp) for successive distances d_i.

    `scale` is median(w) when the precision weights are supplied, which makes
    lam0 DIMENSIONLESS: lam0 ~ 1 means "one probe's worth of precision resists
    one unit of jump", so the same lam0 grid is usable on M-values, beta
    values, or a different cohort. Without this scaling lam0 has to be retuned
    whenever the noise level changes, which is the practical reason a fixed
    loess span does not transfer between studies.

    decay_bp sets how fast neighbouring probes stop being neighbours. At
    d = decay_bp the penalty is 1/e of its maximum; by 5*decay_bp it is
    effectively zero and the track is free to break.
    """
    d = np.diff(np.asarray(pos, dtype=float))
    scale = 1.0
    if w is not None:
        pos_w = np.asarray(w, float)
        pos_w = pos_w[pos_w > 0]
        if pos_w.size:
            scale = float(np.median(pos_w))
    return np.maximum(lam0 * scale * np.exp(-d / decay_bp), lam_floor)


# ----------------------------------------------------------------------------
# region calling from a piecewise-constant track
# ----------------------------------------------------------------------------
@dataclass
class Regions:
    table: list = field(default_factory=list)   # list of dicts
    max_abs_z: float = 0.0


def call_regions(chrom, pos, beta, se, theta, min_probes: int = 3,
                 min_effect: float = 0.05, tol: float = 1e-6,
                 cluster: np.ndarray | None = None) -> Regions:
    """
    Turn the denoised track into regions. Segment boundaries are where the
    piecewise-constant solution jumps (or the chromosome/cluster changes), so
    boundaries are estimated, not thresholded off a smooth curve.

    Region statistic is the precision-weighted mean effect divided by its
    standard error, which accounts for the fact that probes in a region carry
    very different amounts of information.
    """
    n = len(beta)
    w = np.where(se > 0, 1.0 / se ** 2, 0.0)
    brk = np.zeros(n, dtype=bool)
    brk[0] = True
    same_chr = chrom[1:] == chrom[:-1]
    jump = np.abs(np.diff(theta)) > tol
    same_cl = np.ones(n - 1, dtype=bool) if cluster is None \
        else (cluster[1:] == cluster[:-1])
    brk[1:] = (~same_chr) | jump | (~same_cl)
    seg = np.cumsum(brk) - 1

    out = []
    max_abs_z = 0.0
    for s in range(seg[-1] + 1):
        idx = np.flatnonzero(seg == s)
        if len(idx) < min_probes:
            continue
        ws = w[idx]
        if ws.sum() <= 0:
            continue
        bbar = float((ws * beta[idx]).sum() / ws.sum())
        se_bar = float(np.sqrt(1.0 / ws.sum()))
        zbar = bbar / se_bar if se_bar > 0 else 0.0
        if abs(bbar) < min_effect:
            continue
        max_abs_z = max(max_abs_z, abs(zbar))
        out.append(dict(chr=str(chrom[idx[0]]), start=int(pos[idx[0]]),
                        end=int(pos[idx[-1]]),
                        width=int(pos[idx[-1]] - pos[idx[0]] + 1),
                        n_probes=int(len(idx)),
                        effect_M=bbar, se=se_bar, z=zbar,
                        p_pointwise=float(2 * norm.sf(abs(zbar))),
                        theta=float(theta[idx[0]]),
                        probe_start=int(idx[0]), probe_end=int(idx[-1])))
    return Regions(table=out, max_abs_z=max_abs_z)


# ----------------------------------------------------------------------------
# inference: the two permutation schemes
# ----------------------------------------------------------------------------
def within_subject_permutation(exposure: np.ndarray, subject: np.ndarray,
                               rng: np.random.Generator) -> np.ndarray:
    """
    Shuffle the exposure among the visits of the same subject.

    This is the null that matches a longitudinal design: subject identity,
    array type, cohort, chip, number of visits and the within-subject
    covariance are all untouched, and only the pairing between visit and
    exposure is broken. A subject with a single usable visit contributes
    nothing under this null, exactly as it contributes nothing to the
    within-subject estimator.
    """
    out = np.array(exposure, dtype=float, copy=True)
    for s in np.unique(subject):
        idx = np.flatnonzero(subject == s)
        if len(idx) > 1:
            out[idx] = rng.permutation(out[idx])
    return out


def naive_permutation(exposure: np.ndarray, subject: np.ndarray,
                      rng: np.random.Generator) -> np.ndarray:
    """
    Free permutation of the exposure across all samples -- the scheme
    bumphunter uses on the design column. Retained only so the pipeline can
    demonstrate that it is mis-calibrated here (it breaks the subject pairing
    and reassigns exposures across array types).
    """
    return rng.permutation(np.asarray(exposure, dtype=float))


# ----------------------------------------------------------------------------
# stability selection
# ----------------------------------------------------------------------------
def stability_selection(fit_fn, subjects: np.ndarray, n_boot: int = 100,
                        frac: float = 0.5, seed: int = 1) -> tuple:
    """
    Meinshausen-Buhlmann stability selection at the SUBJECT level.

    fit_fn(subject_subset) must return an iterable of region keys selected on
    that subsample. Returns (selection_frequency dict, n_boot).
    Subsampling subjects (not samples) is required: resampling samples would
    split a subject's visit series across train and test and leak the
    within-subject effect being tested.
    """
    rng = np.random.default_rng(seed)
    uniq = np.unique(subjects)
    k = max(2, int(np.floor(frac * len(uniq))))
    counts: dict = {}
    for b in range(n_boot):
        sub = rng.choice(uniq, size=k, replace=False)
        for key in fit_fn(sub):
            counts[key] = counts.get(key, 0) + 1
    return {k_: v / n_boot for k_, v in counts.items()}, n_boot


# ----------------------------------------------------------------------------
# replaces cpgCollapse + wide-window loess: distance-aware block HSMM
# ----------------------------------------------------------------------------
@dataclass
class HSMMFit:
    mu: np.ndarray
    sigma: np.ndarray
    pi: np.ndarray
    length_scale: float
    posterior: np.ndarray      # n x 3
    loglik: float
    n_iter: int


def _dist_transitions(d: np.ndarray, pi: np.ndarray, L: float) -> np.ndarray:
    """
    A(d) = p(d) I + (1 - p(d)) 1 pi',  p(d) = exp(-d / L).

    Adjacent clusters 100 bp apart are almost certainly in the same state;
    clusters 5 Mb apart are effectively independent draws from the stationary
    distribution. cpgCollapse + a fixed 250 kb loess window cannot represent
    that, which is why its block boundaries are resolution-limited.
    """
    p = np.exp(-d / L)[:, None, None]
    I = np.eye(3)[None, :, :]
    P = np.broadcast_to(pi, (len(d), 3)).copy()[:, None, :]
    return p * I + (1.0 - p) * P


def fit_block_hsmm(y: np.ndarray, se: np.ndarray, pos: np.ndarray,
                   chrom: np.ndarray, length_scale: float = 250_000.0,
                   n_iter: int = 60, tol: float = 1e-5,
                   seed: int = 1) -> HSMMFit:
    """
    Three-state (hypomethylated / neutral / hypermethylated) hidden Markov
    model over collapsed open-sea cluster effects, with distance-dependent
    transitions and heteroscedastic emissions y_k ~ N(mu_s, sigma_s^2 + se_k^2).

    Fitted by Baum-Welch. Posterior state probabilities give soft block
    boundaries and a per-block posterior that replaces the permutation FWER of
    the loess-based block finder.
    """
    y = np.asarray(y, float); se = np.asarray(se, float)
    pos = np.asarray(pos, float)
    rng = np.random.default_rng(seed)
    q = np.nanquantile(y, [0.05, 0.5, 0.95])
    mu = np.array([min(q[0], -1e-3), 0.0, max(q[2], 1e-3)])
    sd0 = max(np.nanstd(y), 1e-3)
    sigma = np.array([sd0, sd0 / 2, sd0])
    pi = np.array([0.05, 0.90, 0.05])

    # contiguous runs = chromosomes
    bounds = [0] + list(np.flatnonzero(chrom[1:] != chrom[:-1]) + 1) + [len(y)]
    runs = [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)]

    prev_ll = -np.inf
    post = np.zeros((len(y), 3))
    for it in range(n_iter):
        ll = 0.0
        num_mu = np.zeros(3); den_mu = np.zeros(3)
        num_s2 = np.zeros(3); den_s2 = np.zeros(3)
        pi_acc = np.zeros(3)
        for a, b in runs:
            yy, ss, pp = y[a:b], se[a:b], pos[a:b]
            n = len(yy)
            if n == 0:
                continue
            var = sigma[None, :] ** 2 + (ss ** 2)[:, None]
            B = np.exp(-0.5 * (yy[:, None] - mu[None, :]) ** 2 / var) / np.sqrt(2 * np.pi * var)
            B = np.maximum(B, 1e-300)
            if n == 1:
                g = pi * B[0]; ll += np.log(g.sum()); post[a:b] = g / g.sum()
                pi_acc += post[a:b].sum(axis=0)
                num_mu += (post[a:b] * yy[:, None]).sum(axis=0)
                den_mu += post[a:b].sum(axis=0)
                continue
            A = _dist_transitions(np.diff(pp), pi, length_scale)
            alpha = np.zeros((n, 3)); c = np.zeros(n)
            alpha[0] = pi * B[0]; c[0] = alpha[0].sum(); alpha[0] /= c[0]
            for t in range(1, n):
                alpha[t] = (alpha[t - 1] @ A[t - 1]) * B[t]
                c[t] = alpha[t].sum(); alpha[t] /= c[t]
            beta = np.zeros((n, 3)); beta[-1] = 1.0
            for t in range(n - 2, -1, -1):
                beta[t] = A[t] @ (B[t + 1] * beta[t + 1]) / c[t + 1]
            g = alpha * beta
            g /= g.sum(axis=1, keepdims=True)
            post[a:b] = g
            ll += np.log(c).sum()
            pi_acc += g.sum(axis=0)
            num_mu += (g * yy[:, None]).sum(axis=0)
            den_mu += g.sum(axis=0)
        # M step: mu, then sigma with the measurement-error part held out
        mu = np.where(den_mu > 0, num_mu / np.maximum(den_mu, 1e-12), mu)
        mu[1] = 0.0                      # neutral state is pinned at no change
        for a, b in runs:
            yy, ss = y[a:b], se[a:b]
            g = post[a:b]
            num_s2 += (g * ((yy[:, None] - mu[None, :]) ** 2 - (ss ** 2)[:, None])).sum(axis=0)
            den_s2 += g.sum(axis=0)
        sigma = np.sqrt(np.maximum(num_s2 / np.maximum(den_s2, 1e-12), 1e-8))
        pi = pi_acc / pi_acc.sum()
        pi = np.maximum(pi, 1e-6); pi /= pi.sum()
        if abs(ll - prev_ll) < tol * max(1.0, abs(prev_ll)):
            prev_ll = ll
            break
        prev_ll = ll
    order = np.argsort(mu)
    return HSMMFit(mu=mu[order], sigma=sigma[order], pi=pi[order],
                   length_scale=length_scale, posterior=post[:, order],
                   loglik=float(prev_ll), n_iter=it + 1)


def call_blocks(chrom, start, end, post, min_post: float = 0.80,
                min_clusters: int = 3):
    """
    Maximal runs whose posterior for a non-neutral state exceeds `min_post`.
    The reported block posterior is the mean over its clusters, so a block is
    reported with a calibrated confidence instead of a binary permutation call.
    """
    state = np.where(post.max(axis=1) >= min_post, post.argmax(axis=1), 1)
    out = []
    i = 0
    n = len(state)
    while i < n:
        s = state[i]
        if s == 1:
            i += 1
            continue
        j = i
        while (j + 1 < n and state[j + 1] == s and chrom[j + 1] == chrom[i]):
            j += 1
        if (j - i + 1) >= min_clusters:
            out.append(dict(chr=str(chrom[i]), start=int(start[i]),
                            end=int(end[j]),
                            width=int(end[j] - start[i] + 1),
                            n_clusters=int(j - i + 1),
                            direction="hypo" if s == 0 else "hyper",
                            posterior=float(post[i:j + 1, s].mean())))
        i = j + 1
    return out
