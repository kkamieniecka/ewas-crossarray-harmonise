#!/usr/bin/env python
"""
Generate fixtures and reference outputs for tests/test_equivalence.R.

The fixtures are seeded pseudo-random arrays, not methylation data: this test
asks whether bin/ewasml.R computes the same numbers as bin/ewasml.py on
identical input, which is a question about the port and not about any cohort.
Scientific agreement on real data is a separate check -- see
docs/design-rationale.md, "Two implementations of one estimator".

Everything the R side needs is written as plain CSV plus the same
column-major .f64 + JSON header pair that stage 01 writes, so the R test also
exercises read_f64().

Usage:  python tests/gen_equivalence_fixtures.py --out-dir <dir>
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import ewasml as E  # noqa: E402

N_PROBES = 4000
N_SUBJECTS = 12
SEED = 20240607


def write_f64(stem, mat, rownames=None, colnames=None):
    mat = np.asarray(mat, dtype="<f8")
    # tofile() always writes C order, so ravel in F order to get the
    # column-major dump that stage 01 writes and read_f64 expects
    mat.ravel(order="F").tofile(f"{stem}.f64")
    hdr = dict(nrow=int(mat.shape[0]), ncol=int(mat.shape[1]), order="F")
    if rownames is not None:
        hdr["rownames"] = list(rownames)
    if colnames is not None:
        hdr["colnames"] = list(colnames)
    with open(f"{stem}_dims.json", "w") as fh:
        json.dump(hdr, fh)


def build_design(rng):
    """Unbalanced visit counts, as in a real longitudinal series."""
    visits = rng.integers(1, 6, size=N_SUBJECTS)
    subject = np.repeat([f"S{i:02d}" for i in range(N_SUBJECTS)], visits)
    n = len(subject)
    exposure = np.round(rng.uniform(0, 900, size=n) / 100.0, 4)
    covars = np.column_stack([
        rng.normal(0.2, 0.05, n),      # a cell fraction
        rng.normal(-1.0, 0.4, n),      # a smoking score
        rng.normal(0.1, 0.03, n),      # another cell fraction
    ]).round(6)
    array = np.where(np.arange(n) % 3 == 0, "450K", "EPIC")
    return subject, exposure, covars, array


def build_annotation(rng):
    """Two chromosomes, clustered spacing so both dense and isolated probes occur."""
    per = N_PROBES // 2
    rows = []
    for chrom in ("1", "2"):
        pos = 0
        for i in range(per):
            # alternate tight CpG-island-like runs with large gaps
            step = int(rng.integers(20, 220)) if i % 25 else int(rng.integers(5000, 60000))
            pos += step
            rows.append((chrom, pos))
    ann = pd.DataFrame(rows, columns=["chr", "pos"])
    ann["chr_num"] = ann.chr.astype(int)
    ann["probe"] = [f"cg{i:08d}" for i in range(len(ann))]
    order = np.lexsort((ann.pos.values, ann.chr_num.values))
    return ann.iloc[order].reset_index(drop=True)


def build_matrix(rng, ann, subject, exposure, covars):
    """M-value-like matrix with subject intercepts, a real effect, and AR noise."""
    n = len(subject)
    sub_idx = pd.factorize(subject)[0]
    base = rng.normal(0, 1.5, size=(N_PROBES, 1))
    sub_shift = rng.normal(0, 0.8, size=(1, N_SUBJECTS))[:, sub_idx]
    eff = np.zeros(N_PROBES)
    # a handful of contiguous runs carry a genuine exposure effect
    for start in (100, 600, 1200, 2500, 3300):
        eff[start:start + 12] = rng.normal(0.6, 0.1)
    noise = rng.normal(0, 0.45, size=(N_PROBES, n))
    # spatial correlation along the genome, so clustering has something to find
    noise = np.cumsum(noise, axis=0) * 0.08 + noise
    M = base + sub_shift + np.outer(eff, exposure) + covars[:, 0] * 0.3 + noise
    return np.ascontiguousarray(M)


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", required=True)
    args = p.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    O = lambda f: os.path.join(args.out_dir, f)  # noqa: E731

    rng = np.random.default_rng(SEED)
    subject, exposure, covars, array = build_design(rng)
    ann = build_annotation(rng)
    M = build_matrix(rng, ann, subject, exposure, covars)

    # ---- inputs -----------------------------------------------------------
    write_f64(O("mval"), M, rownames=list(ann.probe), colnames=None)
    ann.to_csv(O("annotation.csv"), index=False)
    pd.DataFrame(dict(subject=subject, exposure=exposure, array=array,
                      cov1=covars[:, 0], cov2=covars[:, 1], cov3=covars[:, 2])
                 ).to_csv(O("design.csv"), index=False)

    # ---- 1. probe-level fit, with and without variance moderation ---------
    ref = {}
    for tag, shrink in (("shrunk", True), ("raw", False)):
        f = E.fit_within(M, exposure, subject, covars, shrink_var=shrink)
        pd.DataFrame(dict(beta=f.beta, se=f.se, z=f.z)).to_csv(
            O(f"fit_{tag}.csv"), index=False)
        ref[f"fit_{tag}_df"] = int(f.df)
    fit = E.fit_within(M, exposure, subject, covars)

    # ---- 2. variance moderation in isolation ------------------------------
    #     an independent s2 vector, so the test is not conditional on the fit
    s2 = np.exp(rng.normal(-1.0, 0.7, size=3000)) ** 2
    df_mod = 25
    pd.DataFrame(dict(s2=s2,
                      post_mom=E._moderate_variance(s2, df_mod))).to_csv(
        O("moderation.csv"), index=False)
    ref["moderation_df"] = df_mod

    # ---- 3. co-methylation clusters ---------------------------------------
    resid = E._demean_rows_by_group(M, subject)
    Xw = E._demean_by_group(np.column_stack([exposure, covars]), subject)
    Xw = Xw[:, np.abs(Xw).max(axis=0) > 1e-9]
    resid = resid - (resid @ Xw @ np.linalg.pinv(Xw.T @ Xw)) @ Xw.T
    write_f64(O("resid"), resid)
    cl = E.comethylation_clusters(ann.chr_num.values, ann.pos.values.astype(float),
                                  resid, max_gap=1000, rho_min=0.30)
    pd.DataFrame(dict(cluster=cl)).to_csv(O("clusters.csv"), index=False)

    # ---- 4. penalty + TV track --------------------------------------------
    w = np.where(fit.se > 0, 1.0 / fit.se ** 2, 0.0)
    lam0, decay_bp = 0.5, 1000.0
    theta = np.empty_like(fit.beta)
    lam_all = []
    for c in np.unique(ann.chr_num.values):
        idx = np.flatnonzero(ann.chr_num.values == c)
        lam = E.distance_penalty(ann.pos.values[idx].astype(float), lam0=lam0,
                                 decay_bp=decay_bp, w=w[idx])
        lam_all.append(pd.DataFrame(dict(chr_num=c, i=idx[:-1], lam=lam)))
        theta[idx] = E.tv_denoise(fit.beta[idx], w[idx], lam, n_iter=3000)
    pd.concat(lam_all).to_csv(O("penalty.csv"), index=False)
    pd.DataFrame(dict(theta=theta)).to_csv(O("tv_track.csv"), index=False)
    ref.update(lam0=lam0, decay_bp=decay_bp)

    # ---- 5. region calling ------------------------------------------------
    cluster_arg = np.where(cl >= 0, cl, -np.arange(1, len(cl) + 1))
    reg = E.call_regions(ann.chr.values, ann.pos.values, fit.beta, fit.se, theta,
                         min_probes=3, min_effect=0.05, cluster=cluster_arg)
    pd.DataFrame(reg.table).to_csv(O("regions.csv"), index=False)
    ref["max_abs_z"] = float(reg.max_abs_z)

    # ---- 6. permutation draws --------------------------------------------
    # The two languages have different generators, so R cannot reproduce these
    # draws. The test consumes them instead, which makes the downstream FWER
    # arithmetic exactly comparable.
    prng = np.random.default_rng(SEED + 1)
    draws = {}
    for scheme, fn in (("within", E.within_subject_permutation),
                       ("naive", E.naive_permutation)):
        d = np.column_stack([fn(exposure, subject, prng) for _ in range(5)])
        draws[scheme] = d
        pd.DataFrame(d, columns=[f"rep{i}" for i in range(d.shape[1])]).to_csv(
            O(f"perm_{scheme}.csv"), index=False)
    # max|z| of the null replicate, the quantity the FWER p-value counts
    null_max = []
    for i in range(draws["within"].shape[1]):
        f = E.fit_within(M, draws["within"][:, i], subject, covars)
        th = np.empty_like(f.beta)
        for c in np.unique(ann.chr_num.values):
            idx = np.flatnonzero(ann.chr_num.values == c)
            wl = np.where(f.se > 0, 1.0 / f.se ** 2, 0.0)[idx]
            lam = E.distance_penalty(ann.pos.values[idx].astype(float),
                                     lam0=lam0, decay_bp=decay_bp, w=wl)
            th[idx] = E.tv_denoise(f.beta[idx], wl, lam, n_iter=800)
        r = E.call_regions(ann.chr.values, ann.pos.values, f.beta, f.se, th,
                           min_probes=3, min_effect=0.05)
        null_max.append(r.max_abs_z)
    pd.DataFrame(dict(max_abs_z=null_max)).to_csv(O("null_within.csv"), index=False)

    # ---- 7. identification guard and deterministic fold split ------------
    sub_u = np.unique(subject)
    ref["within_df_full"] = E.within_df(exposure, subject, covars)
    ref["within_df_half"] = E.within_df(
        exposure[np.isin(subject, sub_u[:6])], subject[np.isin(subject, sub_u[:6])],
        covars[np.isin(subject, sub_u[:6])])
    rows = []
    for k in (3, 5, 7):
        for f_i, f in enumerate(E.fold_assign(subject, k, seed=1)):
            rows.extend(dict(k=k, fold=f_i, subject=s) for s in f)
    pd.DataFrame(rows).to_csv(O("folds.csv"), index=False)

    with open(O("reference.json"), "w") as fh:
        json.dump(ref, fh, indent=1)
    print(f"fixtures written to {args.out_dir}: "
          f"{N_PROBES} probes x {len(subject)} samples, "
          f"{len(pd.DataFrame(reg.table))} reference regions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
