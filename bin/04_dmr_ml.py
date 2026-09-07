#!/usr/bin/env python3
"""
04_dmr_ml.py -- differentially methylated REGIONS for a cross-array
longitudinal EWAS. Replaces section 2.6 of Aryee et al. 2014 (clusterMaker +
loess bump hunting + design-column permutation).

Pipeline
  1. within-subject probe-level effects (subject intercepts absorb array,
     cohort and chip, which are perfectly confounded in this design)
  2. co-methylation clustering  -> replaces the fixed 300 bp maximum gap
  3. weighted total-variation denoising with a distance-decayed penalty,
     smoothness chosen by held-out-SUBJECT cross-validation
                                 -> replaces loess with a fixed span
  4. region calling at the TV breakpoints, precision-weighted region effects
  5. family-wise error by WITHIN-SUBJECT permutation
                                 -> replaces free permutation of the exposure
  6. stability selection over subject subsamples (a robustness score that a
     single p-value cannot express)
  7. cross-array generalisation: 450K-trained regions scored on the EPIC
     cohort and vice versa, which is the property a harmonisation pipeline
     actually has to demonstrate

Usage
  python 04_dmr_ml.py --in-dir results/harmonised --out-dir results/dmr
"""
from __future__ import annotations

import argparse, json, os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ewasml as E

T0 = time.time()
def log(*a):
    print(f"[{time.time()-T0:7.1f}s]", *a, flush=True)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-dir", required=True,
                   help="directory holding mval.f64/mval_dims.json, pheno_used.csv, probe_annotation.csv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--exposure", default="days_on_clozapine")
    p.add_argument("--exposure-scale", type=float, default=100.0,
                   help="divide the exposure by this so effects read per unit (default per 100 days)")
    p.add_argument("--subject", default="Subject_ID")
    p.add_argument("--array-col", default="Array_Type")
    p.add_argument("--covars", default="smoking_score,cd8t,cd4t,bcell,gran,mono,nk",
                   help="TIME-VARYING covariates only; time-invariant ones are removed by the within transform")
    p.add_argument("--max-gap", type=int, default=1000)
    p.add_argument("--rho-min", type=float, default=0.30)
    p.add_argument("--decay-bp", type=float, default=1000.0)
    p.add_argument("--lam-grid", default="0.05,0.1,0.25,0.5,1.0,2.0,4.0")
    p.add_argument("--cv-folds", type=int, default=5)
    p.add_argument("--min-probes", type=int, default=3)
    p.add_argument("--min-effect", type=float, default=0.05,
                   help="minimum |region effect| in M-value units")
    p.add_argument("--n-perm", type=int, default=200)
    p.add_argument("--n-boot", type=int, default=50)
    p.add_argument("--perm-iter", type=int, default=800,
                   help="ADMM iterations for permutation replicates (the null does not need the observed precision)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--also-naive-perm", action="store_true",
                   help="additionally run the bumphunter-style free permutation, to quantify its mis-calibration")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
def load_inputs(args):
    M, probes, samples = E.read_f64(os.path.join(args.in_dir, "mval"))
    ph = pd.read_csv(os.path.join(args.in_dir, "pheno_used.csv"), dtype={args.subject: str})
    ann = pd.read_csv(os.path.join(args.in_dir, "probe_annotation.csv"))
    if samples is not None:
        key = "Sample_Name" if "Sample_Name" in ph.columns else ph.columns[0]
        ph = ph.set_index(key).loc[samples].reset_index()
    if probes is not None:
        ann = ann.set_index(ann.columns[0]).loc[probes].reset_index()
        ann.columns = ["probe"] + list(ann.columns[1:])
    else:
        ann = ann.rename(columns={ann.columns[0]: "probe"})
    chrcol = next(c for c in ann.columns if c.lower() in ("chr", "seqnames", "chromosome"))
    poscol = next(c for c in ann.columns if c.lower() in ("pos", "mapinfo", "start", "position"))
    ann = ann.rename(columns={chrcol: "chr", poscol: "pos"})
    # genomic order, autosomes only (sex chromosomes need a sex-stratified model)
    ann["chr"] = ann.chr.astype(str).str.replace("^chr", "", regex=True)
    aut = ann.chr.isin([str(i) for i in range(1, 23)])
    ann = ann[aut].copy()
    M = M[np.flatnonzero(aut.values), :]
    ann["chr_num"] = ann.chr.astype(int)
    order = np.lexsort((ann.pos.values, ann.chr_num.values))
    ann = ann.iloc[order].reset_index(drop=True)
    M = M[order, :]
    log(f"loaded {M.shape[0]} autosomal probes x {M.shape[1]} samples")
    return M, ann, ph


def build_design(ph, args):
    expo = ph[args.exposure].astype(float).values / args.exposure_scale
    subj = ph[args.subject].astype(str).values
    arr = ph[args.array_col].astype(str).values
    cov_names = [c for c in args.covars.split(",") if c and c in ph.columns]
    cov = ph[cov_names].astype(float).values if cov_names else None
    return expo, subj, arr, cov, cov_names


def tv_track(ann, beta, se, lam0, decay_bp, n_iter):
    """Chromosome-wise TV denoising of the precision-weighted effect track."""
    w = np.where(se > 0, 1.0 / se ** 2, 0.0)
    theta = np.empty_like(beta)
    cn = ann.chr_num.values
    pos = ann.pos.values.astype(float)
    for c in np.unique(cn):
        idx = np.flatnonzero(cn == c)
        lam = E.distance_penalty(pos[idx], lam0=lam0, decay_bp=decay_bp, w=w[idx])
        theta[idx] = E.tv_denoise(beta[idx], w[idx], lam, n_iter=n_iter)
    return theta, w


def cv_lambda(M, ann, expo, subj, cov, args, rng):
    """
    Choose the smoothness by leaving out whole SUBJECTS.

    Held-out subjects (not held-out samples) are required: a subject's visits
    are the unit the within-subject effect is estimated from, so splitting them
    across train and test would leak the quantity being validated.
    Score = precision-weighted squared error of the training track against the
    held-out subjects' own within-subject effects.
    """
    uniq = np.unique(subj)
    folds = np.array_split(rng.permutation(uniq), args.cv_folds)
    grid = [float(x) for x in args.lam_grid.split(",")]
    rows = []
    for lam0 in grid:
        tot, ntot = 0.0, 0
        for f in folds:
            te = np.isin(subj, f)
            tr = ~te
            if len(np.unique(subj[te])) < 2 or len(np.unique(subj[tr])) < 4:
                continue
            ftr = E.fit_within(M[:, tr], expo[tr], subj[tr],
                               cov[tr] if cov is not None else None)
            fte = E.fit_within(M[:, te], expo[te], subj[te],
                               cov[te] if cov is not None else None)
            th, _ = tv_track(ann, ftr.beta, ftr.se, lam0, args.decay_bp, args.perm_iter)
            wte = np.where(fte.se > 0, 1.0 / fte.se ** 2, 0.0)
            tot += float((wte * (fte.beta - th) ** 2).sum())
            ntot += int((wte > 0).sum())
        rows.append(dict(lam0=lam0, wsse=tot, n=ntot,
                         score=tot / ntot if ntot else np.nan))
        log(f"  CV lam0={lam0:<6g} score={rows[-1]['score']:.6g}")
    cv = pd.DataFrame(rows)
    best = float(cv.loc[cv.score.idxmin(), "lam0"])
    return best, cv


# ---------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    M, ann, ph = load_inputs(args)
    expo, subj, arr, cov, cov_names = build_design(ph, args)
    n_vis = pd.Series(subj).value_counts()
    log(f"{len(np.unique(subj))} subjects, visits/subject "
        f"min={n_vis.min()} median={int(n_vis.median())} max={n_vis.max()}")
    log(f"array split: {pd.Series(arr).value_counts().to_dict()}")
    if (n_vis < 2).any():
        log(f"note: {(n_vis < 2).sum()} subject(s) with a single visit contribute "
            "nothing to a within-subject contrast and are carried but not informative")

    # --- 1. observed probe-level fit ---------------------------------------
    fit = E.fit_within(M, expo, subj, cov)
    log(f"probe fit: df={fit.df}, median |z|={np.median(np.abs(fit.z)):.3f}, "
        f"max |z|={np.abs(fit.z).max():.2f}")

    # --- 2. co-methylation clusters ----------------------------------------
    resid = E._demean_rows_by_group(M, subj)
    Xw = E._demean_by_group(np.column_stack([expo] + ([cov] if cov is not None else [])), subj)
    keep = np.abs(Xw).max(axis=0) > 1e-9
    Xw = Xw[:, keep]
    resid = resid - (resid @ Xw @ np.linalg.pinv(Xw.T @ Xw)) @ Xw.T
    cl = E.comethylation_clusters(ann.chr_num.values, ann.pos.values.astype(float),
                                  resid, max_gap=args.max_gap, rho_min=args.rho_min)
    n_cl = int(cl.max() + 1) if cl.max() >= 0 else 0
    log(f"co-methylation clusters: {n_cl} clusters covering "
        f"{int((cl >= 0).sum())} probes ({(cl >= 0).mean()*100:.1f}%)")
    del resid

    # --- 3. smoothness by held-out-subject CV ------------------------------
    log("cross-validating the TV smoothness over held-out subjects")
    lam0, cv = cv_lambda(M, ann, expo, subj, cov, args, rng)
    cv.to_csv(os.path.join(args.out_dir, "lambda_cv.csv"), index=False)
    log(f"selected lam0={lam0}")

    # --- 4. observed regions ----------------------------------------------
    theta, w = tv_track(ann, fit.beta, fit.se, lam0, args.decay_bp, n_iter=3000)
    obs = E.call_regions(ann.chr.values, ann.pos.values, fit.beta, fit.se, theta,
                         min_probes=args.min_probes, min_effect=args.min_effect,
                         cluster=np.where(cl >= 0, cl, -np.arange(1, len(cl) + 1)))
    reg = pd.DataFrame(obs.table)
    log(f"observed regions: {len(reg)}; max |z| = {obs.max_abs_z:.2f}")
    if reg.empty:
        reg.to_csv(os.path.join(args.out_dir, "dmr_ml.csv"), index=False)
        json.dump(dict(n_regions=0, lam0=lam0), open(os.path.join(args.out_dir, "run_config.json"), "w"), indent=1)
        log("no regions passed the effect-size floor; stopping")
        return 0

    # --- 5. permutation FWER ----------------------------------------------
    def null_max(scheme, n_perm):
        out = []
        fn = (E.within_subject_permutation if scheme == "within" else E.naive_permutation)
        for b in range(n_perm):
            pe = fn(expo, subj, rng)
            f = E.fit_within(M, pe, subj, cov)
            th, _ = tv_track(ann, f.beta, f.se, lam0, args.decay_bp, args.perm_iter)
            r = E.call_regions(ann.chr.values, ann.pos.values, f.beta, f.se, th,
                               min_probes=args.min_probes, min_effect=args.min_effect)
            out.append(r.max_abs_z)
            if (b + 1) % 25 == 0:
                log(f"  {scheme} permutation {b+1}/{n_perm}, running max|z| "
                    f"95th pct = {np.quantile(out, 0.95):.2f}")
        return np.array(out)

    log(f"within-subject permutation null, {args.n_perm} replicates")
    null_w = null_max("within", args.n_perm)
    reg["p_fwer_within"] = [(1 + (null_w >= abs(z)).sum()) / (1 + len(null_w))
                            for z in reg.z]
    nulls = pd.DataFrame({"within": null_w})
    if args.also_naive_perm:
        log(f"naive (bumphunter-style) permutation null, {args.n_perm} replicates")
        null_n = null_max("naive", args.n_perm)
        reg["p_fwer_naive"] = [(1 + (null_n >= abs(z)).sum()) / (1 + len(null_n))
                               for z in reg.z]
        nulls["naive"] = null_n
        log(f"null max|z| 95th pct: within={np.quantile(null_w,0.95):.2f} "
            f"naive={np.quantile(null_n,0.95):.2f}")
    nulls.to_csv(os.path.join(args.out_dir, "permutation_null.csv"), index=False)

    # --- 6. stability selection over subject subsamples --------------------
    log(f"stability selection, {args.n_boot} subject subsamples")
    keys = [(r.chr, r.start, r.end) for r in reg.itertuples()]
    iv = {}
    for c, s, e in keys:
        iv.setdefault(c, []).append((s, e))

    def fit_subset(sub):
        sel = np.isin(subj, sub)
        if len(np.unique(subj[sel])) < 3:
            return []
        f = E.fit_within(M[:, sel], expo[sel], subj[sel],
                         cov[sel] if cov is not None else None)
        th, _ = tv_track(ann, f.beta, f.se, lam0, args.decay_bp, args.perm_iter)
        r = E.call_regions(ann.chr.values, ann.pos.values, f.beta, f.se, th,
                           min_probes=args.min_probes, min_effect=args.min_effect)
        found = []
        for row in r.table:
            for (s, e) in iv.get(row["chr"], []):
                # count as re-discovered when the two intervals overlap
                if row["start"] <= e and s <= row["end"]:
                    found.append((row["chr"], s, e))
        return set(found)

    freq, nb = E.stability_selection(fit_subset, subj, n_boot=args.n_boot,
                                     frac=0.5, seed=args.seed + 7)
    reg["stability"] = [freq.get(k, 0.0) for k in keys]
    log(f"stability: median={np.median(reg.stability):.2f}, "
        f"{(reg.stability >= 0.6).sum()} regions selected in >=60% of subsamples")

    # --- 7. cross-array generalisation ------------------------------------
    log("cross-array generalisation")
    rep = {}
    for a in np.unique(arr):
        sel = arr == a
        if len(np.unique(subj[sel])) < 3:
            continue
        f = E.fit_within(M[:, sel], expo[sel], subj[sel],
                         cov[sel] if cov is not None else None)
        rep[a] = f
        pe, se_ = [], []
        for row in reg.itertuples():
            idx = np.arange(row.probe_start, row.probe_end + 1)
            ws = np.where(f.se[idx] > 0, 1.0 / f.se[idx] ** 2, 0.0)
            pe.append((ws * f.beta[idx]).sum() / ws.sum() if ws.sum() > 0 else np.nan)
            se_.append(np.sqrt(1.0 / ws.sum()) if ws.sum() > 0 else np.nan)
        reg[f"effect_{a}"] = pe
        reg[f"se_{a}"] = se_
    arrs = list(rep)
    if len(arrs) == 2:
        x = reg[f"effect_{arrs[0]}"].values
        y = reg[f"effect_{arrs[1]}"].values
        ok = np.isfinite(x) & np.isfinite(y)
        r_pearson = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 2 else np.nan
        sign_conc = float((np.sign(x[ok]) == np.sign(y[ok])).mean()) if ok.sum() else np.nan
        reg["sign_concordant"] = np.where(np.isfinite(x) & np.isfinite(y),
                                          np.sign(x) == np.sign(y), False)
        log(f"cross-array region effects: r={r_pearson:.3f}, "
            f"sign concordance={sign_conc*100:.1f}% over {int(ok.sum())} regions")
    else:
        r_pearson = sign_conc = np.nan

    reg = reg.sort_values("p_fwer_within").reset_index(drop=True)
    reg.to_csv(os.path.join(args.out_dir, "dmr_ml.csv"), index=False)

    cfg = dict(vars(args))
    cfg.update(n_probes=int(M.shape[0]), n_samples=int(M.shape[1]),
               n_subjects=int(len(np.unique(subj))), covars_used=cov_names,
               lam0_selected=lam0, n_clusters=n_cl, n_regions=int(len(reg)),
               n_fwer_05=int((reg.p_fwer_within <= 0.05).sum()),
               cross_array_r=r_pearson, cross_array_sign_concordance=sign_conc,
               # --n-perm 0 is a legitimate mode (smoothness diagnostics
               # without the expensive null), so the null may be empty
               null_within_q95=(float(np.quantile(null_w, 0.95))
                                if len(null_w) else None),
               probe_fit_df=int(fit.df), seed=args.seed,
               runtime_s=round(time.time() - T0, 1))
    if args.also_naive_perm and len(nulls.get("naive", [])):
        cfg["null_naive_q95"] = float(np.quantile(nulls["naive"].values, 0.95))
    json.dump(cfg, open(os.path.join(args.out_dir, "run_config.json"), "w"),
              indent=1, default=str)
    log(f"wrote {len(reg)} regions; {int((reg.p_fwer_within <= 0.05).sum())} at FWER<=0.05")
    return 0


if __name__ == "__main__":
    sys.exit(main())
