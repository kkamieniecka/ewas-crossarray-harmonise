#!/usr/bin/env python3
"""
05_blocks_hsmm.py -- large-scale methylation BLOCKS for a cross-array
longitudinal EWAS. Replaces section 2.7 of Aryee et al. 2014 (cpgCollapse with
a fixed 500 bp gap / 1500 bp width rule, followed by bump hunting with a
>=250 kb loess window).

Two things make the original unusable after 450K/EPIC harmonisation:

  * It operates on open-sea probes, and open sea is exactly the content that
    differs most between the arrays -- 78% of EPIC-only probes are open sea
    versus 36% of shared probes. A fixed-width collapse therefore produces
    units whose composition depends on which array a sample came from.
  * The block boundary is wherever a >=250 kb loess curve crosses a cutoff.
    The paper itself notes this makes boundaries resolution-limited, and it
    gives no measure of how confident any particular boundary is.

Replacement
  * clusters are formed where open-sea probes are close AND empirically
    co-methylated, so a cluster is a property of the locus, not of the panel
  * cluster-level effects are estimated by the within-subject model, carrying
    each cluster's own standard error
  * a 3-state hidden Markov model (hypo / neutral / hyper) runs along the
    chromosome with a transition matrix that is an explicit function of
    genomic distance, A(d) = e^{-d/L} I + (1 - e^{-d/L}) 1 pi'
  * posterior decoding gives a soft, calibrated block boundary and a per-block
    posterior probability instead of a thresholded smooth curve

Usage
  python 05_blocks_hsmm.py --in-dir results/harmonised \\
      --probe-map results/crossarray_probe_map.csv.gz --out-dir results/blocks
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
    p.add_argument("--in-dir", required=True)
    p.add_argument("--probe-map", required=True,
                   help="crossarray_probe_map.csv.gz (probe, status, island context)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--exposure", default="days_on_clozapine")
    p.add_argument("--exposure-scale", type=float, default=100.0)
    p.add_argument("--subject", default="Subject_ID")
    p.add_argument("--array-col", default="Array_Type")
    p.add_argument("--covars", default="smoking_score,cd8t,cd4t,bcell,gran,mono,nk")
    p.add_argument("--max-gap", type=int, default=1500,
                   help="maximum gap when forming open-sea clusters (cpgCollapse used 500 bp gap / 1500 bp width)")
    p.add_argument("--rho-min", type=float, default=0.20,
                   help="minimum within-subject co-methylation to join two open-sea probes")
    p.add_argument("--length-scale", type=float, default=250_000.0,
                   help="HSMM distance length scale L in bp; the 250 kb default matches the loess window it replaces")
    p.add_argument("--min-post", type=float, default=0.80)
    p.add_argument("--min-clusters", type=int, default=3)
    p.add_argument("--min-fit-clusters", type=int, default=200,
                   help="refuse to fit the block model on fewer open-sea "
                        "clusters than this; the default is a production "
                        "floor, lower it only for small test panels")
    p.add_argument("--fixed-collapse", action="store_true",
                   help="also report the legacy fixed-width collapse, for comparison of the resulting units")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    M, probes, samples = E.read_f64(os.path.join(args.in_dir, "mval"))
    ph = pd.read_csv(os.path.join(args.in_dir, "pheno_used.csv"), dtype={args.subject: str})
    ann = pd.read_csv(os.path.join(args.in_dir, "probe_annotation.csv"))
    if samples is not None:
        key = "Sample_Name" if "Sample_Name" in ph.columns else ph.columns[0]
        ph = ph.set_index(key).loc[samples].reset_index()
    ann = ann.rename(columns={ann.columns[0]: "probe"})
    if probes is not None:
        ann = ann.set_index("probe").loc[probes].reset_index()
    chrcol = next(c for c in ann.columns if c.lower() in ("chr", "seqnames", "chromosome"))
    poscol = next(c for c in ann.columns if c.lower() in ("pos", "mapinfo", "start", "position"))
    ann = ann.rename(columns={chrcol: "chr", poscol: "pos"})
    ann["chr"] = ann.chr.astype(str).str.replace("^chr", "", regex=True)

    pm = pd.read_csv(args.probe_map, usecols=["probe", "status", "Relation_to_UCSC_CpG_Island"])
    ann = ann.merge(pm, on="probe", how="left")

    sel = (ann.Relation_to_UCSC_CpG_Island.fillna("OpenSea") == "OpenSea") \
        & ann.chr.isin([str(i) for i in range(1, 23)]) \
        & (ann.status == "shared_450K_EPIC")
    log(f"open-sea, shared, autosomal probes: {int(sel.sum())} of {len(ann)} harmonised")
    if sel.sum() < 1000:
        raise SystemExit("too few open-sea probes survived filtering to look for blocks")
    ann = ann[sel.values].copy()
    M = M[np.flatnonzero(sel.values), :]
    ann["chr_num"] = ann.chr.astype(int)
    order = np.lexsort((ann.pos.values, ann.chr_num.values))
    ann = ann.iloc[order].reset_index(drop=True)
    M = M[order, :]

    expo = ph[args.exposure].astype(float).values / args.exposure_scale
    subj = ph[args.subject].astype(str).values
    arr = ph[args.array_col].astype(str).values
    cov_names = [c for c in args.covars.split(",") if c and c in ph.columns]
    cov = ph[cov_names].astype(float).values if cov_names else None

    # --- data-driven open-sea clusters -------------------------------------
    resid = E._demean_rows_by_group(M, subj)
    cl = E.comethylation_clusters(ann.chr_num.values, ann.pos.values.astype(float),
                                  resid, max_gap=args.max_gap, rho_min=args.rho_min)
    del resid
    n_cl = int(cl.max() + 1) if cl.max() >= 0 else 0
    log(f"co-methylated open-sea clusters: {n_cl} covering {int((cl>=0).sum())} probes")
    # Relaxing the clustering merges probes further and yields FEWER clusters,
    # so the way out of this is to split them, not to relax them.
    if n_cl < args.min_fit_clusters:
        raise SystemExit(f"too few open-sea clusters ({n_cl} < "
                         f"{args.min_fit_clusters}) to fit a block model; "
                         "raise --rho-min or lower --max-gap to split "
                         "clusters, supply a larger panel, or lower "
                         "--min-fit-clusters if a coarse fit is intended")

    # collapse each cluster to its mean M-value profile (the cpgCollapse step,
    # but over data-defined clusters)
    keepc = cl >= 0
    idx = np.flatnonzero(keepc)
    cid = cl[idx]
    order2 = np.argsort(cid, kind="stable")
    idx, cid = idx[order2], cid[order2]
    bounds = np.searchsorted(cid, np.arange(n_cl + 1))
    Mc = np.empty((n_cl, M.shape[1]))
    cchr = np.empty(n_cl, dtype=object); cstart = np.empty(n_cl, dtype=np.int64)
    cend = np.empty(n_cl, dtype=np.int64); cnp = np.empty(n_cl, dtype=int)
    for k in range(n_cl):
        rows = idx[bounds[k]:bounds[k + 1]]
        Mc[k] = M[rows].mean(axis=0)
        cchr[k] = ann.chr.values[rows[0]]
        cstart[k] = ann.pos.values[rows].min()
        cend[k] = ann.pos.values[rows].max()
        cnp[k] = len(rows)
    cmid = ((cstart + cend) / 2.0)
    o = np.lexsort((cmid, np.array([int(c) for c in cchr])))
    Mc, cchr, cstart, cend, cnp, cmid = Mc[o], cchr[o], cstart[o], cend[o], cnp[o], cmid[o]
    log(f"cluster widths: median {int(np.median(cend-cstart+1))} bp, "
        f"median {int(np.median(cnp))} probes/cluster")

    # --- cluster-level within-subject effects ------------------------------
    fit = E.fit_within(Mc, expo, subj, cov)
    log(f"cluster effects: median |z|={np.median(np.abs(fit.z)):.3f}, "
        f"max |z|={np.abs(fit.z).max():.2f}")

    # --- distance-aware block HSMM -----------------------------------------
    hs = E.fit_block_hsmm(fit.beta, fit.se, cmid, cchr,
                          length_scale=args.length_scale, seed=args.seed)
    log(f"HSMM converged in {hs.n_iter} iterations, loglik={hs.loglik:.1f}")
    names = E.state_labels(hs.neutral)
    log(f"  state means (M-value per {args.exposure_scale:g} units): "
        + ", ".join(f"{names[k]}={hs.mu[k]:.4f}" for k in range(3)))
    log(f"  state sd: {np.round(hs.sigma,4).tolist()}; stationary pi: {np.round(hs.pi,4).tolist()}")
    if hs.neutral != 1:
        # Every cluster effect fell on one side of zero, so the state pinned at
        # mu = 0 sorted to an end. Direction is labelled relative to that
        # column, not to the middle one: with neutral at column 0 no block can
        # be hypomethylated, and the middle column is a real effect state that
        # would go unreported if it were treated as no change.
        log(f"  NOTE: the zero-effect state is column {hs.neutral}, not the "
            f"middle one -- all cluster effects lie on one side of zero, so "
            f"only {'hyper' if hs.neutral == 0 else 'hypo'}methylated blocks "
            f"can be called")

    blocks = E.call_blocks(cchr, cstart, cend, hs.posterior,
                           min_post=args.min_post, min_clusters=args.min_clusters,
                           neutral=hs.neutral)
    bl = pd.DataFrame(blocks)
    log(f"blocks called: {len(bl)}")

    clus = pd.DataFrame(dict(chr=cchr, start=cstart, end=cend, n_probes=cnp,
                             effect_M=fit.beta, se=fit.se, z=fit.z))
    for k in range(3):
        clus[f"post_{names[k]}"] = hs.posterior[:, k]

    # --- cross-array check on the called blocks ----------------------------
    per_arr = {}
    for a in np.unique(arr):
        s = arr == a
        # Subject count alone does not make the within fit identified: after
        # differencing out subject means, the exposure and the time-varying
        # covariates must still leave residual degrees of freedom. Same guard
        # as stage 04; without it an unidentified arm returns numbers rather
        # than declining to fit.
        df_a = E.within_df(expo[s], subj[s], cov[s] if cov is not None else None)
        if len(np.unique(subj[s])) < 3 or df_a <= 0:
            log(f"  {a}: no within-subject information to fit (df={df_a}); skipped")
            continue
        per_arr[a] = E.fit_within(Mc[:, s], expo[s], subj[s],
                                  cov[s] if cov is not None else None)
        clus[f"effect_{a}"] = per_arr[a].beta
    arrs = list(per_arr)
    r_pearson = sign_conc = np.nan
    if len(arrs) == 2 and len(bl):
        ka, kb = arrs
        ea, eb = [], []
        for row in bl.itertuples():
            m = (clus.chr.values == row.chr) & (clus.start.values >= row.start) \
                & (clus.end.values <= row.end)
            wa = 1.0 / per_arr[ka].se[m] ** 2
            wb = 1.0 / per_arr[kb].se[m] ** 2
            ea.append((wa * per_arr[ka].beta[m]).sum() / wa.sum())
            eb.append((wb * per_arr[kb].beta[m]).sum() / wb.sum())
        bl[f"effect_{ka}"] = ea
        bl[f"effect_{kb}"] = eb
        ok = np.isfinite(ea) & np.isfinite(eb)
        if ok.sum() > 2:
            r_pearson = float(np.corrcoef(np.array(ea)[ok], np.array(eb)[ok])[0, 1])
        sign_conc = float((np.sign(ea) == np.sign(eb)).mean()) if len(ea) else np.nan
        log(f"cross-array block effects: r={r_pearson:.3f}, "
            f"sign concordance={sign_conc*100:.1f}%")

    # --- optional legacy fixed-width collapse, for unit comparison ---------
    legacy = None
    if args.fixed_collapse:
        cn = ann.chr_num.values
        pos = ann.pos.values.astype(float)
        brk = np.ones(len(pos), dtype=bool)
        brk[1:] = (cn[1:] != cn[:-1]) | (np.diff(pos) > 500)
        gid = np.cumsum(brk) - 1
        widths = pd.DataFrame(dict(g=gid, pos=pos)).groupby("g").pos.agg(["min", "max", "count"])
        widths["width"] = widths["max"] - widths["min"] + 1
        widths = widths[(widths["count"] >= 2) & (widths.width <= 1500)]
        # median of an empty selection is NaN, and json.dump would write a bare
        # NaN token that no strict JSON parser accepts; report null instead.
        legacy = dict(n_units=int(len(widths)),
                      median_width_bp=float(widths.width.median()) if len(widths) else None,
                      median_probes=float(widths["count"].median()) if len(widths) else None)
        mw = legacy["median_width_bp"]
        log(f"legacy fixed-width collapse would give {legacy['n_units']} units "
            f"(median {'NA' if mw is None else format(mw, '.0f')} bp) vs {n_cl} "
            f"co-methylation clusters (median {np.median(cend-cstart+1):.0f} bp)")

    bl.to_csv(os.path.join(args.out_dir, "blocks_hsmm.csv"), index=False)
    clus.to_csv(os.path.join(args.out_dir, "openSea_cluster_effects.csv.gz"),
                index=False, compression="gzip")
    def jsonable(x):
        """NaN/Inf are not JSON; a skipped array or an empty selection must
        serialise as null so any strict parser can read the file."""
        if isinstance(x, dict):
            return {k: jsonable(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [jsonable(v) for v in x]
        if isinstance(x, float) and not np.isfinite(x):
            return None
        return x

    json.dump(jsonable(dict(mu=hs.mu.tolist(), sigma=hs.sigma.tolist(), pi=hs.pi.tolist(),
                   length_scale=hs.length_scale, loglik=hs.loglik,
                   n_iter=hs.n_iter, neutral_state=int(hs.neutral),
                   state_labels=list(names), implementation="python",
                   n_clusters=n_cl,
                   n_openSea_probes=int(len(ann)), n_blocks=int(len(bl)),
                   cross_array_r=r_pearson, cross_array_sign_concordance=sign_conc,
                   legacy_fixed_collapse=legacy, covars_used=cov_names,
                   args=vars(args), runtime_s=round(time.time() - T0, 1))),
              open(os.path.join(args.out_dir, "hsmm_params.json"), "w"),
              indent=1, default=str)
    log("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
