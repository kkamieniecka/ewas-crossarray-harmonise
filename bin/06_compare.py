#!/usr/bin/env python3
"""
06_compare.py -- benchmark the replacement region/block finders against the
legacy minfi implementations they replace, on the identical harmonised matrix.

The comparison is deliberately NOT "which method calls more regions". On a
cross-array longitudinal design more calls is the failure mode, because the
legacy model puts between-subject (= between-array, between-cohort,
between-chip) variation in the residual and permutes the exposure freely. The
metrics reported are therefore:

  n_regions            how many were called at all
  n_fwer05             how many survive family-wise correction
  median_width_bp      the unit size the method actually resolves
  median_n_probes      probes supporting a call
  cross_array_r        correlation of the region effect estimated in the 450K
                       cohort alone with the same region estimated in the EPIC
                       cohort alone -- the property a harmonisation pipeline
                       has to demonstrate
  sign_concordance     fraction of regions with the same direction in both
  null_q95             95th percentile of the permutation null of max|z|;
                       an inflated null is the signature of a permutation
                       scheme that manufactures variance the design cannot

Usage
  python 06_compare.py --dmr-dir results/dmr --blocks-dir results/blocks \\
      --baseline-dir results/baseline --out-dir results/comparison
"""
from __future__ import annotations

import argparse, json, os, sys
import numpy as np
import pandas as pd


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dmr-dir", required=True)
    p.add_argument("--blocks-dir", required=True)
    p.add_argument("--baseline-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--probe-model-dir", default=None,
                   help="02_probe_model.R output; enables the common-footing "
                        "cross-array columns (shared within-subject estimator)")
    p.add_argument("--no-figure", action="store_true")
    return p.parse_args(argv)


def jload(path):
    return json.load(open(path)) if os.path.exists(path) else {}


def cload(path, **kw):
    return pd.read_csv(path, **kw) if os.path.exists(path) else pd.DataFrame()


def _num(s):
    return pd.to_numeric(s, errors="coerce") if s is not None else pd.Series(dtype=float)


def summarise(name, level, tab, width_cols, probe_col, fwer_col,
              eff_cols, null_q95=np.nan, cfg=None):
    cfg = cfg or {}
    row = dict(method=name, level=level, n_regions=int(len(tab)))
    if len(tab):
        st, en = width_cols
        if st in tab and en in tab:
            w = _num(tab[en]) - _num(tab[st]) + 1
            row["median_width_bp"] = float(np.nanmedian(w))
            row["total_bp"] = float(np.nansum(w))
        if probe_col in tab:
            row["median_n_probes"] = float(np.nanmedian(_num(tab[probe_col])))
        if fwer_col and fwer_col in tab:
            row["n_fwer05"] = int((_num(tab[fwer_col]) <= 0.05).sum())
        if eff_cols and all(c in tab for c in eff_cols):
            a, b = _num(tab[eff_cols[0]]).values, _num(tab[eff_cols[1]]).values
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() > 2:
                row["cross_array_r"] = float(np.corrcoef(a[ok], b[ok])[0, 1])
            if ok.sum():
                row["sign_concordance"] = float((np.sign(a[ok]) == np.sign(b[ok])).mean())
                row["n_cross_array"] = int(ok.sum())
    row["null_q95_max_abs_z"] = null_q95
    for k in ("runtime_s", "n_clusters"):
        if k in cfg:
            row[k] = cfg[k]
    return row


def common_footing(tab, per_array, chr_col="chr"):
    """Cross-array agreement for one region set under a SHARED estimator.

    Each method reports its own per-array effects, but those come from
    different per-array models (the legacy baseline fits no subject term, the
    replacement fits a within-subject model), so the two methods'
    cross_array_r values are not comparable: the difference conflates which
    units were selected with how the per-array effect was estimated.

    This re-estimates every method's per-array region effect from the SAME
    within-subject per-probe fits (02_probe_model.R), as an inverse-variance
    weighted mean over the probes falling inside the region interval. The only
    remaining difference between methods is then the unit set each one chose,
    which is the thing being compared.

    per_array: {array_name: DataFrame(probe, chr, pos, beta_M, se_M)}
    """
    if not len(tab) or len(per_array) < 2:
        return {}
    names = list(per_array)[:2]
    idx = {}
    for a in names:
        d = per_array[a].dropna(subset=["beta_M", "se_M"])
        d = d[np.isfinite(d.se_M) & (d.se_M > 0)]
        byc = {}
        for c, g in d.groupby("chr", sort=False):
            g = g.sort_values("pos")
            byc[str(c)] = (g.pos.values.astype(float),
                           g.beta_M.values.astype(float),
                           g.se_M.values.astype(float))
        idx[a] = byc

    def region_effect(byc, c, s, e):
        g = byc.get(c)
        if g is None:
            return np.nan
        pos, eff, se = g
        lo = np.searchsorted(pos, s, side="left")
        hi = np.searchsorted(pos, e, side="right")
        if hi <= lo:
            return np.nan
        w = 1.0 / se[lo:hi] ** 2
        return float(np.sum(w * eff[lo:hi]) / np.sum(w))

    def as_chr(v):
        # iterrows() upcasts a mixed numeric frame, so an integer chromosome
        # arrives as 5.0 and would build the label "chr5.0"
        s = str(v)
        if s.startswith("chr"):
            return s
        try:
            return "chr" + str(int(float(s)))
        except ValueError:
            return "chr" + s

    out = {a: [] for a in names}
    for _, r in tab.iterrows():
        c = as_chr(r[chr_col])
        s, e = float(r["start"]), float(r["end"])
        for a in names:
            out[a].append(region_effect(idx[a], c, s, e))
    a = np.asarray(out[names[0]], float)
    b = np.asarray(out[names[1]], float)
    ok = np.isfinite(a) & np.isfinite(b)
    res = {"n_cross_array_common": int(ok.sum())}
    if ok.sum() > 2:
        res["cross_array_r_common"] = float(np.corrcoef(a[ok], b[ok])[0, 1])
        res["sign_concordance_common"] = float(
            (np.sign(a[ok]) == np.sign(b[ok])).mean())
    return res


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    dmr_cfg = jload(os.path.join(args.dmr_dir, "run_config.json"))
    blk_cfg = jload(os.path.join(args.blocks_dir, "hsmm_params.json"))
    base_cfg = jload(os.path.join(args.baseline_dir, "baseline_summary.json"))

    dmr = cload(os.path.join(args.dmr_dir, "dmr_ml.csv"))
    nulls = cload(os.path.join(args.dmr_dir, "permutation_null.csv"))
    blk = cload(os.path.join(args.blocks_dir, "blocks_hsmm.csv"))
    bh = cload(os.path.join(args.baseline_dir, "bumphunter_regions.csv"))
    bf = cload(os.path.join(args.baseline_dir, "blockfinder_blocks.csv"))

    def arr_cols(tab):
        c = [x for x in tab.columns if x.startswith("effect_")
             and x not in ("effect_M",)]
        return c[:2] if len(c) >= 2 else None

    q95_within = float(np.quantile(nulls["within"], 0.95)) if "within" in nulls else np.nan
    q95_naive = float(np.quantile(nulls["naive"], 0.95)) if "naive" in nulls else np.nan

    rows = [
        summarise("bumphunter (minfi 2.6)", "region", bh,
                  ("start", "end"), "L", "fwer", arr_cols(bh),
                  null_q95=np.nan, cfg=base_cfg),
        summarise("TV + within-subject perm (2.6 replacement)", "region", dmr,
                  ("start", "end"), "n_probes", "p_fwer_within", arr_cols(dmr),
                  null_q95=q95_within, cfg=dmr_cfg),
        summarise("blockFinder (minfi 2.7)", "block", bf,
                  ("start", "end"), "L", "fwer", arr_cols(bf),
                  null_q95=np.nan, cfg=base_cfg),
        summarise("distance-aware HSMM (2.7 replacement)", "block", blk,
                  ("start", "end"), "n_clusters", None, arr_cols(blk),
                  null_q95=np.nan, cfg=blk_cfg),
    ]
    cmp = pd.DataFrame(rows)

    # ---- common-footing cross-array agreement ----------------------------
    if args.probe_model_dir:
        per_array = {}
        for a in ("450K", "EPIC"):
            f = os.path.join(args.probe_model_dir, f"probe_stats_{a}.csv")
            d = cload(f, usecols=["probe", "chr", "pos", "beta_M", "se_M"])
            if len(d):
                per_array[a] = d
        if len(per_array) == 2:
            for i, tab in ((0, bh), (1, dmr), (2, bf), (3, blk)):
                for k, v in common_footing(tab, per_array).items():
                    cmp.loc[i, k] = v
    # a naive-permutation column on the SAME observed regions isolates the
    # permutation scheme from every other change
    if "p_fwer_naive" in dmr.columns:
        cmp.loc[cmp.method.str.contains("replacement") & (cmp.level == "region"),
                "n_fwer05_naive_perm"] = int((_num(dmr.p_fwer_naive) <= 0.05).sum())
        cmp.loc[cmp.method.str.contains("replacement") & (cmp.level == "region"),
                "null_q95_naive_perm"] = q95_naive

    front = ["method", "level", "n_regions", "n_fwer05", "median_width_bp",
             "median_n_probes", "cross_array_r_common",
             "sign_concordance_common", "cross_array_r", "sign_concordance",
             "null_q95_max_abs_z"]
    cmp = cmp[[c for c in front if c in cmp.columns]
              + [c for c in cmp.columns if c not in front]]
    cmp.to_csv(os.path.join(args.out_dir, "method_comparison.csv"), index=False)
    print(cmp.to_string(index=False))

    # ---- short written report -------------------------------------------
    lines = ["# Region and block detection: legacy vs replacement", "",
             f"Harmonised matrix: {dmr_cfg.get('n_probes','?')} autosomal probes x "
             f"{dmr_cfg.get('n_samples','?')} samples, "
             f"{dmr_cfg.get('n_subjects','?')} subjects.", "",
             cmp.to_markdown(index=False), ""]
    if np.isfinite(q95_within) and np.isfinite(q95_naive):
        lines += [
            "## Permutation scheme",
            "",
            f"Null 95th percentile of max|z|: within-subject **{q95_within:.2f}**, "
            f"free (bumphunter-style) **{q95_naive:.2f}**. Both nulls were computed "
            "from the identical observed statistic, so the difference is attributable "
            "to the permutation scheme alone. Because subject, cohort, chip and array "
            "type are nested here, a free permutation of the exposure creates "
            "between-array contrasts that no real reassignment of visit times could "
            "produce; the resulting null is therefore not a null for this design.", ""]
    if "cross_array_r" in cmp.columns:
        lines += ["## Cross-array generalisation", "",
                  "Both columns correlate a region's effect estimated in the 450K "
                  "cohort alone against the same region estimated in the EPIC cohort "
                  "alone -- the only metric here that a method cannot improve by "
                  "simply calling more regions.", "",
                  "`cross_array_r` uses each method's own per-array estimator and is "
                  "therefore **not comparable between methods**: the legacy per-array "
                  "fit has no subject term while the replacement's does, so the "
                  "difference mixes unit selection with estimator choice.", "",
                  "`cross_array_r_common` re-estimates every method's per-array "
                  "region effect from the same within-subject per-probe fits "
                  "(inverse-variance weighted over the probes inside each interval). "
                  "The only remaining difference is which units each method chose, "
                  "so this is the column to read when comparing methods.", "",
                  "Read `n_cross_array_common` alongside it. A correlation over "
                  "the handful of units a block method returns carries almost no "
                  "information, and should not be compared against one computed "
                  "over two hundred regions.", ""]
    open(os.path.join(args.out_dir, "comparison.md"), "w").write("\n".join(lines))

    # ---- figure ----------------------------------------------------------
    if not args.no_figure and len(cmp):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        try:
            from figure_style import apply_figure_style  # optional
            apply_figure_style()
        except Exception:
            pass
        fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6))
        short = {"bumphunter (minfi 2.6)": "bumphunter\n(2.6)",
                 "TV + within-subject perm (2.6 replacement)": "TV + within\nperm (new)",
                 "blockFinder (minfi 2.7)": "blockFinder\n(2.7)",
                 "distance-aware HSMM (2.7 replacement)": "HSMM\n(new)"}
        lbl = [short.get(m, m) for m in cmp.method]
        col = ["#8c8c8c" if "new" not in l else "#1f6fb4" for l in lbl]
        for ax, (c, ttl) in zip(axes, [
                ("n_regions", "regions / blocks called"),
                ("cross_array_r_common" if "cross_array_r_common" in cmp
                 else "cross_array_r",
                 "cross-array effect correlation\n(shared within-subject estimator)"),
                ("median_width_bp", "median unit width (bp)")]):
            v = _num(cmp.get(c, pd.Series([np.nan] * len(cmp)))).values
            ax.bar(range(len(cmp)), np.nan_to_num(v), color=col)
            ax.set_xticks(range(len(cmp))); ax.set_xticklabels(lbl, fontsize=7)
            ax.set_title(ttl, fontsize=9)
            if c == "median_width_bp":
                ax.set_yscale("log")
            for i, x in enumerate(v):
                if np.isfinite(x):
                    ax.text(i, x, f"{x:,.2f}" if abs(x) < 10 else f"{x:,.0f}",
                            ha="center", va="bottom", fontsize=7)
            if c.startswith("cross_array_r"):
                # a correlation over 4 blocks is not the same evidence as one
                # over 200 regions, so the unit count belongs on the panel
                nn = _num(cmp.get("n_cross_array_common",
                                  cmp.get("n_cross_array"))).values
                ax.set_xticklabels(
                    [f"{l}\nn={int(k)}" if np.isfinite(k) else l
                     for l, k in zip(lbl, nn)], fontsize=7)
                ax.axhline(0, color="k", lw=0.6)
        fig.tight_layout()
        fig.savefig(os.path.join(args.out_dir, "fig_method_comparison.png"), dpi=200)
    return 0


if __name__ == "__main__":
    sys.exit(main())
