#!/usr/bin/env python
"""Figure 1 for the cross-array EWAS applications note.

Four panels, 178 mm wide, written as vector PDF and 350 dpi TIFF.

Inputs, all produced by fig1/01_probe_sets.R and the inventory step, none of
them recomputed here:
    spacing_baseline.csv            per-probe-set spacing summary
    gaps.npz                        the gap vectors behind panel B
    panelC_probes_per_region.csv    probes per region, both methods
    panelD_bumphunter.csv           per-region 450K and EPIC effects,
    panelD_distance-aware.csv       shared within-subject estimator

Usage:  python fig1/02_figure1.py <dir>
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SEED = 1
MM = 1 / 25.4
WIDTH = 178 * MM

# Okabe-Ito, reserved for the two methods and used for nothing else
BH_C, TV_C = "#E69F00", "#0072B2"
BH_M, TV_M = "o", "^"
# probe sets are a different entity class, so they take greys plus dash codes
SET_STYLE = {"native 450K": ("#4d4d4d", (0, (4, 1.6))),
             "native EPIC": ("#9a9a9a", (0, (1.2, 1.4))),
             "shared panel": ("#111111", "solid")}
MAXGAP = 300
# --min-probes as recorded in results/GSE237561/run_config.json for the 207 regions
MIN_PROBES = 3


def panel_a(ax):
    """Workflow schematic, in the same words the data panels use."""
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    box = dict(boxstyle="round,pad=0.34", linewidth=0.8,
               facecolor="#f5f5f5", edgecolor="#555555")
    arrow = dict(arrowstyle="-|>", linewidth=0.8, color="#555555",
                 shrinkA=1.5, shrinkB=1.5)

    def put(x, y, text, **kw):
        return ax.text(x, y, text, ha="center", va="center", fontsize=6.6,
                       bbox={**box, **kw}, zorder=3)

    a1 = put(0.22, 0.90, "450K IDATs\n92 arrays")
    a2 = put(0.70, 0.90, "EPIC IDATs\n34 arrays")
    qc = put(0.46, 0.685, "per-array QC and normalisation")
    sp = put(0.46, 0.475, "shared-probe panel\n422,898 probes", facecolor="#e8eef5")
    wm = put(0.46, 0.272, "within-subject model\n38 subjects, 126 arrays")
    rd = put(0.46, 0.075, "region detection", facecolor="#e8eef5")
    for src in (a1, a2):
        ax.annotate("", xy=(0.46, 0.745), xytext=(src.get_position()[0], 0.845),
                    arrowprops=arrow)
    for y0, y1 in [(0.625, 0.545), (0.415, 0.34), (0.212, 0.145)]:
        ax.annotate("", xy=(0.46, y1), xytext=(0.46, y0), arrowprops=arrow)
    ax.text(0.46, -0.03, "bumphunter   |   distance-aware segmentation",
            ha="center", va="center", fontsize=6.2, color="#333333")
    return rd


def panel_b(ax, gaps, summary):
    # gaps below 10 bp go into one wide bin spanning 1-10 bp. Because the
    # histogram is a density, that bin's height is divided by its own width, so
    # pooling reads as a low wide bar instead of a spike at the left edge.
    edges = np.concatenate([[0.0], np.linspace(1.0, 5.2, 57)])
    for label, key in [("native 450K", "native_450K"), ("native EPIC", "native_EPIC"),
                       ("shared panel", "analysed_panel")]:
        g = gaps[key]
        g = g[g > 0]
        colour, dash = SET_STYLE[label]
        ax.hist(np.clip(np.log10(g), 0.0, None), bins=edges, density=True,
                histtype="step",
                color=colour, linestyle=dash, linewidth=1.1, label=label)
    ax.axvline(np.log10(MAXGAP), color="#555555", linewidth=0.8)
    ax.text(np.log10(MAXGAP) + 0.07, 0.035, "clustering threshold\n300 bp",
            transform=ax.get_xaxis_transform(), ha="left", va="bottom",
            fontsize=6.2, color="#333333")
    ax.set_xlim(-0.05, 5.2)
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(["1 bp", "10 bp", "100 bp", "1 kb", "10 kb", "100 kb"])
    ax.set_xlabel("gap to the neighbouring probe")
    ax.set_ylabel("density")
    med = {r["set"].split(" (")[0]: r["median_gap_bp"] for _, r in summary.iterrows()}
    ax.legend(frameon=False, fontsize=6.2, loc="upper left",
              handlelength=2.1, borderaxespad=0.1)
    return med


def panel_c(ax, pc, min_probes):
    x = np.arange(len(pc.index))
    w = 0.38
    ax.bar(x - w / 2, pc["bumphunter"], w, color=BH_C, edgecolor="white",
           linewidth=0.4, label="bumphunter")
    ax.bar(x + w / 2, pc["distance-aware segmentation"], w, color=TV_C,
           edgecolor="white", linewidth=0.4, hatch="///",
           label="distance-aware segmentation")
    for xi, v in zip(x - w / 2, pc["bumphunter"]):
        if v == 0:
            ax.plot([xi], [0], marker="_", color=BH_C, markersize=5)
    for xi, v in zip(x + w / 2, pc["distance-aware segmentation"]):
        if v == 0:
            ax.plot([xi], [0], marker="_", color=TV_C, markersize=5)
    ax.set_xticks(x)
    ax.set_xticklabels(pc.index)
    ax.set_xlabel("probes per region")
    ax.set_ylabel("regions")
    ax.set_ylim(0, float(np.max(pc.to_numpy())) * 1.32)   # headroom for the legend
    blank = mpl.lines.Line2D([], [], linestyle="none")
    h, lab = ax.get_legend_handles_labels()
    ax.legend(h + [blank], lab + ["minimum region size = %d probes" % min_probes],
              frameon=False, fontsize=6.2, loc="upper right", borderaxespad=0.1,
              handlelength=1.6, handletextpad=0.6)


def panel_d(axes, d_bh, d_tv, stats):
    for ax, (name, d, colour, marker) in zip(
            axes, [("bumphunter", d_bh, BH_C, BH_M),
                   ("distance-aware segmentation", d_tv, TV_C, TV_M)]):
        a = d.effect_450K_common.to_numpy()
        b = d.effect_EPIC_common.to_numpy()
        ok = np.isfinite(a) & np.isfinite(b)
        ax.axhline(0, color="#bbbbbb", linewidth=0.6, zorder=1)
        ax.axvline(0, color="#bbbbbb", linewidth=0.6, zorder=1)
        ax.scatter(a[ok], b[ok], s=7, facecolor=colour, edgecolor="none",
                   marker=marker, alpha=0.75, zorder=3)
        s = stats[name]
        ax.set_title("%s\nr = %.3f, n = %d\nsign concordance %.1f%%"
                     % (name.replace(" segmentation", "\nsegmentation"),
                        s["r"], s["n"], 100 * s["sign"]),
                     loc="left", fontsize=6.2, color="#222222", linespacing=1.5, pad=3.5)
        lim = 0.42
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
        ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4])
        ax.set_xlabel("450K effect (M)")
    axes[0].set_ylabel("EPIC effect (M)")
    axes[1].set_yticklabels([])


def gap_vectors(d):
    """Load gaps.npz, or rebuild it from the probe coordinate files.

    gaps.npz is 4.9 MB and regenerable, so it is not committed. The three
    native/shared sets come from 01_probe_sets.R; the analysed panel comes from
    results/GSE237561/panel_probe_positions.csv.gz, which is the probe
    coordinate column of the run's own probe_stats_joint.csv.
    """
    npz = os.path.join(d, "gaps.npz")
    if os.path.exists(npz):
        return np.load(npz)
    srcs = {"native_450K": "probes_450k_native.csv",
            "native_EPIC": "probes_epic_native.csv",
            "analysed_panel": "panel_probe_positions.csv.gz"}
    out = {}
    for key, fn in srcs.items():
        t = pd.read_csv(os.path.join(d, fn))
        t["chr_i"] = pd.to_numeric(t.chr.astype(str).str.replace("chr", "", regex=False),
                                   errors="coerce")
        t = t[t.chr_i.between(1, 22)].sort_values(["chr_i", "pos"], kind="mergesort")
        g = np.diff(t.pos.to_numpy(np.int64))
        out[key] = g[np.diff(t.chr_i.to_numpy()) == 0]
    np.savez_compressed(npz, **out)
    return np.load(npz)


def main(d):
    summary = pd.read_csv(os.path.join(d, "spacing_baseline.csv"))
    gaps = gap_vectors(d)
    pc = pd.read_csv(os.path.join(d, "panelC_probes_per_region.csv"), index_col=0)
    d_bh = pd.read_csv(os.path.join(d, "panelD_bumphunter.csv"))
    d_tv = pd.read_csv(os.path.join(d, "panelD_distance-aware.csv"))
    stats = {}
    for name, dd in [("bumphunter", d_bh), ("distance-aware segmentation", d_tv)]:
        a, b = dd.effect_450K_common.to_numpy(), dd.effect_EPIC_common.to_numpy()
        ok = np.isfinite(a) & np.isfinite(b)
        stats[name] = dict(n=int(ok.sum()),
                           r=float(np.corrcoef(a[ok], b[ok])[0, 1]),
                           sign=float((np.sign(a[ok]) == np.sign(b[ok])).mean()))

    fig = plt.figure(figsize=(WIDTH, 5.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.92], width_ratios=[1, 1.08],
                          hspace=0.42, wspace=0.26,
                          left=0.075, right=0.985, top=0.945, bottom=0.085)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    sub = gs[1, 1].subgridspec(1, 2, wspace=0.10)
    ax_d1 = fig.add_subplot(sub[0, 0])
    ax_d2 = fig.add_subplot(sub[0, 1])

    panel_a(ax_a)
    panel_b(ax_b, gaps, summary)
    panel_c(ax_c, pc, min_probes=MIN_PROBES)
    panel_d([ax_d1, ax_d2], d_bh, d_tv, stats)

    for ax, letter in [(ax_a, "A"), (ax_b, "B"), (ax_c, "C"), (ax_d1, "D")]:
        ax.text(-0.02 if ax is ax_a else -0.16, 1.06, letter, transform=ax.transAxes,
                fontsize=9, fontweight="bold", va="bottom", ha="left")

    fig.savefig(os.path.join(d, "figure1.pdf"), facecolor="white")
    fig.savefig(os.path.join(d, "figure1.tiff"), dpi=350, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(os.path.join(d, "figure1.png"), dpi=200, facecolor="white")
    json.dump(stats, open(os.path.join(d, "figure1_stats.json"), "w"), indent=1)
    return fig, stats


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fig1")
