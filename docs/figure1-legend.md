# Figure 1 legend and alt text

## Legend

**Figure 1. Harmonisation of 450K and EPIC arrays, probe spacing, and the
regions each detection method calls.** (**A**) Workflow. Raw intensity files
from both platforms pass per-array quality control and normalisation, are
reduced to the probes the two designs share, and are analysed under a
within-subject model before region detection by bumphunter and by
co-methylation segmentation. (**B**) Distance to the neighbouring probe within
a chromosome, for the two native designs and for the analysed shared panel
(422,898 probes); the vertical line marks the 300 bp threshold at which
`clusterMaker` starts a new cluster. Gaps below 10 bp are pooled into a single
wide bin, whose height is divided by that width as for every other bin. Median
gaps are 307 bp on native 450K, 514 bp on native EPIC and 343 bp on the shared
panel. Probe sets are filtered identically: autosomes only, loci with a SNP at
the CpG or single-base-extension position at minor allele frequency 0.01 or
above removed. (**C**) Regions called by each method and setting, split by the
number of probes they contain, with the total above each bar: bumphunter at
`maxGap` 300, 500 and 1000 bp, and co-methylation segmentation, whose minimum
region size is three probes (`--min-probes 3`, with `--min-effect 0.05`), so it
calls no one- or two-probe region by construction. Bar position names the
method and setting; the greyscale value and hatch name the size class.
(**D**) Region effect estimated separately in the 450K and EPIC cohorts, one
point per region, restricted to regions of three or more probes so the two
methods are compared on the same footing: bumphunter at `maxGap` 1000 (42
regions) and co-methylation segmentation (207 regions). Both axes use the same
within-subject per-probe fits, aggregated over each region by inverse-variance
weighting, so the only difference between the two scatters is which regions the
method selected. The number of regions, the Pearson correlation and the
proportion agreeing in direction are given above each scatter. Effects are on
the M-value scale per 100 days of exposure. No region reached family-wise
significance under either method at any setting shown.

## Alt text

A four-panel figure. Panel A is a flow diagram: 450K IDATs from 92 arrays and
EPIC IDATs from 34 arrays feed into per-array quality control and
normalisation, then into a shared-probe panel of 422,898 probes, then into a
within-subject model of 38 subjects and 126 arrays, then into region detection
by bumphunter and by co-methylation segmentation. Panel B is a density plot of
the distance between neighbouring probes on a logarithmic axis from 1 base pair
to 100 kilobases, with three curves distinguished by greyscale value and dash
pattern: native 450K, native EPIC and the shared panel. All three rise to a
peak near 100 base pairs and have a second broader peak between 1 and 10
kilobases; the EPIC curve sits lower at the first peak and higher between 1 and
10 kilobases, indicating wider spacing. A vertical line at 300 base pairs marks
the clustering threshold and falls between the two peaks. Gaps shorter than 10
base pairs are pooled into one wide bin at the left, which reads as a low flat
shelf. Panel C is a stacked bar chart with four bars: bumphunter at 300, 500
and 1000 base pairs, and co-methylation segmentation. Each bar is split into
regions of one probe, two probes, and three or more probes, shown in three
greyscale values with different hatching, and the total is printed above. The
bumphunter bars total 214, 337 and 532 regions and are dominated by the
one-probe class, which holds 157, 247 and 400 of them. The segmentation bar
totals 207 regions, all of which are in the three-or-more class. Panel D is a
pair of scatter plots of the EPIC effect against the 450K effect for each
region, with reference lines at zero on both axes and the same scale from minus
0.4 to plus 0.4; each carries its statistics as a heading above the plotting
area rather than inside it. The left plot, bumphunter at 1000 base pairs in
orange circles, shows 42 regions with correlation 0.581 and 90.5 per cent
agreeing in direction. The right plot, co-methylation segmentation in blue
triangles, shows 207 regions with correlation 0.388 and 82.1 per cent agreeing
in direction. Panels C and D share a baseline and height.

## Where the numbers come from

| Quantity | File |
|---|---|
| Panel B medians and the probe-set counts | `results/GSE237561/spacing_baseline.csv` |
| Panel C counts | `results/GSE237561/bumphunter_regions_maxgap{300,500,1000}.csv` and the 207-region co-methylation table |
| Panel D points and statistics | the same region tables, re-estimated per array from `results/probe_model/probe_stats_450K.csv` and `..._EPIC.csv` by inverse-variance weighting, as `bin/06_compare.py` does |
| Intervals for the same statistics | `results/GSE237561/table1_cross_array_ci.csv` |

## Specification check

- Width 178 mm, vector PDF plus 350 dpi LZW TIFF, and a screen PNG.
- The methods are named "bumphunter" and "co-methylation segmentation"
  throughout the figure, the legend and the alt text.
- Panel C distinguishes its three size classes by greyscale value and by hatch;
  panel D distinguishes its two methods by colour (Okabe-Ito orange and blue)
  and by marker shape. No distinction in the figure depends on hue alone.
- Panel D prints n, r and sign concordance above each scatter, clear of the
  points: 42 regions with r 0.581 and 90.5%, and 207 regions with r 0.388 and
  82.1%.
- Panels C and D have the same height and the same baseline.
