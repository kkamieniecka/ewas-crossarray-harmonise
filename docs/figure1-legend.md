# Figure 1 legend and alt text

## Legend

**Figure 1. Harmonisation of 450K and EPIC arrays and the effect of
distance-aware region detection.** (**A**) Workflow. Raw intensity files from
both platforms pass per-array quality control and normalisation, are reduced to
the probes the two designs share, and are analysed under a within-subject model
before region detection. (**B**) Distance to the neighbouring probe within a
chromosome, for the two native designs and for the analysed shared panel
(422,898 probes); the vertical line marks the 300 bp threshold at which
`clusterMaker` starts a new cluster. Gaps below 10 bp are pooled into a single
wide bin, whose height is divided by that width as for every other bin. Median
gaps are 307 bp on native 450K, 514 bp on native EPIC and 343 bp on the shared
panel; 52.2% of shared-panel gaps exceed the threshold, and 70.4% of the
resulting 220,883 clusters contain a single probe. Probe sets are filtered identically: autosomes only, loci with
a SNP at the CpG or single-base-extension position at minor allele frequency
0.01 or above removed. (**C**) Probes per called region. bumphunter returns 157
single-probe regions of 214; distance-aware segmentation returns none, because
its minimum region size is three probes (`--min-probes 3`, with
`--min-effect 0.05`), as the key records. (**D**) Region effect estimated
separately in the 450K and EPIC cohorts, one point per region. Both axes use the
same within-subject per-probe fits, aggregated over each region by
inverse-variance weighting, so the only difference between the two scatters is
which units the method selected. Pearson correlation, the number of regions and
the proportion agreeing in direction are given above each scatter. Effects are on the M-value
scale per 100 days of exposure. No region reached family-wise significance in
either method.

## Alt text

A four-panel figure. Panel A is a flow diagram: boxes for 450K intensity files
from 92 arrays and EPIC intensity files from 34 arrays both feed into per-array
quality control and normalisation, then into a shared-probe panel of 422,898
probes, then into a within-subject model of 38 subjects across 126 arrays, and
finally into region detection by either bumphunter or distance-aware
segmentation. Panel B plots three density curves of the distance between
neighbouring probes on a logarithmic axis from about 3 base pairs to 100
kilobases: native 450K as a dashed dark grey line, native EPIC as a dotted light
grey line, and the shared panel as a solid black line. All three are bimodal,
with a peak near 100 base pairs and a second broader peak between 1 and 10
kilobases; the EPIC curve sits lower at the first peak and higher between 1 and
10 kilobases, indicating wider spacing. A vertical line at 300 base pairs marks
the clustering threshold and falls between the two peaks, so roughly half of all
gaps lie on either side of it. Gaps shorter than 10 base pairs are pooled into
one wide bin at the left, which reads as a low flat shelf. Panel C is a grouped bar chart of the number of
regions by probe count, with categories one, two, three, four and five or more.
bumphunter, in orange, is concentrated at one probe with about 157 regions and
declines steeply. Distance-aware segmentation, in blue with diagonal hatching,
has nothing at one or two probes, about 141 regions at three probes, 45 at four
and 21 at five or more. Panel D is a pair of scatter plots of the EPIC effect
against the 450K effect for each region, with reference lines at zero on both
axes and equal scales from minus 0.4 to plus 0.4; each carries its statistics as
a heading above the plotting area rather than inside it. The left plot,
bumphunter in orange circles, shows 214 regions in a diffuse cloud with
correlation 0.397 and 72.9 per cent agreeing in direction. The right plot,
distance-aware segmentation in blue triangles, shows 207 regions in a tighter
cloud with correlation 0.388 and 82.1 per cent agreeing in direction. Panels C
and D share a baseline and height.

## Specification check

- Width 178 mm, vector PDF plus 350 dpi LZW TIFF, and a screen PNG.
- Methods are named "bumphunter" and "distance-aware segmentation" throughout.
- Methods carry both a colour (Okabe-Ito orange and blue) and a shape or fill
  cue (circles against triangles in D, plain against hatched bars in C); the
  three probe sets in B are distinguished by greyscale value and dash pattern,
  so no distinction in the figure depends on hue alone.
- The two method colours are used for nothing else in the figure.
- Panel D prints the values the manuscript reports: r = 0.397 and 72.9% for
  bumphunter, r = 0.388 and 82.1% for distance-aware segmentation.
