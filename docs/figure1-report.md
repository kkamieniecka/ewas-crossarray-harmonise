# Figure 1 task: inventory, spacing baseline, and two stop items

Run against the recorded GSE237561 pipeline outputs. No preprocessing was
rerun, nothing was simulated, and no manuscript number was re-derived or
re-rounded.

## 1. Files found

The workspace was cleaned between sessions, so the run outputs were restored
from the archived result bundles rather than from disk.

| What the task asked for | Where it is | Count | Agrees with the task |
|---|---|---|---|
| Harmonised matrix probe IDs | `results/harmonised/mval_dims.json` (in `ewas-harmonise-results.tar.gz`) | 422,898 probes x 126 samples | yes |
| Probe filter record | `results/harmonised/probe_filter.csv` | 452,567 rows | pre-filter set |
| bumphunter regions | `results/baseline/bumphunter_regions.csv` | 214 | yes |
| Distance-aware regions | `results/dmr/dmr_ml.csv` | 207 | yes |
| Per-array region effects | `effect_450K` / `effect_EPIC` columns in both region tables | 214 and 207 | yes |
| Per-array probe effects | `results/probe_model/probe_stats_450K.csv`, `_EPIC.csv`, `_joint.csv` | 422,898 each | yes |
| blockFinder blocks | `results/baseline/blockfinder_blocks.csv` | 4 | n/a |
| HSMM blocks | `results/blocks/blocks_hsmm.csv` | 8 (comparison run); 38 in the full run | see note |

**The harmonised M-value matrix itself (`harmonised.rds`) is not in any
archive.** Only its dimensions and row names survive. The task allowed "the
harmonised matrix or its probe IDs", so the inventory condition is met, but
nothing that needs the values can be recomputed.

## 2. The exact bumphunter call

From `bin/03_baseline_bumphunter.R` (lines 40, 80, 85) and confirmed in
`results/baseline.log`:

```r
cl <- clusterMaker(chr = ann$chr_s, pos = ann$pos, maxGap = 300)   # --max-gap default 300L
bh <- bumphunter(M, design = design, chr = ann$chr_s, pos = ann$pos, cluster = cl,
                 coef = 2, cutoff = 0.10, B = 200, smooth = TRUE,
                 smoothFunction = loessByCluster, nullMethod = "permutation")
```

So the comparison used `bumphunter::clusterMaker` with `maxGap = 300`, which is
the value in the EWASGalaxy section 2.6 workflow, not minfi's `bumphunter()`
default. The log records `clusterMaker(maxGap=300): 220883 clusters; median
probes/cluster = 1`.

Two things the log shows that are worth carrying into the paper. minfi itself
warned that `nullMethod = "permutation"` is not recommended with more than two
design columns, and this design has nine; the baseline was run that way
deliberately, because it is what the existing workflow does, but a reviewer who
runs it will see the same warning. And the blockFinder path reports `collapsed
units: 202407 (median width 1 bp)`, which is where the "median width 1 bp"
figure comes from: it describes the collapsed units, not the four called
blocks, whose median width is 2,460.75 bp.

## 3. Where 24,552 clusters over 56,623 probes comes from

`results/dmr.log`: `co-methylation clusters: 24552 clusters covering 56623
probes (13.4%)`. Stage 04 clusters by correlation, not by distance alone:
`max_gap = 1000`, `rho_min = 0.3`, `decay_bp = 1000` (`results/dmr/run_config.json`).
A probe joins a cluster only if a neighbour within 1 kb is correlated with it
above 0.3, so 366,275 of the 422,898 probes are left out and the clusters that
do form average 2.3 probes.

The relation to 220,883 is that the two numbers count different things under
different rules on the same panel. `clusterMaker(maxGap = 300)` partitions
**all** 422,898 probes by distance alone and returns 220,883 clusters, 70.4% of
which hold a single probe. The co-methylation rule returns 24,552 clusters over
the 13.4% of probes that have a correlated neighbour. Neither is a subset of
the other, and the paper should not present them as comparable cluster counts.

## 4. Probe spacing baseline

One filter rule for all three sets, as the task specifies: autosomes only, and
loci with a SNP at the CpG or single-base-extension position at MAF >= 0.01
removed, which is the `dropLociWithSnps(snps = c("SBE", "CpG"), maf = 0.01)`
call `bin/01_harmonise.R` makes. Gaps are distances between neighbouring probes
within a chromosome; clusters are `clusterMaker`'s rule at `maxGap = 300`,
reimplemented in Python and validated against the run (see below).

| Probe set | Probes | Median gap | IQR | Gaps > 300 bp | Clusters | Single-probe clusters |
|---|---|---|---|---|---|---|
| Native 450K | 456,741 | 307 bp | 80–2,471 | 50.4% | 229,977 | 70.1% |
| Native EPIC | 816,578 | 514 bp | 108–2,517 | 58.6% | 478,852 | 75.4% |
| Shared, both manifests | 425,939 | 342 bp | 86–2,687 | 52.2% | 222,351 | 70.4% |
| Analysed panel, as run | 422,898 | 343 bp | 87–2,696 | 52.2% | **220,883** | 70.4% |

The cluster count on the analysed panel reproduces the run exactly (220,883),
which validates the reimplementation of the clustering rule.

Two filter notes. Detection-p filtering (4,735 probes) is sample-dependent and
cannot be applied to a native manifest without IDATs, which no longer exist, so
the native rows carry the two filters the task names and not that one.
Cross-reactive filtering was **not active in the recorded run**:
`probe_filter.csv` has `is_crossreactive = 0` for all 452,567 probes, because
`bin/01_harmonise.R` only applies that list when `EWAS_CROSSREACTIVE` names a
file and the variable was unset. If the methods section claims Chen or Pidsley
cross-reactive probes were removed, it is claiming something the run did not do.

## 5. Stop item 1: median gap 335 bp is not reproducible as stated

The task requires the shared-panel row to give median gap 335 bp and 220,883
clusters. Those two figures do not describe the same probe set.

| Value | Probe set that produces it |
|---|---|
| 335 bp | Shared manifest, autosomes, **SNP probes retained** (443,483 probes, from `crossarray_probe_map.csv.gz`) |
| 342 bp | Shared manifest, autosomes, SNP-at-CpG/SBE removed (425,939 probes) |
| 343 bp | The analysed panel as run (422,898 probes), which is the set carrying the 220,883 clusters |

So 335 bp is measured before SNP filtering and 220,883 clusters after it. Each
figure is correct for its own set; the sentence that pairs them is not. Nothing
was substituted. Asked which set the row should denote, you left the choice to
me, so panel B and the table above use **the analysed panel (422,898 probes,
median gap 343 bp)**, on the grounds that it is the only set that carries both
quantities and the only one the paper actually analyses. The manuscript's
335 bp should be corrected to 343 bp, or qualified as the pre-SNP-filter figure
for a 443,483-probe set that was never analysed.

## 6. Stop item 2: the 100,000-probe test set does not exist

The equivalence fixtures are generated by `tests/gen_equivalence_fixtures.py`,
which writes `N_PROBES = 4000`, `N_SUBJECTS = 12`, `SEED = 20240607`, and its
own docstring states the fixtures are seeded pseudo-random arrays rather than
methylation data. `tests/test-data/` holds only a README, so nothing is stored.
There is no 100,000-probe set anywhere in the repository or the result bundles,
and the harmonised matrix that could supply 100,000 real probes is gone.

Running the comparison at 4,000 probes is possible and would exercise both
implementations on identical input. Generating a 100,000-probe set would mean
simulating data, which the task forbids. Analysis 2 therefore stops here.

## 7. What was produced

| File | Contents |
|---|---|
| `spacing_baseline.csv` | The four-row table in section 4 |
| `probes_450k_native.csv`, `probes_epic_native.csv`, `probes_shared_manifest.csv` | The filtered probe sets with positions |
| `gaps.npz` | The gap vectors for panel B, one array per set |
| `panelC_probes_per_region.csv` | Probes per region, 1/2/3/4/5+, both methods |
| `panelD_bumphunter.csv`, `panelD_distance-aware.csv` | Per-region 450K and EPIC effects on the shared estimator |
| `provenance.json` | Every figure above with the file and line it came from |
| `01_probe_sets.R` | The probe-set script |
| `sessionInfo_R.txt`, `versions_python.txt` | Environment records |

Panels C and D are computed and reproduce the manuscript exactly: bumphunter
r = 0.397 and sign concordance 72.9% over 214 regions, distance-aware
segmentation r = 0.388 and 82.1% over 207. Both come from re-estimating each
region's per-array effect from the same within-subject per-probe fits, as
`bin/06_compare.py` does, which is the `*_common` pair in
`method_comparison.csv`. Note that the `effect_450K` and `effect_EPIC` columns
stored in the region tables are a different quantity, each method's own
per-array refit, and they give r = 0.380 with 67.3% concordance for bumphunter
and r = 0.191 with 64.7% for segmentation. Panel D must use the shared
estimator or it will not match the manuscript.

## 8. Figure 1

Drawn as `figure1.pdf` (vector, 178 mm wide), `figure1.tiff` (350 dpi, LZW) and
`figure1.png` (screen copy), by `02_figure1.py` from the files in section 7.
Panel B uses the analysed panel, per section 5. The legend and the alt-text
description are in `figure1_legend.md`, together with a point-by-point check
against the specification. Panel D prints the manuscript's own figures, which
the computation reproduced to three decimal places.

## 9. Not done

Analysis 2, for the reason in section 6: no 100,000-probe test set exists and
the harmonised matrix is gone. The equivalence comparison can be run at 4,000
probes on the existing seeded fixture if that is acceptable; it cannot be run at
100,000 without generating data.
