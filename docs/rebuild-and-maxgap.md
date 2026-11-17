# Rebuild of the harmonised matrix, and bumphunter at three gap values

Steps 1 to 5 of the brief, in order, stopping at nothing. No extra
permutations and no simulation were run. The results are reported without a
judgement about which method is preferable.

## 1. IDAT reader test

minfi 1.56.0 under R 4.5.3 read one pair from each platform:

| Platform | Addresses | Detection-p probes | noob |
|---|---|---|---|
| 450K (GSM7625568) | 622,399 | 485,512 | 485,512 x 1 |
| EPIC (GSM7625664) | 1,051,815 | 866,091 | 866,091 x 1 |

## 2. Fetch

`tools/fetch_idats.py --jobs 6` against the URLs in
`data/sample_sheets/sample_sheet_GSE237561.csv`: **252 of 252 files, 0
failures, 1.3 GB**, each verified as a complete gzip stream. Measured sizes are
8.65 MB per 450K sample and 14.65 MB per EPIC sample.

## 3. Rebuild and verification

`bin/01_harmonise.R` with the recorded parameters (`--detp 0.01 --detp_frac
0.05 --drop_sex TRUE --bmiq FALSE --keep_gset TRUE --seed 20240717`):
**6 min 47 s** wall, peak vector memory 10.3 GB, 452,567 probes combined and
422,898 after filtering, 0 QC failures and 0 sex mismatches.

Against the recorded run:

| Check | Result |
|---|---|
| Matrix dimensions | 422,898 x 126, as recorded |
| Row names against `mval_dims.json` | identical, and in the same order |
| Column names against `mval_dims.json` | identical, and in the same order |
| `probe_filter.csv`, all four flags over 452,567 probes | identical; zero disagreements |
| `qc_samples.csv` | 126 samples, 0 failures, 0 sex mismatches |

Checksum of the regenerated M-value matrix:
`sha256 3845f277e9651ec6a9acdb9fe3ab5dc4acde5d7ce15825d5d41f99c33a620ed9`
(`rebuild/harmonised/mval.f64`, 407 MB).

The rebuild ran under minfi 1.56.0 where the recorded run used 1.52.1. The
probe set and every filter decision are unchanged across that difference.

## 4 and 5. bumphunter at maxGap 300, 500 and 1000

`bin/03_baseline_bumphunter.R` on the rebuilt matrix, everything except
`--max-gap` as recorded: `--covars smoking_score,cd8t,cd4t,bcell,gran,mono,nk`,
`--cutoff 0.10`, `--n-perm 200`, `--seed 1`, blocks skipped. maxGap 300 is
included as the control.

Cross-array agreement is the shared estimator: each region's per-array effect
re-estimated from the same within-subject per-probe fits
(`results/probe_model/probe_stats_450K.csv` and `..._EPIC.csv`) by
inverse-variance weighting, which is what `bin/06_compare.py` computes as
`cross_array_r_common`. Those per-probe fits are the recorded ones, which the
matrix verification above licenses.

Intervals are 95%: Fisher z for the correlation, back-transformed, and Wilson
for the proportion agreeing in direction. Both treat regions as independent
observations, which they are not exactly — overlapping regions share probes —
so these are the nominal intervals rather than a design-aware bound.

| Method | maxGap | Clusters | Region set | Regions | Median width | Median probes | Cross-array r (95% CI) | Sign concordance (95% CI) |
|---|---|---|---|---|---|---|---|---|
| bumphunter | 300 | 220,883 | all | 214 | 1 bp | 1 | 0.397 (0.278, 0.505) | 0.729 (0.666, 0.784) |
| bumphunter | 500 | 186,782 | all | 337 | 1 bp | 1 | 0.484 (0.398, 0.562) | 0.745 (0.696, 0.788) |
| bumphunter | 1000 | 153,440 | all | 532 | 1 bp | 1 | 0.444 (0.373, 0.510) | 0.735 (0.696, 0.771) |
| bumphunter | 300 | 220,883 | >= 3 probes | 27 | 165 bp | 4 | 0.353 (-0.031, 0.646) | 0.741 (0.553, 0.868) |
| bumphunter | 500 | 186,782 | >= 3 probes | 35 | 186 bp | 4 | 0.464 (0.154, 0.690) | 0.800 (0.641, 0.900) |
| bumphunter | 1000 | 153,440 | >= 3 probes | 42 | 235.5 bp | 3.5 | 0.581 (0.337, 0.752) | 0.905 (0.779, 0.962) |
| co-methylation segmentation | n/a | 24,552 | all (>= 3 probes by construction) | 207 | 84 bp | 3 | 0.388 (0.266, 0.498) | 0.821 (0.763, 0.867) |

The segmentation row is the recorded 207-region run, recomputed here with the
same shared estimator; it reproduces the published 0.388 and 0.821. Its
clusters come from co-methylation clustering rather than from `clusterMaker`,
so that count is not comparable with the three above it.

Every correlation interval in the table overlaps every other, and the one at
maxGap 300 restricted to three or more probes includes zero. Among the
concordance intervals, bumphunter at maxGap 1000 restricted to three or more
probes (0.779, 0.962) and the segmentation row (0.763, 0.867) overlap each
other and both overlap the all-region bumphunter rows.

No region reached family-wise significance at any of the three values
(`at fwer<=0.05: 0` in all three logs).

The maxGap 300 row reproduces the published bumphunter row of
`method_comparison.csv` exactly: 214 regions, median width 1 bp, median 1
probe, r 0.3973, sign concordance 0.7290.

Probes per region, by count:

| maxGap | 1 probe | 2 | 3 | 4 | >= 3 |
|---|---|---|---|---|---|
| 300 | 157 | 30 | 10 | 6 | 27 |
| 500 | 247 | 55 | 17 | 7 | 35 |
| 1000 | 400 | 90 | 21 | 7 | 42 |

Raising maxGap raises the region count and leaves the median region at one
probe at every value tested; the share of multi-probe regions moves from 27 of
214 to 42 of 532. Within the subset of regions carrying at least three probes,
both cross-array statistics rise with maxGap.

## Files

| File | Contents |
|---|---|
| `bumphunter_maxgap_sweep.csv` | the seven-row table above, with the interval bounds and the counts they were computed from (`n_cross_array`, `n_same_sign`) |
| `baseline_300/`, `baseline_500/`, `baseline_1000/` | each run's `bumphunter_regions.csv` and `baseline_summary.json` |
| `baseline_300.log`, `baseline_500.log`, `baseline_1000.log` | the run logs, including the clusterMaker line |
| `harmonise.log` | the rebuild log |
| `harmonised/` | the regenerated matrix and its side files (not committed; 1.9 GB) |
