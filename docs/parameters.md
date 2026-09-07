# Parameter reference

Two ways to set anything: as a Nextflow `--param` (the workflow forwards it to
the stage), or directly on a stage script, each of which takes `--help`.
Where a workflow default differs from a script default, **the workflow default
is what runs** and is the one listed here.

## Required inputs

| Nextflow param | Meaning |
|---|---|
| `--sheet` | sample sheet, CSV (see below) |
| `--idat_dir` | directory containing the IDAT pairs named in `Basename` |
| `--probe_map` | cross-array probe map from the manifest comparison, `.csv.gz` |
| `--outdir` | output root (default `results`) |

### Sample sheet columns

A standard minfi sheet plus the longitudinal columns:

| column | required | notes |
|---|---|---|
| `Sample_Name` | yes | unique per array |
| `Basename` | yes | IDAT path stem, without `_Grn.idat` |
| `Sentrix_ID`, `Sentrix_Position` | yes | chip and position |
| `Array_Type` | yes | `450K` or `EPIC` |
| `Subject_ID` | yes | **the identifying variable** — repeated across visits |
| `visit` | yes | visit index or label |
| *exposure column* | yes | named by `--exposure` |
| `sex` | recommended | used for the predicted-sex QC check |
| time-varying covariates | optional | named by `--covars` |

## Design

| param | default | notes |
|---|---|---|
| `--exposure` | `days_on_clozapine` | **dataset-specific default** — set this for your own study |
| `--exposure_scale` | `100` | exposure divided by this, so coefficients read per 100 units |
| `--time_scale` | `per100days` | label carried into the probe-model output |
| `--subject` | `Subject_ID` | the within-subject grouping |
| `--covars` | `smoking_score,cd8t,cd4t,bcell,gran,mono,nk` | **must be time-varying**; constant covariates are annihilated by the within transform |
| `--seed` | `1` | seeds permutation, bootstrap and CV folds |

Age and sex belong in QC, not in `--covars`: they do not change within a
subject over the study, so the within-subject contrast removes them exactly.
There is a unit test asserting this.

## Harmonisation (`01_harmonise.R`)

| param | default | notes |
|---|---|---|
| `--detp` | `0.01` | detection p threshold for calling a probe failed in a sample |
| `--detp_frac` | `0.05` | a probe is dropped if it fails in more than this fraction of samples |
| `--drop_sex` | `true` | drop chrX/chrY probes |
| `--bmiq` | `false` | off deliberately: every shared probe has the same Infinium design type, so there is no type-I/II imbalance introduced by harmonisation |
| `--keep_gset` | `true` | also write a `GenomicRatioSet`, needed by the legacy baseline |

## Region detection, §2.6 replacement (`04_dmr_ml.py`)

| param | default | notes |
|---|---|---|
| `--max-gap` | `1000` | maximum bp between probes considered for the same cluster — an upper bound, not the cluster rule |
| `--rho-min` | `0.30` | probes join a cluster only if residual co-methylation exceeds this. This, not distance, is the cluster rule |
| `--decay-bp` | `1000` | length scale of the `exp(-d/decay)` fusion penalty |
| `--lam-grid` | `0.05,…,4.0` | smoothness grid, dimensionless (scaled by the median precision weight), selected by held-out-subject CV |
| `--cv-folds` | `5` | folds hold out **whole subjects** |
| `--min-probes` | `3` | minimum probes for a reported region |
| `--min-effect` | `0.05` | minimum segment effect, M-value units |
| `--n-perm` | `200` | within-subject permutations for the FWER null |
| `--n-boot` | `50` | subject subsamples for stability selection |
| `--also-naive-perm` | off | additionally run the legacy free permutation on the same statistic, to quantify its mis-calibration |

Setting `--n-perm 0` is a legitimate diagnostic mode (segmentation only, no
inference); the run record omits the null quantiles rather than failing.

## Block detection, §2.7 replacement (`05_blocks_hsmm.py`)

| param | default | notes |
|---|---|---|
| `--max-gap` | `1500` | upper bound on gap within an open-sea cluster |
| `--rho-min` | `0.20` | lower than the region threshold: open-sea correlation is weaker |
| `--length-scale` | `250000` | HMM transition length scale `L` in `A(d) = e^{-d/L} I + (1-e^{-d/L}) 1π'` |
| `--min-post` | `0.80` | posterior threshold for calling a block |
| `--min-clusters` | `3` | minimum clusters per block |
| `--fixed-collapse` | off | also report what the legacy fixed 500/1500 bp collapse would have produced on the same probes |

## Legacy baseline (`03_baseline_bumphunter.R`)

Run for comparison only. Defaults reproduce minfi's documented workflow.

| param | default | notes |
|---|---|---|
| `--max-gap` | `300` | `clusterMaker` default |
| `--cutoff` | `0.10` | probe-level bump cutoff |
| `--block-cutoff` | `0.02` | **separate, lower** cutoff for `blockFinder`. Required: at the probe-level cutoff no collapsed bump clears the threshold and minfi errors instead of returning an empty result |
| `--block-window` | `250000` | loess window for block smoothing |
| `--n-perm` | `200` | free permutations, matching the region stage's count |
| `--adjust-array` | off | add `Array_Type` to the design. Produces a rank-deficient design in a nested two-cohort study; exposed only to demonstrate that |
| `--skip-blocks` | off | skip §2.7 |

## Comparison (`06_compare.py`)

| param | notes |
|---|---|
| `--probe-model-dir` | **pass this.** It enables the `*_common` columns, which re-estimate every method's per-array effect from one shared within-subject estimator. Without it the cross-array columns are not comparable between methods |
| `--no-figure` | skip the figure |
