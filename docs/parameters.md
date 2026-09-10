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

## Region detection, §2.6 replacement (`04_dmr_ml.R`, `04_dmr_ml.py`)

Two implementations, identical flags, proved equivalent in
`tests/test_equivalence.R` (see design-rationale §7). R is the default;
`--dmr_impl python` runs the Python one. R costs about 2x runtime and 2x
memory.


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
| `--var-method` | `limma` | **R only.** `limma` uses `squeezeVar`; `mom` uses the ported method-of-moments moderation. The two agree to 2.7e-14, so this is an audit switch, not a modelling choice |
| `--seed` | `1` | seeds the permutation and subsample draws. It does **not** seed the CV fold split, which is a deterministic hash of `(seed, subject id)` so that R and Python select the same smoothness; changing the seed still changes the split |

Setting `--n-perm 0` is a legitimate diagnostic mode (segmentation only, no
inference); the run record omits the null quantiles rather than failing.

`dmr_ml.csv` reports `p_fwer_within_mcse`, the Monte Carlo standard error
`sqrt(p(1-p)/n_perm)` of each permutation p-value, and the stage logs how many
regions lie within two standard errors of 0.05. When correlated regions share
a similar |z| near the threshold, the count passing at 0.05 moves with the
null draw and not with the data: at `--n-perm 200` on a test input the two
implementations reported 8 and 1 significant regions from indistinguishable
nulls. Treat that log line as an instruction to raise `--n-perm` (1000 or
more for publication) rather than as a result.

The Nextflow parameters that choose the implementation are `--blocks_impl`,
`--dmr_impl`
(`r` | `python`) and `--var_method` (`limma` | `mom`); the Galaxy wrapper runs
R only, and exposes the moderation choice as *Variance moderation*.

## Block detection, §2.7 replacement (`05_blocks_hsmm.R`, `05_blocks_hsmm.py`)

Two implementations, identical flags. Baum-Welch from fixed starting values
uses no RNG, so the two agree exactly — identical clusters, blocks and
directions, posteriors to 5e-16, log-likelihood to 2e-15 relative
(`tests/test_stage05_equivalence.py`). R is the default; `--blocks_impl python`
runs the Python one, which is how the agreement is re-checked on real data.

| param | default | notes |
|---|---|---|
| `--max-gap` | `1500` | upper bound on gap within an open-sea cluster |
| `--rho-min` | `0.20` | lower than the region threshold: open-sea correlation is weaker |
| `--length-scale` | `250000` | HMM transition length scale `L` in `A(d) = e^{-d/L} I + (1-e^{-d/L}) 1π'` |
| `--min-post` | `0.80` | posterior threshold for calling a block |
| `--min-clusters` | `3` | minimum clusters per block |
| `--min-fit-clusters` | `200` | precondition on the whole panel: below this many open-sea clusters the stage refuses to fit rather than returning an unidentifiable three-state model. Distinct from `--min-clusters`, which is per called block. Lower it only for small test panels — the Galaxy tool test sets `100` for a 4000-probe fixture |
| `--fixed-collapse` | off | also report what the legacy fixed 500/1500 bp collapse would have produced on the same probes |
| `--var-method` | `limma` | R only, as in stage 04: `limma` (squeezeVar) or `mom`, the ported method-of-moments path the Python uses. The stage-level equivalence check runs `mom` on both sides, which is also what makes it runnable without limma |
| `--array-col` | `Array_Type` | column splitting the per-array check. An arm with no residual degrees of freedom after the within transform is declined rather than fitted, and `cross_array_r` is then `null` |

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
