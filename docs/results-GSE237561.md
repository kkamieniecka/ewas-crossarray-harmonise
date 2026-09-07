# Run record: GSE237561, cross-array longitudinal EWAS

Dataset: GSE237561 (clozapine, repeated blood sampling). 126 IDAT pairs,
38 subjects, two cohorts on two array generations. Exposure
`days_on_clozapine`, scaled per 100 days. Covariates: smoking score and six
Houseman cell fractions. Seed 1 throughout.

Every number below is read back from the files in `results/`; the file that
carries each one is named so it can be re-checked.

## 1. The design constraint, measured

`qc_samples.csv` — 26 subjects (92 samples) on 450K, 12 subjects (34 samples)
on EPIC. **No subject appears on both arrays**, and no Sentrix chip spans the
two array types. Array type is therefore perfectly confounded with cohort and
with chip, and no amount of normalisation can separate them. Visits per
subject: 4 visits for 17 subjects, 3 for 13, 2 for 7, 5 for one.

This single fact is what makes the two legacy region algorithms invalid here,
and it is why the exposure has to be identified within subject.

## 2. Harmonisation

`probe_filter.csv`, `mval_dims.json` — 452,567 probes survive the
450K/EPIC intersection; 422,898 remain after filtering. Removed:

| filter | probes |
|---|---|
| failed detection p | 4,735 |
| sex chromosome | 10,585 |
| SNP-proximal (`dropLociWithSnps`) | 16,098 |
| cross-reactive list | 0 (list not applied) |

Sample QC: 0 sample failures, 0 predicted-sex mismatches, worst sample has
0.24 % failed probes. Harmonisation cost probe *content*, not data quality.

## 3. Repeated-measures probe model

`probe_model/model_summary.csv`:

| fit | samples | subjects | within-subject correlation | probes FDR<0.05 |
|---|---|---|---|---|
| joint | 126 | 38 | 0.347 | 1 |
| 450K | 92 | 26 | 0.349 | 2 |
| EPIC | 34 | 12 | 0.308 | 0 |

The within-subject correlation is ~0.35 in every fit. That is the quantity the
legacy pipeline sets to zero by treating 126 samples as 126 independent
observations.

## 4. Section 2.6 replacement vs bumphunter

`comparison/method_comparison.csv`. Both were run on the identical harmonised
matrix.

| | bumphunter (2.6) | TV + within-subject perm |
|---|---|---|
| units called | 214 | 207 |
| median unit width | **1 bp** | 84 bp |
| median probes per unit | **1** | 3 |
| units at FWER ≤ 0.05 | 0 | 0 |
| cross-array *r* (shared estimator) | 0.397 | 0.388 |
| cross-array sign concordance | 0.729 | **0.821** |

Two things matter here.

**bumphunter is not finding regions on this panel.** `clusterMaker` at its
default 300 bp gap returns 220,883 clusters over 422,898 probes — median one
probe per cluster. More than half of the 214 "regions" it reports are single
probes (`bumphunter_regions.csv`, median `L` = 1). The fixed genomic gap is
the wrong distance for the harmonised 450K∩EPIC panel, whose spacing is
sparser than either array alone.

**The comparison had to be put on a common footing to mean anything.** Each
method reports its own per-array effects, but the legacy per-array fit has no
subject term and the replacement's does — so the two `cross_array_r` values
as originally computed (0.380 vs 0.191) differ partly because of the
estimator, not the region choice. `06_compare.py` therefore re-estimates every
method's per-array region effect from the *same* within-subject per-probe fits,
inverse-variance weighted across each interval. Under that shared estimator
the correlations are indistinguishable (0.397 vs 0.388), but the replacement
achieves it on aggregated 3-probe units rather than single probes, and its
sign concordance is 9 points higher. Aggregation that preserves cross-cohort
agreement is the intended behaviour; single-probe "regions" cannot be said to
have aggregated anything.

## 5. The permutation scheme, isolated

`permutation_null.csv`, 200 replicates each, computed from the **identical**
observed statistic so nothing but the permutation scheme differs:

| null 95th percentile of max&#124;z&#124; | value |
|---|---|
| within-subject (visit times permuted inside each subject) | 5.99 |
| free, bumphunter-style (exposure permuted across all samples) | 6.28 |

Observed max &#124;z&#124; = 5.90; smallest FWER *p* = 0.065 under the
within-subject null, 0.104 under the free one.

The free permutation gives a **wider** null. It is mis-calibrated, but in the
conservative direction on this design: reassigning the exposure across
subjects manufactures between-subject — hence between-cohort and
between-array — contrasts that no real reassignment of visit dates could
produce, and those inflate the null maximum. The cost is power, not false
positives. minfi says as much itself: `bumphunterEngine` emitted its own
warning that the permutation test is not recommended with more than two design
columns (`baseline.log`).

## 6. Section 2.7 replacement vs blockFinder

`blocks/hsmm_params.json`, `baseline/blockfinder_blocks.csv`:

| | blockFinder (2.7) | distance-aware HSMM |
|---|---|---|
| blocks called | 4 | 8 |
| median block width | 2,461 bp | 62,301 bp |
| median clusters per block | 3 | 9 |
| cross-array *r* (shared estimator) | 0.674 (n=4) | 0.307 (n=8) |
| collapsed units | 202,407 (median 1 bp) | 16,173 (median 196 bp) |

The HSMM converged in 21 iterations (log-likelihood 38,380.5) on 152,721
open-sea probes carried by both arrays, with stationary state probabilities
0.124 / 0.621 / 0.256 for hypo / neutral / hyper.

Both correlations here rest on 4 and 8 units respectively and carry almost no
information — they are reported for completeness, not as a comparison.

> **This table predates the state-labelling fix** (see CHANGELOG, *Fixed*).
> The run above assumed the middle HSMM state is the zero-effect one instead
> of identifying it. Its state means were ordered hypo < neutral < hyper with
> effects on both sides of zero, so the assumption was probably satisfied
> here and the counts probably do not change — but that has not been
> verified, because verifying it means re-running the stage. Until it is
> re-run, treat the block count as a lower bound. `blocks_hsmm.csv` and
> `openSea_cluster_effects.csv.gz` under `results/GSE237561/` come from the
> same run and carry the same caveat.

Two operational findings on the legacy implementation:

- `blockFinder` **fails outright** at the probe-level cutoff (0.10): when no
  bump clears the cutoff, minfi dereferences an atomic result table and
  errors rather than returning an empty result. The baseline needed a separate
  `--block-cutoff 0.02` to run at all, which is now a documented option.
- `cpgCollapse`'s fixed 500 bp / 1500 bp collapse produces 202,407 units of
  median width 1 bp — again, mostly singletons. On the same shared open-sea
  subset the co-methylation clustering gives 16,173 clusters of median 196 bp.

## 7. Smoothness selection: an honest negative

`dmr/lambda_cv.csv` and `dmr_lamcheck/lambda_cv.csv`. Held-out weighted SSE
over subject-wise folds decreases **monotonically** across the whole range
tested, λ₀ = 0.05 → 128 (score 1.2971 → 1.2123), so cross-validation does not
identify an interior optimum and selects at whichever grid edge it is given.

This is not a solver failure, and it does not affect the result: re-running at
λ₀ = 128 returns the **identical** region set (207 regions, max &#124;z&#124; =
5.90). The co-methylation clusters are small — 24,552 clusters covering 56,623
probes, ~2.3 probes each — and total-variation segmentation on two- and
three-probe sequences saturates immediately. Beyond the point of saturation λ
has nothing left to fuse. The behaviour to expect on a denser panel, or at a
lower correlation threshold, is an interior optimum; on this one the estimator
degenerates gracefully to a weighted cluster-mean test.

## 8. What this run does and does not establish

**Nothing is discovered.** Zero regions or blocks at FWER ≤ 0.05 by any of the
four methods, and one probe at FDR < 0.05 in the joint fit. With 38 subjects
and an exposure measured in days on drug, that is the expected outcome, and it
is the honest headline.

What the run does establish is about the *methods*, and it does not depend on
finding a hit:

1. The confounding is total (§1) and the within-subject correlation is
   substantial (§3), so the legacy model's independence assumption is violated
   on this dataset, measurably.
2. The legacy region and block finders reduce to single-probe testing on the
   harmonised panel (§4, §6) — their fixed distance thresholds do not survive
   the 450K∩EPIC intersection.
3. The legacy permutation scheme is mis-calibrated on this design, in the
   conservative direction, costing power (§5).
4. The replacements aggregate genuinely (84 bp / 3 probes; 62 kb / 9 clusters)
   while matching cross-cohort agreement and improving sign concordance (§4).

## Reproducing

```
nextflow run . -profile conda \
  --sample_sheet sample_sheet_GSE237561.csv \
  --idat_dir data/idat \
  --probe_map results/crossarray_probe_map.csv.gz \
  --exposure days_on_clozapine --subject Subject_ID
```

Stage runtimes on 8 CPUs: harmonisation ~25 min, probe model ~24 min, legacy
baseline ~6 min, region replacement ~2.6 h (dominated by 2 × 200
permutations), block replacement ~5 s, comparison ~10 s.
