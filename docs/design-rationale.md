# Design rationale

Why this pipeline is shaped the way it is. Read this before changing an
estimator: the choices below are not stylistic, and several of them are the
only options that leave the exposure effect identified.

## 1. In a two-cohort 450K + EPIC study, array type is not adjustable

Illumina arrays are processed in chip batches, and a chip carries one array
type. So in a study that recruited one cohort onto 450K and a later cohort
onto EPIC:

```
Subject  ──nested in──>  Cohort  ──nested in──>  Array type
                            │
                            └──nested in──>  Sentrix chip
```

Array type, cohort and chip are **mutually collinear**. There is no design
matrix that separates them, and putting `Array_Type` in the model produces a
rank-deficient design — `03_baseline_bumphunter.R --adjust-array` exists only
to demonstrate this rather than to be used.

The usual reflexes do not rescue the situation:

- **Adjusting for array type** is the collinearity above.
- **Batch correction** (ComBat and relatives) requires the batch and the
  effect of interest to be separable. Here "batch" *is* the cohort, so
  removing it removes any between-cohort component of the exposure effect too,
  silently.
- **Restricting to one array** discards the study.

### What is estimable

Each subject sits wholly on one array, so a **subject intercept absorbs array
type, cohort and chip together and exactly** — no residual array term is left
to adjust for. The exposure coefficient is then identified from repeated
visits of the same person, and the estimator never compares a 450K
measurement with an EPIC measurement.

This is the identifying assumption of the whole pipeline. Everything below is
a consequence of it.

## 2. Three consequences, enforced in code

**Only time-varying covariates are usable.** Age at baseline, sex, genotype
and any other subject-constant variable is annihilated by the within-subject
transform. Including them is not conservative, it is a rank problem.
`tests/test_ewasml.py` asserts the annihilation directly. Age and sex still
belong in QC — predicted-sex concordance is a per-sample check — but not in
`--covars`.

**Permutation must shuffle the exposure only among visits of the same
subject.** Freely permuting the exposure across all samples, as
`bumphunter`'s design-column permutation does, reassigns a subject's exposure
value to a different subject — hence to a different cohort, chip and array.
The resulting null contains between-array contrasts that no real reassignment
of visit dates could produce. It is not a null for this design. §5 of
[`results-GSE237561.md`](results-GSE237561.md) measures the size and the
direction of the error.

**Cross-validation must hold out whole subjects.** Splitting one subject's
visits across folds puts the within-subject contrast being validated on both
sides of the split, so the held-out score is optimistically biased by exactly
the quantity under selection.

## 3. What harmonisation costs: probe content, not chemistry

Measured directly from the two GEO platform manifests (`GPL13534`,
`GPL21145`); see `results/GSE237561/manifest_concordance.csv`.

| | probes |
|---|---|
| shared by CpG name | 454,181 |
| 450K only | 31,396 |
| EPIC only | 413,745 |

For all 454,181 shared probes, Infinium design type, both probe sequences,
colour channel and hg19 position are **100 % concordant** — zero discordant on
every one of those fields. Every single `AddressA_ID` differs, so an
address-level join returns nothing: the intersection has to be on CpG name,
and the manifest resolves addresses per array.

Two conclusions follow. First, **no type-I/II re-typing correction is
required**, because no shared probe changes design type between arrays; BMIQ
is therefore off by default rather than mandatory. Second, the cost of
harmonisation is *which* CpGs you can study, not how well you measure them —
and that cost falls unevenly, which is what breaks the legacy region and block
finders.

## 4. Why the legacy distance thresholds do not survive harmonisation

**Spacing.** On the harmonised autosomal set the median inter-probe gap is
335 bp, and only **48.2 %** of neighbouring pairs fall within 300 bp — the
`clusterMaker` default. A fixed 300 bp gap therefore fragments the harmonised
panel, and a "cluster" becomes a property of which array panel you happen to
be using rather than of the locus. Measured on GSE237561: `clusterMaker` at
its default returns 220,883 clusters over 422,898 probes, **median one probe
per cluster**, and more than half of `bumphunter`'s 214 reported regions are
single probes. The region finder is not aggregating.

**Open sea.** 78.2 % of EPIC-only probes are open sea, against 36.4 % of
shared probes. Harmonisation discards 323,589 EPIC open-sea probes and leaves
162,471 shared autosomal open-sea probes. Open sea is thus the content that
differs *most* between the two arrays — and it is exactly the content that
block finding relies on. A fixed-width collapse over it produces units whose
definition depends on which array the sample came from. Measured:
`cpgCollapse`'s 500 bp / 1500 bp rule gives 202,407 units of median width
1 bp on this panel.

The replacements therefore make distance an **upper bound** rather than the
rule, and let residual co-methylation decide what belongs together. A cluster
then reflects the locus's correlation structure, which is a property of the
biology, instead of the panel's spacing, which is a property of the product.

## 5. Why piecewise-constant rather than smoothed

`loessByCluster` uses a span measured in *probes*, which treats a 50 bp step
and a 50 kb step as equally close neighbours, and it returns a smooth curve —
an object with no breakpoints. Thresholding a curve gives boundaries that are
resolution-limited artefacts of the smoother's bandwidth, and the minfi paper
says as much about block boundaries.

Weighted total-variation denoising instead returns an explicitly
piecewise-constant track, so the breakpoints are the solver's own sparse
variable rather than a threshold crossing. Two design details matter:

- the fusion penalty decays as `exp(-d/decay)`, so genomic distance enters the
  penalty rather than the neighbour count;
- the penalty is scaled by the **median precision weight**, which makes the
  smoothness parameter dimensionless. One grid transfers across M-values, beta
  values and cohorts without retuning.

Each segment is reported as its precision-weighted mean (the relaxed-lasso
convention), so a region effect is an unbiased weighted mean in M-value units
and is directly comparable with a per-probe coefficient.

## 6. Why the comparison needed a shared estimator

The first old-vs-new comparison appeared to favour the legacy region finder on
cross-array agreement (r = 0.380 vs 0.191). That comparison was invalid, for
a reason worth recording because it is easy to repeat: each method reported
its own per-array effects, and the legacy per-array fit has no subject term
while the replacement's does. The gap mixed **which units were selected** with
**how the per-array effect was estimated**.

`06_compare.py` now re-estimates every method's per-array region effect from
the same within-subject per-probe fits, inverse-variance weighted across each
region interval. Under that shared estimator the two correlations are
indistinguishable (0.397 vs 0.388), and the difference that remains is in
sign concordance (0.729 vs 0.821) and in what the units *are*: single probes
versus 3-probe, 84 bp regions.

Read `n_cross_array_common` alongside any of these correlations. A
correlation over the 4 blocks a block method returns is not the same evidence
as one over 200 regions, and the figure annotates the unit count for that
reason.
