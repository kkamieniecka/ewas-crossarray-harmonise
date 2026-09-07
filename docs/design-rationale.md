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

## 7. Why the region finder exists twice, in R and in Python

§2.6 was first written in Python because the estimator is banded sparse linear
algebra and a total-variation solve. Nothing in it needs Python: the whole
numerical dependency surface is dense/banded linear algebra, a trigamma
inverse and a group-wise mean, all of which R has natively or through
`Matrix`. Since the rest of the suite — and the group that maintains it — is
R, the region finder was ported (`bin/ewasml.R`, `bin/04_dmr_ml.R`) and R is
now the default. The Python implementation stays in the tree so agreement can
be re-checked at any time with `--dmr_impl python`.

### What "equivalent" was made to mean

`tests/test_equivalence.R` generates inputs and reference outputs by importing
the Python core directly, then requires R to reproduce them. Every
deterministic quantity is held to a numerical tolerance, and every
integer-valued or identity-valued one to exact equality:

| quantity | agreement |
|---|---|
| within-subject effect, SE, z | 1e-16 |
| variance moderation (`mom` path) | 2e-15 |
| co-methylation cluster ids | exact |
| distance penalty, denoised track | 1e-15 |
| segmentation breakpoints | identical set |
| region boundaries, probe index ranges | exact |
| region effect, SE, z | 1e-15 |
| residual df, fold assignment | exact |

Two things cannot be compared draw for draw, because the two languages have
different pseudo-random generators: the permutation null and the
stability-selection subsamples. These are handled two ways. The arithmetic
downstream of the draws is made comparable by feeding R the *reference* draws
from the fixture, and the sampling schemes themselves are checked by the
properties that define them — within-subject permutation preserves each
subject's own multiset of exposures and leaves single-visit subjects fixed;
free permutation preserves the overall multiset; stability selection draws
subjects without replacement at the intended size.

On a 100k-probe input the two implementations return the same 98 clusters, the
same 200 regions and the same max |z| = 144.97.

### limma is the default moderation, and that is now a measurement

The Python core hand-rolled an empirical-Bayes variance moderation. Reading it
suggested it was a reimplementation of `limma::squeezeVar`; the R port made
that testable, and the two agree to 2.7e-14. The R implementation therefore
calls limma by default (`--var-method limma`) and keeps the ported
method-of-moments path (`--var-method mom`) for the equivalence check.

### What the port cost, and what it bought

Cost, measured on a 100k-probe x 126-sample input and linear in probe count,
so roughly 4.2x these numbers at EPIC scale:

| | runtime | peak RSS |
|---|---|---|
| Python | 37 s | 0.96 GiB |
| R | 80 s | 1.78 GiB |

R is about twice as slow and takes about twice the memory; `r_heavy` is sized
for that. What the port bought was three defects that two independent
implementations of the same estimator made visible and that one implementation
had hidden:

1. **The cross-validation guard did not check identification.** Both
   implementations failed identically on a fold whose subjects were mostly
   single-visit: the guard counted subjects but never asked whether the fold
   retained residual degrees of freedom after subject intercepts. The
   condition `n_samples - n_subjects - rank(within design) > 0` is now a
   first-class core function (`within_df`) applied at all three places a
   subset is fitted — CV folds, stability subsamples, per-array replication —
   with skipped folds counted in the log rather than silently dropped.
2. **The fold split was not portable.** Each language's RNG produced a
   different assignment, so the selected smoothness differed between
   implementations on identical input. Folds are now assigned by an explicit
   hash of `(seed, subject id)` (`fold_assign`), which is reproducible across
   languages, R versions and NumPy versions. The seed still gives an
   independent split.
3. **The FWER decision was less stable than it looked.** At 200 permutations
   the two implementations reported 8 and 1 regions at FWER <= 0.05 from
   statistically indistinguishable nulls (KS p = 0.47). The cause is not a
   bug: seven regions had |z| between 14.0 and 15.6, exactly where the 0.05
   cut falls, so an ordinary Monte Carlo wobble in the null tail flips them
   together. Both implementations now report `p_fwer_within_mcse`, the Monte
   Carlo standard error of each permutation p-value, and log how many regions
   sit within two of them of the threshold. A bootstrap over the null
   replicates confirms the per-region binomial error is the right scale.

Point 3 is a property of the method, not of the port: any permutation FWER on
correlated regions has it. It was invisible until two implementations
disagreed. Publication runs should use `--n-perm 1000` or more and read the
MCSE column before treating a region near the threshold as significant.
