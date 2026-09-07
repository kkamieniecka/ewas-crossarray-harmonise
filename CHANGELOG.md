# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [0.1.0] - 2026-09-05

First working version. Validated end to end on GSE237561 (126 arrays,
38 subjects, two cohorts across 450K and EPIC); see
[`docs/results-GSE237561.md`](docs/results-GSE237561.md).

### Added

- `01_harmonise.R` — per-array QC in native probe space, then merge through
  `minfi::combineArrays()`. Exports a raw `float64` M-value matrix so the
  numerical stages need no R dependency, and optionally a `GenomicRatioSet`
  so the legacy baseline runs on the identical matrix.
- `02_probe_model.R` — repeated-measures probe-level model with a subject
  term, variance moderation, and per-array fits for cross-array checks.
- `03_baseline_bumphunter.R` — the legacy minfi §2.6/§2.7 path, run on the
  same harmonised matrix as the replacements so the comparison is controlled.
- `04_dmr_ml.py` — §2.6 replacement: co-methylation clustering, weighted
  total-variation segmentation, held-out-subject cross-validation,
  within-subject permutation, stability selection, cross-array generalisation.
- `05_blocks_hsmm.py` — §2.7 replacement: distance-aware 3-state HMM over
  open-sea co-methylation clusters with soft boundaries and per-block
  posteriors.
- `06_compare.py` — controlled old-vs-new benchmark. Re-estimates every
  method's per-array region effect from a single shared within-subject
  estimator, so the cross-array comparison reflects unit selection rather
  than estimator choice.
- `ewasml.py` — numerical core, with 17 property-based tests in
  `tests/test_ewasml.py`.
- Nextflow DSL2 workflow with `conda`, `docker`, `singularity` and `test`
  profiles.
- Galaxy wrappers for the harmonisation, region and block tools, sharing a
  `macros.xml` that defines the longitudinal design inputs once.

### Known limitations

- Cross-validation over the smoothness parameter finds no interior optimum on
  the GSE237561 panel; the region set is invariant across the whole range
  tested, so this does not affect the result. See §7 of the results document.
- The cross-reactive probe list is not yet applied (the filter column is
  present and always `FALSE`).
- No Galaxy wrapper for `02_probe_model.R` yet; it runs as a standalone
  `Rscript` stage inside Nextflow.
- Galaxy tool tests have no committed fixtures — see `tests/test-data/`.
