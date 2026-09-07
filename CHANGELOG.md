# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [Unreleased]

### Added
- `bin/ewasml.R` and `bin/04_dmr_ml.R`: the §2.6 region finder ported to R,
  now the pipeline default (`--dmr_impl r`). The Python implementation stays
  in the tree and is selectable with `--dmr_impl python`.
- `tests/test_equivalence.R` + `tests/gen_equivalence_fixtures.py`: the R core
  must reproduce the Python core on generated inputs — effects and region
  statistics to 1e-15, cluster ids, breakpoints and region boundaries exactly.
  Run in CI without limma, whose comparison then skips.
- `within_df()` in both cores: residual degrees of freedom of the
  within-subject fit, `n_samples - n_subjects - rank(within design)`.
- `fold_assign()` in both cores: cross-validation folds assigned by a hash of
  `(seed, subject id)` instead of an RNG shuffle, so the selected smoothness
  is identical across languages and library versions.
- `p_fwer_within_mcse` (and `p_fwer_naive_mcse`) in `dmr_ml.csv`: Monte Carlo
  standard error of each permutation p-value, with a log line counting regions
  within two standard errors of 0.05.
- `--var-method limma|mom` (R only) and `--cv_folds` as a Nextflow parameter.

### Fixed
- Cross-validation, stability selection and per-array replication skipped
  subsets by subject count alone; a fold of mostly single-visit subjects has
  no within-subject information and crashed the stage. All three now test
  identification with `within_df` and skip with a counted log line; the stage
  errors only when no fold is estimable.
- The Galaxy region-finder wrapper ran the Python script with the Python
  requirement set; it now runs R with `requirements_r_ml` (r-base, limma,
  Matrix, optparse, jsonlite — neither minfi nor bumphunter).

## [Unreleased]

### Fixed

- `nextflow.config` parses under the strict config parser introduced in
  Nextflow 25. Three problems, all of which only appeared once the workflow was
  run through Nextflow itself rather than stage by stage: a `def` timestamp
  declaration (rejected — "variable declarations cannot be mixed with config
  statements", the nf-core idiom that worked under 23/24), quoted resource
  literals such as `'8.GB'` (a quoted value is parsed as a memory or duration
  string, where the valid forms are `'8 GB'` and `'8 h'`; the dotted form is
  only valid unquoted), and a top-level `workflow.onComplete { }` in `main.nf`
  (now a config closure, which still reports on a failed run).
- `params.outdir` is declared in `nextflow.config`, which references it for the
  timeline, report, trace and DAG paths. It was declared only in `main.nf`, and
  because the config is evaluated before the script, the audit trail was written
  under `null/pipeline_info/` unless `--outdir` was given explicitly.

### Changed

- CI checks workflow syntax with `nextflow lint` rather than
  `nextflow run . --help`, which executed the workflow and failed on the
  required `--sheet` guard instead of testing syntax. Profile and lint failures
  are emitted as `::error::` annotations so the message is visible in the run
  summary without opening the log.

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
