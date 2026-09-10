# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [Unreleased]

### Fixed
- The Galaxy tool tests could not run at all: `galaxy/ewas_dmr_ml.xml` and
  `galaxy/ewas_blocks_hsmm.xml` referenced five `test-data` files that were
  never committed, so `planemo test` failed at collection. The fixtures are
  now in `galaxy/test-data/`, generated from
  `tests/gen_equivalence_fixtures.py` (4000 synthetic probes, 39 arrays, 12
  subjects — no study data).
- `05_blocks_hsmm.R` refused to fit below a hard-coded 200 open-sea clusters,
  a floor no option could reach and one that no small panel can clear, which
  made the block tool untestable. The floor is now `--min-fit-clusters`
  (default `200`, so production behaviour is unchanged) and is exposed in the
  wrapper; the tool test sets `100`. Its error message also advised relaxing
  `--rho-min`/`--max-gap`, which merges probes and yields *fewer* clusters —
  it now names the direction that actually splits them and reports both counts.
- `ewas_blocks_hsmm.xml` published the open-sea cluster table as `format="csv"`
  while handing Galaxy the gzipped file, so the job failed while setting
  metadata (`UnicodeDecodeError` on the gzip magic byte) and no output was
  collected. The wrapper now decompresses that table before Galaxy picks it
  up; the Nextflow stage still writes it gzipped. Found by `planemo test`, and
  invisible to the serverless runner, which does not check datatypes.
- `-profile test` skipped the comparison stage entirely: `COMPARE` was inside
  the `if (params.run_baseline)` branch, so the smoke test stopped one stage
  short of the table it exists to produce. `COMPARE` now always runs and is
  handed an empty baseline channel when the legacy stage is off;
  `06_compare.py --baseline-dir` became optional and reports the replacement
  methods only in that case, with the row set and the common-footing loop
  driven off one spec list so they cannot fall out of step.
- `COMPARE` was never given the per-probe model directory, so the pipeline as
  wired could not produce the `*_common` columns — the only cross-array
  numbers that are comparable between methods. `PROBE_MODEL.out` is now an
  input and `--probe-model-dir` is passed.
- Published outputs landed one directory too deep
  (`results/dmr/dmr/dmr_ml.csv`): every process already emits its stage
  directory, and `publishDir` appended it again. Publishing is now flat to
  `params.outdir` across all eight processes.
- `workflow.onComplete` in `nextflow.config` called `log.info`, which is not
  bound in that scope; the resulting handler error masked the real failure on
  every failed run. It prints instead.
- The `test` profile inherited the 32 GB `r_heavy` memory request, which the
  local executor refuses outright on a smaller machine, so the run failed
  before submitting any work. The profile now sets its own resource requests.
- `06_compare.py` falls back to a plain pipe table when pandas' optional
  `tabulate` dependency is absent, instead of failing at the last stage.

### Added
- `planemo lint` and `planemo test` run in CI (`galaxy-tool-tests` job) on the
  two wrappers that carry `<tests>`, with dependency resolution off against a
  micromamba environment holding the R requirements, and the planemo report
  uploaded as a build artifact. Linting fails on errors only, so the
  `TestsMissing` warning on `ewas_harmonise` (it needs IDATs) stays visible
  without failing the build.
- `tests/run_galaxy_tool_tests.py`: runs a tool's `<tests>` without a Galaxy
  server, for sandboxes where planemo cannot bind a local port. It renders the
  `<command>` with Cheetah from the test values plus the XML defaults, executes
  it, and checks `expect_num_outputs` and `<assert_contents>`. It does not
  replace planemo — no datatype, metadata, output-format or dependency-
  resolution checks.
- `bin/05_blocks_hsmm.R` and the block model in `bin/ewasml.R`
  (`dist_transitions`, `fit_block_hsmm`, `call_blocks`, `state_labels`): the
  §2.7 block finder ported to R, now the pipeline default (`--blocks_impl r`).
  The Python implementation stays in the tree (`--blocks_impl python`). The
  HSMM draws no random numbers, so the two agree exactly, not in distribution:
  identical clusters, blocks and directions, posteriors to 5e-16.
- `tests/test_stage05_equivalence.py` and `--stage-dir` in
  `tests/gen_equivalence_fixtures.py`: a stage-level check that runs both
  stage-05 drivers end to end and compares every output file, including the
  case where an array arm is unidentified. Both defects below were in driver
  code, out of reach of the function-level fixtures.
- `state_labels()` in both cores, plus `state_labels` and `implementation` in
  `hsmm_params.json`, so the two records compare field for field.
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
- `--var-method limma|mom` (R only) on both R stages, and `--cv_folds` as a
  Nextflow parameter. Stage 05 initially shipped without the flag, which made
  it silently require limma; the stage-level check runs `mom` on both sides,
  which is the comparison that isolates the port and keeps CI limma-free.

### Fixed
- Block direction was labelled by assuming the zero-effect HSMM state is the
  middle of the three. When all cluster effects fall on one side of zero the
  state pinned at mu = 0 sorts to an end, and the middle state is then a
  genuine effect state reported as "no change". On a fixture forcing that
  ordering the fix takes the call from 10 blocks to 22. Both implementations
  now label relative to the identified neutral state and log which single
  direction remains callable when it sorts to an end. `blocks_hsmm.csv` and
  `openSea_cluster_effects.csv.gz` produced before this fix understate the
  block count; see `docs/results-GSE237561.md`.
- `hsmm_params.json` was not always valid JSON: a declined array arm or an
  empty legacy fixed-width selection put a bare `NaN` (Python) or the string
  `"NA"` (R) in the file, so no strict parser could read it and the two
  implementations disagreed on how an absent value is spelled. Both now
  serialise every non-finite value as `null`.
- Stage 05's per-array check tested subject count alone before fitting an arm;
  it now requires residual degrees of freedom after the within transform, the
  same guard stage 04 uses, and declines the arm with a log line instead of
  returning numbers from an unidentified fit.
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
