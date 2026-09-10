# Contributing

## Scope

This repository is a methods pipeline for longitudinal, cross-array EWAS. The
design rests on one identifying assumption — the exposure effect is estimated
**within subject** — and contributions that break it will not be merged.
Read [`docs/design-rationale.md`](docs/design-rationale.md) first; it explains
why array type cannot be adjusted for as a covariate in this design.

## Before opening a pull request

1. **Run the numerical tests.** `python tests/test_ewasml.py`. All 17 checks
   must pass. Each states its property and tolerance, so a failure should
   localise to one estimator.
2. **Add a test with any estimator change.** The tests are property-based
   (recovery of a known effect, annihilation of a time-invariant covariate,
   the flat limit of the penalty, permutation calibration), not
   value-regression. Assert the property your change is supposed to preserve.
3. **Keep the legacy baseline runnable.** `03_baseline_bumphunter.R` exists so
   claims about the replacements are measured on the same matrix rather than
   asserted. Do not "fix" the legacy path to behave better — its behaviour on
   a harmonised panel is the finding.
4. **Do not commit data.** No IDATs, no full probe-level matrices, no
   `harmonised.rds`. `.gitignore` covers the usual paths. Result tables under
   `results/` are documentation exhibits and must stay small. The one binary
   in the tree is `galaxy/test-data/test_mval.f64` (1.2 MB), a synthetic
   4000-probe fixture the tool tests cannot generate at run time — regenerate
   it rather than editing it, and do not grow it.
5. **Run the Galaxy tool tests if you touched `galaxy/`.** The reference check
   is planemo, which CI runs on every push:

   ```sh
   planemo lint --fail_level error galaxy/ewas_*.xml
   planemo test --no_dependency_resolution \
       galaxy/ewas_dmr_ml.xml galaxy/ewas_blocks_hsmm.xml
   ```

   `planemo test` starts a Galaxy instance on a local port, which some
   development sandboxes forbid. Where it cannot run, use

   ```sh
   python tests/run_galaxy_tool_tests.py galaxy/ewas_dmr_ml.xml \
       galaxy/ewas_blocks_hsmm.xml
   ```

   which renders the same `<command>` from the same `<test>` values plus the
   XML defaults, runs it, and checks `expect_num_outputs` and
   `<assert_contents>` — no server, and no substitute for planemo: it does not
   check datatypes, metadata, output formats or dependency resolution. Both
   paths need `Rscript` with `limma`, `Matrix`, `optparse` and `jsonlite` on
   the path. Fixtures live in `galaxy/test-data/` and are regenerated with
   `python tests/gen_equivalence_fixtures.py --stage-dir <dir>`; they are
   synthetic by design, so no study data enters the repository.
6. **Record parameters, not just outputs.** Every stage writes a run record
   (`run_config.json`, `hsmm_params.json`, `baseline_summary.json`) with
   resolved parameters, counts, seed and runtime. New parameters belong there
   too.

## Style

- R: `optparse` for arguments, `data.table` for tabular work, explicit
  `set.seed()`.
- Python: standard library plus numpy/scipy/pandas. `ewasml.py` must stay free
  of R and of Bioconductor annotation dependencies so it is unit-testable in
  isolation.
- Comments should explain *why* a choice was made, especially where it differs
  from minfi. The reasons are the contribution.

## Reporting a problem

Open an issue with the stage that failed, the resolved run record if one was
written, and the relevant log from `results/*/logs/`. If the problem is
statistical rather than mechanical, say what the design is — number of
subjects, visits per subject, and how array type maps onto cohort and chip.
