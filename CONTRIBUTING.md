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
   `results/` are documentation exhibits and must stay small.
5. **Record parameters, not just outputs.** Every stage writes a run record
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
