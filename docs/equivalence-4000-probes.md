# Analysis 2: R versus Python, on the existing fixture

Run as `tests/gen_equivalence_fixtures.py` followed by
`tests/test_equivalence.R`, both unmodified. Fixture: 4,000 probes, 39 samples,
12 subjects, seed 20240607 — seeded pseudo-random arrays, not methylation data,
as the generator's own docstring states. No 100,000-probe set exists; see
`docs/figure1-report.md` section 6.

| Quantity | Measure | Value | Component giving it | Tolerance | Checks |
|---|---|---|---|---|---|
| Probe-level effect estimates | max absolute difference | 2.22e-16 | `fit_within [raw] beta` | 1e-11 | 2 |
| Standard errors | max absolute difference | 8.73e-16 | `moderate_variance [mom vs python]` | 1e-10 | 3 |
| Cluster assignments | mismatching probes | 0 | `comethylation_clusters ids` | 0 | 1 |
| Segment boundaries | mismatching breakpoints | 0 | `tv_denoise breakpoints` | 0 | 1 |
| Region-level statistics | max absolute difference | 1.885e-14 | `call_regions p_pointwise` | 1e-09 | 6 |

Each row reports the largest value among the checks in that category, and the
component that produced it. The per-category detail:

- **Probe-level effects**: `beta` 2.22e-16 raw and 2.22e-16 shrunk.
- **Standard errors**: `se` 8.327e-17 raw and 8.327e-17 shrunk; the variance
  moderation agrees with the Python method-of-moments path to 8.728e-16, and
  with limma to a maximum relative difference of 2.677e-14, which the test
  records as informational rather than as a tolerance check.
- **Cluster assignments**: 181 clusters over 3,954 clustered probes, no probe
  assigned differently.
- **Segment boundaries**: 169 segments, no breakpoint placed differently; the
  denoised track itself agrees to 1.665e-16 and the distance penalty to
  1.421e-14.
- **Region-level statistics**: 141 regions in both; `chr`, `start`, `end`,
  `n_probes`, `probe_start` and `probe_end` identical; `effect_M` 1.665e-16,
  `se` 2.082e-17, `z` 1.066e-14, `theta` 1.665e-16, `p_pointwise` 1.885e-14,
  `max_abs_z` 1.066e-14.

Fifty-five checks ran in total and none failed. The full console output is
`equivalence_run.txt`; the five-row table is `equivalence_table.csv`.

## Scope of the comparison

The test compares `bin/ewasml.R` against `bin/ewasml.py` on identical input,
which is a question about the port. It is not a comparison on methylation data
and not a check of either implementation against a reference answer: both sides
could implement the same mistake and every difference above would still be at
machine precision. What the two implementations agree on is what the paper may
claim; whether the estimator is right is a separate question, which the design
rationale and the permutation calibration address.

The differences are at the scale of double-precision rounding rather than of
algorithmic divergence, and the two discrete quantities — which probe belongs to
which cluster, and where a segment breaks — agree exactly.
