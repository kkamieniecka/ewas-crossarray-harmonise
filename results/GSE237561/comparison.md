# Region and block detection: legacy vs replacement

Harmonised matrix: 422898 autosomal probes x 126 samples, 38 subjects.

| method                                     | level   |   n_regions |   n_fwer05 |   median_width_bp |   median_n_probes |   cross_array_r_common |   sign_concordance_common |   cross_array_r |   sign_concordance |   null_q95_max_abs_z |   total_bp |   n_cross_array |   runtime_s |   n_clusters |   n_cross_array_common |   n_fwer05_naive_perm |   null_q95_naive_perm |
|:-------------------------------------------|:--------|------------:|-----------:|------------------:|------------------:|-----------------------:|--------------------------:|----------------:|-------------------:|---------------------:|-----------:|----------------:|------------:|-------------:|-----------------------:|----------------------:|----------------------:|
| bumphunter (minfi 2.6)                     | region  |         214 |          0 |              1    |                 1 |               0.397339 |                  0.728972 |        0.380489 |           0.672897 |            nan       |       7575 |             214 |     337.991 |       220883 |                    214 |                   nan |              nan      |
| TV + within-subject perm (2.6 replacement) | region  |         207 |          0 |             84    |                 3 |               0.387972 |                  0.821256 |        0.190718 |           0.647343 |              5.98947 |      31633 |             207 |    9399.2   |        24552 |                    207 |                     0 |                6.2769 |
| blockFinder (minfi 2.7)                    | block   |           4 |          0 |           2460.75 |                 3 |               0.67386  |                  1        |      nan        |         nan        |            nan       |      11468 |             nan |     337.991 |       220883 |                      4 |                   nan |              nan      |
| distance-aware HSMM (2.7 replacement)      | block   |           8 |        nan |          62300.5  |                 9 |               0.307095 |                  0.875    |       -0.117672 |           0.625    |            nan       |     785131 |               8 |       4.4   |        16173 |                      8 |                   nan |              nan      |

## Permutation scheme

Null 95th percentile of max|z|: within-subject **5.99**, free (bumphunter-style) **6.28**. Both nulls were computed from the identical observed statistic, so the difference is attributable to the permutation scheme alone. Because subject, cohort, chip and array type are nested here, a free permutation of the exposure creates between-array contrasts that no real reassignment of visit times could produce; the resulting null is therefore not a null for this design.

## Cross-array generalisation

Both columns correlate a region's effect estimated in the 450K cohort alone against the same region estimated in the EPIC cohort alone -- the only metric here that a method cannot improve by simply calling more regions.

`cross_array_r` uses each method's own per-array estimator and is therefore **not comparable between methods**: the legacy per-array fit has no subject term while the replacement's does, so the difference mixes unit selection with estimator choice.

`cross_array_r_common` re-estimates every method's per-array region effect from the same within-subject per-probe fits (inverse-variance weighted over the probes inside each interval). The only remaining difference is which units each method chose, so this is the column to read when comparing methods.

Read `n_cross_array_common` alongside it. A correlation over the handful of units a block method returns carries almost no information, and should not be compared against one computed over two hundred regions.
