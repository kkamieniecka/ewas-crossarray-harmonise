# Galaxy tool test fixtures

The Galaxy wrappers in `galaxy/` declare `<tests>` but no fixtures are
committed here yet, because a useful fixture for these tools is not small: a
harmonisation test needs real IDAT pairs from both array types, and a region
test needs a matrix with enough subjects for a within-subject permutation to
be meaningful.

The intended fixtures, in the order they are worth adding:

1. `mval_tiny.f64` + `pheno_tiny.csv` + `anno_tiny.csv` — a simulated matrix
   of a few thousand probes on one chromosome arm, 8 subjects x 3 visits,
   with a planted piecewise-constant region. This exercises
   `ewas_dmr_ml.xml` and `ewas_blocks_hsmm.xml` end to end in seconds.
   `tests/test_ewasml.py` already builds equivalent data in memory; the
   fixture is that generator's output written to disk with a fixed seed.
2. A two-sample-per-array IDAT set for `ewas_harmonise.xml`. GEO's smallest
   450K/EPIC pairs are still ~8 MB each, so these belong in a
   `test-data` release asset or a Zenodo deposit referenced by URL, not in
   git history.

Until then, `python tests/test_ewasml.py` is the real test suite and runs in
CI on every push.
