---
name: Bug report
about: A stage failed, or produced something that looks wrong
labels: bug
---

**Stage**
Which script or Nextflow process (e.g. `04_dmr_ml.py`, `PROBE_MODEL`).

**Command**
The exact command or `nextflow run` invocation.

**Run record**
Attach the stage's run record if one was written (`run_config.json`,
`hsmm_params.json`, `baseline_summary.json`) and the tail of the log.

**Design**
Number of subjects, visits per subject, and how `Array_Type` maps onto cohort
and Sentrix chip. Most statistical surprises in this pipeline trace back to
the design rather than to the code.

**Environment**
OS and architecture, `nextflow -version`, and whether you used the `conda`,
`docker` or `singularity` profile. On Apple Silicon, say whether the
Bioconductor layer came from `conf/install_bioc.R`.
