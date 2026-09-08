# Installation

Three routes. Use containers if you can; the conda route is for development,
and the Apple Silicon route exists because Bioconda has no builds for that
platform.

## 1. Containers (recommended)

```bash
nextflow run . -profile docker      # or -profile singularity
```

Nothing else to install beyond Nextflow and a container runtime. This is the
route to use for anything you intend to publish, because the software stack is
pinned by image digest rather than by solver outcome.

## 2. Conda, on Linux or Intel macOS

```bash
conda env create -f conf/env-r.yml
conda env create -f conf/env-py.yml
nextflow run . -profile conda
```

## 3. Apple Silicon (osx-arm64)

Bioconda publishes no `bioconductor-*` builds for `osx-arm64`, so route 2
fails at the solver. The Bioconductor project's own arm64 binaries do not
help either: they are built against the CRAN R *framework*, so their compiled
objects resolve `libR` at a path holding a different R version than a conda R.

Build the CRAN layer from conda-forge and compile only the Bioconductor layer
from source:

```bash
conda create -n ewas -c conda-forge \
    r-base=4.4 r-biocmanager r-data.table r-optparse r-jsonlite r-nlme \
    c-compiler cxx-compiler fortran-compiler make
conda activate ewas
Rscript conf/install_bioc.R            # writes into ./.r-libs/ewas
export R_LIBS_USER="$PWD/.r-libs/ewas"
```

Expect 30–60 minutes for the compile. `conf/install_bioc.R`:

- pins Bioconductor to the release matching the environment's R version rather
  than the newest release;
- installs into a **workspace-local** library, because a conda environment's
  own `site-library` is read-only under some managed installs;
- builds `preprocessCore` with `--disable-threading`. Its OpenMP path
  deadlocks on macOS builds that link a different `libomp` than the compiler
  used, and the failure presents as a hang rather than an error;
- guards against `parallel::detectCores()` returning `NA`, which otherwise
  propagates into `MAKEFLAGS=-jNA` and stops the build;
- ends with a presence-and-version audit of every required package, so a
  partial install fails loudly here instead of three stages later.

Set `R_LIBS_USER` in every shell that runs the R stages, or add it to
`nextflow.config` under the `conda` profile's `env` block.

## Python layer

The numerical core (`bin/ewasml.py`, and the two Python stages) needs only
numpy, scipy and pandas, and imports nothing from R or Bioconductor. That is
deliberate: it means the estimators are unit-testable without the
Bioconductor stack.

```bash
python tests/test_ewasml.py     # 17 property checks, seconds
```

## Verifying an install

```bash
nextflow run . -profile test,conda --sheet <sheet> --idat_dir <dir> \
    --probe_map <map>
```

The `test` profile runs the identical graph with 10 permutations, 5 bootstrap
resamples and the legacy baseline switched off. It reduces the resampling
counts only — it does not subset the input arrays, so it still needs the full
IDAT set and the runtime is dominated by the per-probe stages. Measured on
GSE237561 (126 arrays, 8 CPU cores, no container engine): harmonisation
6 min 16 s, per-probe model 17 min 12 s, region finder 14 min 52 s, block
finder 16 s, comparison 10 s, about 40 min wall in total.

With `--run_baseline false` the comparison stage still runs and writes a table
containing the two replacement methods only; the legacy rows appear when the
baseline stage is on. Every process in the graph therefore executes under the
test profile.

The test profile also lowers the `r_heavy` memory request. That request caps
what the executor will admit, not what a stage uses: the default 32 GB request
is refused outright by the local executor on a 16 GB machine.

### Running from conda environments without a container engine

The processes call bare `Rscript` and `python`, so they inherit the launching
shell's `PATH`. Two things bite when the environments are activated by hand
rather than through `-profile conda`:

* Nextflow's own environment ships a Python without numpy/pandas. Put the
  analysis environments *before* it on `PATH`, or the `py_*` processes get the
  wrong interpreter.
* `openjdk` from conda-forge keeps the runtime outside the environment's
  `bin/`, so `JAVA_HOME` has to be set explicitly for the launcher to find it.

`06_compare.py` writes its Markdown table through `tabulate` when it is
installed (it is declared in `conf/env-py.yml`) and falls back to a plain pipe
table when it is not, so an ad-hoc environment does not fail at the last
stage.

## What each stage needs

`04_dmr_ml.R` (the default region finder) needs only `r-base`, `Matrix`,
`optparse`, `jsonlite` and `limma` — not minfi or bumphunter, which are
required by the harmonisation and legacy-baseline stages. That matters if you
want to run only the region finder: the light dependency set installs from
conda-forge plus one Bioconductor package, and `--var-method mom` drops limma
too, at the cost of losing the default moderation.

`05_blocks_hsmm.R` (the default block finder) needs the same light set minus
`Matrix`: `r-base`, `optparse`, `jsonlite`, and `limma` only for the default
variance moderation. It additionally needs the cross-array probe map that
stage 01 writes, since it restricts to open-sea probes carried by both arrays.

Memory: the region finder holds a few copies of the M-value matrix. Measured
peak RSS is 1.78 GiB for R and 0.96 GiB for Python at 100k probes x 126
samples, scaling linearly in probes — so roughly 7.5 GiB (R) or 4.1 GiB
(Python) on a full EPIC panel. The `r_heavy` label requests 24 GB, which
covers it with headroom for the harmonisation stage.
