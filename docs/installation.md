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

The `test` profile runs the identical graph with 10 permutations and the legacy
baseline switched off, which exercises every process and file contract in a
few minutes rather than hours.
