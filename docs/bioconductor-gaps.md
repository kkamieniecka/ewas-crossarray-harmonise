# crossarrayEWAS: what stands between the skeleton and a submission

The package in `pkg/crossarrayEWAS` is generated from `bin/ewasml.R` by
`tools/build_pkg.py` and checked by `tools/run_checks.sh`. This file records
what the two checkers actually said, measured on version 0.99.0, so the
remaining work is a list rather than a reading of the guidelines.

Toolchain used: R 4.5.3 (aarch64-apple-darwin20), roxygen2 8.1.0,
testthat 3.3.1, BiocCheck 1.46.3, plus knitr, rmarkdown, BiocStyle and pandoc
for the vignette.

The submission mechanics — the 2026 BiocContributions/R-universe process, and
the structural change it forces — are in `docs/bioconductor-submission.md`.

## Where it stands

| check | verdict |
|---|---|
| `R CMD check --no-manual` (vignette built and re-built) | **OK** — no errors, warnings or notes |
| `testthat` (`pkg/crossarrayEWAS/tests`) | 179 expectations in 54 blocks, 0 failures |
| `tests/test_pkg_identity.R` | 27 checks, 0 failures — generated bodies identical to `bin/ewasml.R`, hand-written sources present |
| `BiocCheck` | 1 error, 1 warning, 5 notes — itemised below |

`SummarizedExperiment`, `GenomicRanges`, `S4Vectors` and `IRanges` are
`Suggests`. `R/classes.R`, `tests/testthat/test-classes.R` and one vignette
section use them, each behind a `requireNamespace()` guard or a chunk-level
`eval`, so the package installs, checks and builds its vignette without them —
which is what the conda build beside the Galaxy wrappers needs. They are
installed here, so the run above exercised that code rather than skipping it;
`_R_CHECK_FORCE_SUGGESTS_=false` is still set in `tools/run_checks.sh` so the
check keeps working where they are absent. Vignette building is no longer
suppressed: `R CMD build` takes about 11 s with it and the tarball is 429 KB,
an order of magnitude larger than the code-only 31 KB and still far inside the
5 MB limit.

## The one remaining BiocCheck error

**"Unable to find your email in the Support Site: HTTP 404 Not Found."** With
the sandbox allowed to reach support.bioconductor.org, this is now a real
answer rather than a connection failure: `kkamieni@bradford.ac.uk` is not
registered there. Registration under exactly the `DESCRIPTION` address, and
subscription to the bioc-devel mailing list, are account tasks that no build
can do — and the mailing-list check reports "cannot determine" for everyone,
because it needs list-admin credentials, so it is not evidence either way.

Four findings that stood here before are closed:

- **No `vignettes` directory.** `pkg/vignettes-src/crossarrayEWAS.Rmd` is now
  written, copied into the generated tree by `tools/build_pkg.py`, and built
  by `R CMD build`. It simulates a two-array longitudinal cohort rather than
  shipping data, and its planted signals are recovered by the estimators it
  demonstrates — see `docs/bioconductor-submission.md` for the author
  paragraphs still marked in it.
- **"Use double colon for qualified imports" at `R/clusters.R`.** This was a
  false positive on `cand[st:min(...)]`, which the linter read as `st::min`.
  `bin/ewasml.R` now writes it as `seq.int(st, min(...))`, which is identical
  for `chunk >= 1L`. Pure R syntax, so `bin/ewasml.py` needed no mirror.
- **No ORCID in `Authors@R`.** `DESCRIPTION` now carries
  `person("Katarzyna", "Kamieniecka", role = c("aut", "cre"), comment =
  c(ORCID = "0009-0004-2454-5950"))` and
  `person("Krzysztof", "Poterlowicz", role = "aut", comment = c(ORCID =
  "0000-0001-6173-5674"))`. Both also appear in `CITATION.cff` and on the
  vignette. `Authors@R` is written by `tools/build_pkg.py`, so a further
  author or the `fnd` role goes there, not into the generated `DESCRIPTION`.
- **No class entry points.** `R/classes.R` now takes a `SummarizedExperiment`
  (or anything extending it, including minfi's `GenomicRatioSet`) into
  `fit_within()` and returns region and block calls as `GRanges`. It is the
  one part of the package not generated from `bin/ewasml.R` — the core stays
  base R plus limma so the stage drivers can source it without a Bioconductor
  stack — so it lives in `pkg/R-src/`, is copied in verbatim, and is excluded
  from the body comparison in `tests/test_pkg_identity.R` while its presence
  in the tree is asserted. The four functions are plain functions, not S4
  methods: `setMethod()` would need the generic's package at install time,
  which `Suggests` does not give.

## The warning

**`set.seed()` inside `stability_selection()` (`R/resample.R`).**
Bioconductor forbids a package function from resetting the caller's random
number stream, and it is right to: the function currently overwrites the
session seed on every call. The convention is to drop the `seed` argument and
document that the caller seeds, which makes the bootstrap reproducible from
outside rather than from inside. That changes the signature, so it is a
decision for the pipeline rather than a mechanical fix — the two stage drivers
pass `--seed` through to it, and `tests/test_pkg_identity.R` plus the
cross-language equivalence test both have to be updated in the same commit.
`withr::with_seed()` preserves the current behaviour without touching the
global stream, at the cost of a dependency.

## The notes worth acting on

- **Funding**: add the `fnd` role if the work is grant-supported — owner-only,
  it goes into the `Authors@R` block in `tools/build_pkg.py`, which is also
  where a further author would be added.
- **Row names out of `demean_by_group()`**: not a BiocCheck finding, but the
  package's new `linalg` tests turned it up. When the input has no row names,
  the result carries the integer group codes as row names, leaked from the
  grouped-sum step; with row names on the input they are preserved correctly.
  Numerically harmless, and the Python twin cannot show it, so it is pinned by
  a test and documented rather than fixed — the fix is one line in
  `bin/ewasml.R`.
- **Function length**: `fit_block_hsmm()` is 102 lines and `tv_denoise()` is
  55, against a recommended 50. Both are single numerical routines (a
  forward-backward pass and a fused-lasso solver); splitting them to satisfy
  the note would make them harder to compare against `bin/ewasml.py`.
- **Line length and indentation**: 6 lines over 80 characters, 497 lines
  (17%) not at a multiple of four. The package inherits the pipeline's
  two-space style; restyling would break body identity with the core. Three of
  the long lines are unavoidable: the vignette title and its
  `\VignetteIndexEntry` have to match each other, and the third is a URL in
  generated `man/`.

## Not gaps, but missing capability

The package exposes the numerical core plus the class layer over it. One
thing the guidelines care about is still absent:

- **The harmonisation stage.** `bin/01_harmonise.R` is still straight-line
  script code, so the array-combining step that gives the package its name is
  not in it. Extracting it means separating the `minfi::combineArrays()` call
  and the probe-space bookkeeping from the file layout and provenance record
  around them.

## Consuming the package instead of vendoring the core

`bin/04_dmr_ml.R` and `bin/05_blocks_hsmm.R` now prefer the installed package
and fall back to sourcing `bin/ewasml.R` when it is absent, recording which
route ran as `core_source` in the stage's run record. On the stage-05 fixture
the two routes are interchangeable: with the package installed and no
`ewasml.R` in reach, `blocks_hsmm.csv` and `openSea_cluster_effects.csv.gz`
come out byte-identical to the sourced run and all 17 compared fields of
`hsmm_params.json` agree, `core_source` being the only difference
(`crossarrayEWAS 0.99.0` against `bin/ewasml.R`).

`pinv()` is exported for this reason — the stage-04 driver uses it to
residualise the within-transformed design — which is why the package exports
18 functions rather than the 9 the original plan projected. Only
`trigamma_inv()`, `fnv1a()`, `soft_threshold()` and `dist_transitions()` stay
internal.

`conf/conda-recipe/r-crossarrayewas/` builds the package as a conda package.
It is `noarch: generic`, since the package is pure R, so one build serves every
platform. The name is `r-crossarrayewas` rather than
`bioconductor-crossarrayewas` deliberately: bioconda reserves that prefix for
packages in a Bioconductor release and generates those recipes from the
release manifest, so on acceptance this recipe should be deleted rather than
renamed.

**The Galaxy wrappers still vendor the core, and must until the package is in
a channel.** A `<requirement type="package">r-crossarrayewas</requirement>`
that resolves nowhere would break dependency resolution for every user of the
tool, so the switch — two lines in `galaxy/macros.xml`, adding the requirement
to `requirements_r_ml` and dropping `scripts/ewasml.R` from the tool
repository's sync list — is deliberately not made yet. It is unblocked by
either Bioconductor acceptance (after which bioconda's automation publishes
the package) or by building this recipe into a channel the Galaxy instance can
see. Nothing breaks in the meantime: the drivers' fallback is the vendored
core.

What did change in the wrappers is the pin set. `requirements_r_ml` now asks
for R 4.5 with limma 3.66.0, which is both the floor the package declares and
the only one of the three tools' requirement sets that resolves on Apple
silicon — verified by solving and installing it here. `ewas_harmonise` cannot
follow: `bioconductor-minfi` requires `bioconductor-illuminaio`, which has no
osx-arm64 build at any version, so that tool remains x86-only and its pins are
left alone.

## Regenerating

```sh
python3 tools/build_pkg.py                       # R/, man/ inputs, metadata
Rscript -e 'roxygen2::roxygenise("pkg/crossarrayEWAS", clean = TRUE)'
bash tools/run_checks.sh                         # build, check, BiocCheck
Rscript tests/test_pkg_identity.R                # bodies match bin/ewasml.R
```

`tools/build_pkg.py` reads `bin/ewasml.R` by default; set `EWASML_SRC` to
point it elsewhere. It rewrites `R/`, `DESCRIPTION`, `NAMESPACE`, `NEWS.md`,
the package `README.md`, `LICENSE`, `LICENSE.note` and the `testthat` files,
so hand edits to those are lost — change the generator, or the core, instead.
