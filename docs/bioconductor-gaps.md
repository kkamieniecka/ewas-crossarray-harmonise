# crossarrayEWAS: what stands between the skeleton and a submission

The package in `pkg/crossarrayEWAS` is generated from `bin/ewasml.R` by
`tools/build_pkg.py` and checked by `tools/run_checks.sh`. This file records
what the two checkers actually said, measured on version 0.99.0, so the
remaining work is a list rather than a reading of the guidelines.

Toolchain used: R 4.5.3 (aarch64-apple-darwin20), roxygen2 8.1.0,
testthat 3.3.1, BiocCheck 1.46.3.

## Where it stands

| check | verdict |
|---|---|
| `R CMD check --no-build-vignettes --no-manual` | **OK** — no errors, warnings or notes |
| `testthat` (`pkg/crossarrayEWAS/tests`) | 126 expectations, 0 failures |
| `tests/test_pkg_identity.R` | 26 checks, 0 failures — package bodies identical to `bin/ewasml.R` |
| `BiocCheck` | 3 errors, 1 warning, 8 notes — itemised below |

`R CMD check` was run with `_R_CHECK_FORCE_SUGGESTS_=false`, because the
suggested packages that only the unwritten vignette and class entry points
need (`BiocStyle`, `rmarkdown`, `SummarizedExperiment`, `GenomicRanges`) are
not installed on the machine that ran it. Nothing in `R/` or `tests/` imports
them.

## The three BiocCheck errors

**1. No `vignettes` directory.** Bioconductor will not review a package
without a vignette, and it has to be a narrative that runs, not a manual page
index. This is the largest remaining piece of work and it needs a small
committable example object, because no study data can go in the package: the
synthetic fixtures in `tests/gen_equivalence_fixtures.py` are the natural
seed, wrapped as a `SummarizedExperiment` in `inst/extdata` or built in the
vignette itself.

**2. "Use double colon for qualified imports" at `R/clusters.R` line 67.**
This is a false positive. The flagged column is the `:` in

```r
sl <- cand[st:min(st + chunk - 1L, length(cand))]
```

which the linter reads as a single-colon namespace access (`pkg:foo()`).
Writing it as `seq(st, min(st + chunk - 1L, length(cand)))` silences the
check and is exactly equivalent here, since `chunk >= 1L` makes the sequence
increasing. The edit belongs in `bin/ewasml.R`, not in `R/clusters.R` —
`tests/test_pkg_identity.R` fails if the package copy diverges — and once made
it should be mirrored into `bin/ewasml.py` to keep the two cores in step.

**3. "Unable to find your email in the Support Site."** Two separate things
are behind this. The maintainer must be registered at
support.bioconductor.org under the address in `DESCRIPTION`, and subscribed to
the bioc-devel mailing list; both are account tasks that cannot be done from a
build. The check also could not reach the site from this sandbox, so even a
registered address would have reported an error here. The `Authors@R` email is
still the placeholder `ENTER.YOUR.ADDRESS@example.org` — it must be a real
address that reaches the maintainer for as long as the package is in
Bioconductor.

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

- **`biocViews`**: BiocCheck suggests adding `ChipOnChip`. The current terms
  are `DNAMethylation`, `DifferentialMethylation`, `Epigenetics`,
  `MethylationArray`, `Regression`, `Software`.
- **ORCID**: add `comment = c(ORCID = "…")` to `Authors@R`.
- **Funding**: add the `fnd` role if the work is grant-supported.
- **`1:n` idiom**: 4 occurrences, all in `R/blocks.R` (two each on lines 25
  and 220, inside the transition-matrix construction). Both are guarded by a
  fixed three-state dimension, so neither is a live bug, but the same argument
  as the double-colon error applies — fix in `bin/ewasml.R` and regenerate, never in
  the package copy.
- **Function length**: `fit_block_hsmm()` is 102 lines and `tv_denoise()` is
  55, against a recommended 50. Both are single numerical routines (a
  forward-backward pass and a fused-lasso solver); splitting them to satisfy
  the note would make them harder to compare against `bin/ewasml.py`.
- **Line length and indentation**: 4 lines over 80 characters, 354 lines not
  at a multiple of four. The package inherits the pipeline's two-space style;
  restyling would break body identity with the core.

## Not gaps, but missing capability

The package currently exposes the numerical core only. Two things the
guidelines care about are still absent by choice:

- **Class entry points.** Nothing accepts a `SummarizedExperiment` or
  `GenomicRatioSet`, and `call_regions()`/`call_blocks()` return data frames
  rather than `GRanges`. Reviewers ask for this specifically; it is the second
  substantive piece of work after the vignette.
- **The harmonisation stage.** `bin/01_harmonise.R` is still straight-line
  script code, so the array-combining step that gives the package its name is
  not in it. Extracting it means separating the `minfi::combineArrays()` call
  and the probe-space bookkeeping from the file layout and provenance record
  around them.

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
