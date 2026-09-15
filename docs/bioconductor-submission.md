# Submitting crossarrayEWAS to Bioconductor: the route from here

`docs/bioconductor-gaps.md` lists what the checkers say about the package.
This file is the other half: the submission mechanics, and the order the
remaining work has to happen in. Written against the process as it stands in
September 2026 — it changed in June 2026 and most tutorials still describe the
old one.

## The process now

Submissions go through an issue on
[Bioconductor/BiocContributions](https://github.com/Bioconductor/BiocContributions/issues),
not the old `Bioconductor/Contributions` repository, and builds are produced by
GitHub Actions plus R-universe rather than by Bioconductor's own build
machines. The build report that lands on the issue covers, at time of writing,
`bioc-checks` under R 4.5.3 and `R CMD check` under R 4.6.0 on
linux-devel-x86_64, macos-release-arm64 and windows-release; a new submission is
judged on the R version current for Bioconductor devel. Every platform should be
at NOTE or better, and anything left has to be justified in a comment.

## Blocking: the package needs its own repository

Precheck validation, which runs the moment the issue is opened, requires that
the DESCRIPTION `Package` field matches the **repository name, case
sensitive**, and that `DESCRIPTION` and a `vignettes` directory exist — at the
repository root. `pkg/crossarrayEWAS` inside this repository satisfies none of
that, so the first concrete step is a second public repository named exactly
`crossarrayEWAS` whose root is the generated package.

This is the same split already in use for the Galaxy suite: the pipeline stays
the source of truth, `tools/build_pkg.py` stays the generator, and a sync
script pushes the generated tree outward the way `sync-from-pipeline.sh` does
for the wrappers. Nothing about the generate-from-`bin/ewasml.R` convention
changes, and `tests/test_pkg_identity.R` keeps guarding it.

The rest of precheck already passes as the package stands: version is
`0.99.0` (incoming packages must be `x.99.y`), there is no `Remotes` and no
`Additional_repositories` field, no file is anywhere near the 5 MB ceiling (the
whole source tarball is 31 KB), and the repository uses no Git LFS.

## Blocking: the vignette

`pkg/crossarrayEWAS/vignettes/` exists but is empty, and `DESCRIPTION` carries
no `VignetteBuilder`. Precheck only tests that the directory exists, so this
would pass validation and then fail the source build — and BiocCheck's
missing-vignette error stands either way. A vignette is also what the reviewer
reads first.

It needs a committable example object, since no study data may go into the
package. The synthetic fixtures from `tests/gen_equivalence_fixtures.py` are
the seed. Adding it means `VignetteBuilder: knitr` in `DESCRIPTION` and the
`knitr`/`rmarkdown`/`BiocStyle` suggests that are already declared.

## Blocking: things no build can supply

These are account and identity tasks for the maintainer, and the submission
cannot be completed without them:

- A real `Authors@R` email — still `ENTER.YOUR.ADDRESS@example.org`. It has to
  reach the maintainer for as long as the package is in Bioconductor.
- An account at support.bioconductor.org registered under that same address,
  which is what BiocCheck's third error is testing.
- A subscription to the bioc-devel mailing list.
- An ORCID in `comment = c(ORCID = ...)`, and the `fnd` role if the work is
  grant-supported.
- A disclosure of AI-assisted code. Bioconductor's own guide (pkgrevdocs,
  `ai-policy-third-party.Rmd`) requires that non-trivially AI-generated or
  copied code be raised in the submission issue with provenance, cited in the
  code, and redistributable under the package licence, with the submitter
  responsible for the result. Much of `bin/ewasml.R` came from assisted
  sessions, so this applies directly. The issue template itself asks only for
  the repository URL; the disclosure goes in as a comment.

## Wanted before review, not by precheck

- **Class entry points.** Nothing accepts a `SummarizedExperiment` or
  `GenomicRatioSet`, and `call_regions()`/`call_blocks()` return data frames
  rather than `GRanges`. Reviewers ask for this specifically.
- **The `set.seed()` warning** in `stability_selection()`. A package function
  must not reset the caller's stream. Changing it changes the signature, so the
  two stage drivers, the identity test and the cross-language equivalence test
  move in the same commit.
- **The two cosmetic core fixes**, both in `bin/ewasml.R` and mirrored into
  `bin/ewasml.py`: `st:min(...)` written as `seq(...)` to silence the
  false-positive double-colon error, and the four `1:n` occurrences in the
  transition-matrix construction.
- **A check under R 4.6.** Everything here has been checked under 4.5.3; three
  of the four report platforms run 4.6.0. The
  [Bioconductor R-universe GitHub Action](https://docs.r-universe.dev/bioconductor/#debugging-the-ci)
  mimics the submission build and can be added to the package repository to see
  those reports before submitting rather than after.

## Then the issue

1. Open an issue on `Bioconductor/BiocContributions` using
   `new_submission_template`, unmodified, titled with the package name; the
   body is just the repository URL.
2. Precheck runs automatically. On failure the issue is closed with a comment
   saying what failed, and a fresh issue repeats the validation.
3. Comment exactly `/accept-policies`.
4. The repository is cloned into a staging organisation and registered with
   R-universe; the issue comments the new remote. **Subsequent builds are
   triggered by pushes to that remote, not to the original GitHub repository**,
   and only by a valid version bump — advance only `z`: 0.99.1, 0.99.2. Reports
   can take up to 24 hours.
5. Clear the errors, justify whatever remains, and a reviewer is assigned;
   review typically happens within three weeks of a clean build, usually over
   several rounds.
6. On acceptance the package is cloned to the canonical Bioconductor location,
   added to the devel manifest, and a BiocCredentials account manages the SSH
   keys for push access.

The sync script matters at step 4: once the staging remote is live, the
generated tree has to be pushed there, so it should take the target remote as
an argument rather than hard-coding the GitHub one.

## What acceptance unblocks downstream

bioconda generates `bioconductor-*` recipes from the Bioconductor release
manifest, so acceptance produces `bioconductor-crossarrayewas` without any
work here — at which point `conf/conda-recipe/r-crossarrayewas/` should be
deleted rather than renamed, and the Galaxy wrappers can stop vendoring the
core: add the requirement to `requirements_r_ml` in `galaxy/macros.xml`, drop
`scripts/ewasml.R` from the tool repository's sync list. The drivers already
prefer the installed package, so they need no change.

Acceptance lands the package in Bioconductor devel, which becomes a release at
the next cycle; releases are twice a year, in April and October. Nothing in the
Galaxy deployment waits on this — the drivers' fallback is the vendored core.
