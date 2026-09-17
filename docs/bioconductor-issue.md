# The submission issue: what to paste, in what order

Not part of the package. Companion to `docs/bioconductor-submission.md`, which
explains the process; this file is only the text to paste and the order to
paste it in. Open the issue on
[Bioconductor/BiocContributions](https://github.com/Bioconductor/BiocContributions/issues)
with the `new_submission_template`, unmodified.

## 1. Issue title

    crossarrayEWAS

## 2. Issue body

The template asks for the repository URL and nothing else. Keep the template's
own checkboxes as they come and put this as the body text:

    https://github.com/bioinfbrad/crossarrayEWAS

That repository's root is the package: `DESCRIPTION` with `Package:
crossarrayEWAS` matching the repository name case-sensitively, `vignettes/`
present, version 0.99.0, no `Remotes` and no `Additional_repositories`, source
tarball 427 KB, no Git LFS. Precheck runs the moment the issue is opened and
checks exactly those things.

## 3. First comment, after precheck passes

Exactly this, on its own:

    /accept-policies

## 4. Second comment: the AI-assistance disclosure

Paste `docs/ai-disclosure-comment.md` from the `---` divider onwards, in your
own words. It is required by the Bioconductor policy on AI-generated and
third-party code (pkgrevdocs, `ai-policy-third-party.Rmd`) and the issue
template does not ask for it, so it has to be volunteered as a comment.

## 5. Optional third comment: the outstanding notes

Only if a reviewer asks why BiocCheck is not at zero notes. The package is at
**0 errors, 0 warnings, 5 notes**, and the justifications are:

    The five remaining BiocCheck notes are deliberate rather than outstanding:

    - No `fnd` role in `Authors@R`. The work is PhD research and has no grant
      funding, so there is no funder to name.
    - `fit_block_hsmm()` is 102 lines and `tv_denoise()` is 55, against the
      recommended 50. Each is a single numerical routine — a forward-backward
      pass and a fused-lasso solver — and the package's core is held
      byte-identical to a Python implementation of the same functions by a
      cross-language equivalence test. Splitting them would make that
      comparison harder to maintain without making either routine clearer.
    - Six lines exceed 80 characters. Three are unavoidable: the vignette
      title and its `\VignetteIndexEntry` have to match each other, and the
      third is a URL in generated `man/`.
    - 17% of lines are not indented at a multiple of four. The package is
      generated from the pipeline's core, which uses two-space indentation;
      restyling the generated tree would break the body-identity test that
      guarantees the package and the validated core are the same code.
    - Row names out of `demean_by_group()` when the input has none. This is
      not a BiocCheck finding but is recorded in the package's own notes: the
      integer group codes leak through as row names in that case. It is
      numerically harmless, pinned by a test, and documented.

## 6. After the issue is open, builds move

The repository is cloned into a staging organisation and registered with
R-universe; the issue comments the new remote. From then on **only pushes to
that remote trigger builds**, and only with a valid version bump: advance `z`
only, 0.99.1 then 0.99.2. Reports can take up to 24 hours.

Add the staging remote to the package checkout and pass its name to the sync
script, which takes it as the second argument:

    tools/sync-to-pkg-repo.sh /path/to/crossarrayEWAS <staging-remote>

Bump `Version` in `tools/build_pkg.py` — not in the generated `DESCRIPTION` —
so the next regeneration carries it.
