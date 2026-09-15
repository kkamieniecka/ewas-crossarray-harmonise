#!/usr/bin/env bash
# Publish the generated package tree to a standalone crossarrayEWAS
# repository, which is what a Bioconductor submission needs: precheck requires
# the DESCRIPTION Package field to match the repository name case-sensitively
# with the package at the repository ROOT, and pkg/crossarrayEWAS inside this
# repository cannot satisfy that. This repository stays the source of truth;
# the package repository is a published artefact of it, exactly as
# ewas-crossarray-harmonise is for the Galaxy wrappers.
#
#   tools/sync-to-pkg-repo.sh /path/to/crossarrayEWAS [remote]
#
# The target is a checkout rather than a URL, so the same script serves both
# destinations the submission needs: your own GitHub repository before the
# issue is opened, and the Bioconductor staging remote afterwards. Give the
# checkout a second remote and pass its name -- after policy acceptance only
# pushes to the staging remote trigger builds, and only a z-level version bump
# produces a new report.
#
# Nothing is committed or pushed here: review the diff, bump Version in
# tools/build_pkg.py when the push is meant to trigger a build, and commit in
# the package checkout with the upstream commit this prints.
set -euo pipefail

dst=${1:?usage: sync-to-pkg-repo.sh /path/to/crossarrayEWAS [remote]}
remote=${2:-origin}
here=$(cd "$(dirname "$0")/.." && pwd)

[ -d "$dst/.git" ] || { echo "not a git checkout: $dst" >&2; exit 2; }
git -C "$dst" remote get-url "$remote" >/dev/null 2>&1 ||
    { echo "no remote '$remote' in $dst" >&2; exit 2; }
[ "$(basename "$(cd "$dst" && pwd)")" = crossarrayEWAS ] ||
    { echo "checkout must be named crossarrayEWAS (precheck is case sensitive)" >&2
      exit 2; }

# Regenerate first: pkg/crossarrayEWAS/R/ is generated from bin/ewasml.R and
# is never edited in place, so syncing without regenerating would publish a
# stale core.
cd "$here"
python3 tools/build_pkg.py
# roxygen2 reports every .Rd it writes on stderr; set -e catches real failures
Rscript -e 'roxygen2::roxygenise("pkg/crossarrayEWAS", clean = TRUE)' >/dev/null 2>&1
Rscript tests/test_pkg_identity.R >/dev/null

# --delete so a function removed from the core disappears downstream too;
# .git and the package repository's own workflows are not ours to manage.
rsync -a --delete --exclude .git --exclude .github \
    "$here/pkg/crossarrayEWAS/" "$dst/"

# The R-universe test workflow is installed once and then left alone: it is
# the package repository's own CI, not a generated file, and re-syncing must
# not revert edits made to it there.
wf=$dst/.github/workflows/r-universe.yml
if [ ! -f "$wf" ]; then
    mkdir -p "$(dirname "$wf")"
    cp "$here/pkg/pkg-repo-template/.github/workflows/r-universe.yml" "$wf"
    echo "installed .github/workflows/r-universe.yml (was missing)"
fi

up=$(git -C "$here" rev-parse --short HEAD 2>/dev/null || echo unknown)
ver=$(awk '/^Version:/ {print $2}' "$dst/DESCRIPTION")
echo
echo "synced crossarrayEWAS $ver from pipeline commit $up"
echo "remote: $remote -> $(git -C "$dst" remote get-url "$remote")"
echo "review: git -C $dst status --short"
