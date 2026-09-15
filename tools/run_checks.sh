#!/usr/bin/env bash
# Document, install, test and BiocCheck the package.
# Vignettes are built and checked: knitr, rmarkdown, BiocStyle and pandoc
# have to be present, which is what the submission build machines have.
set -uo pipefail
export PATH="$HOME/.claude-science/conda/envs/ewas/bin:$PATH"
# The development library (BiocCheck, roxygen2, testthat) is scratch and
# lives outside the repository; override R_LIBS_USER to point at yours.
export R_LIBS_USER="${R_LIBS_USER:-$PWD/.r-libs/pkgdev}"
PKG=pkg/crossarrayEWAS
# SummarizedExperiment and GenomicRanges are Suggests only, and neither has
# an osx-arm64 bioconda build; the submission build machines have them, so
# locally the check runs without requiring suggests.
export _R_CHECK_FORCE_SUGGESTS_=false
mkdir -p check-out

step() { printf '\n===== %s\n' "$1"; }

step "roxygen2: man/ and NAMESPACE"
Rscript -e 'roxygen2::roxygenise("pkg/crossarrayEWAS", clean = TRUE)' 2>&1 | tail -20
ls "$PKG/man" 2>/dev/null | tr '\n' ' '; echo

step "R CMD build"
R CMD build "$PKG" 2>&1 | tail -6
TARBALL=$(ls -t crossarrayEWAS_*.tar.gz 2>/dev/null | head -1)
echo "tarball: $TARBALL $(du -h "$TARBALL" 2>/dev/null | cut -f1)"

step "R CMD check"
/usr/bin/time -p R CMD check --no-manual \
    -o check-out "$TARBALL" 2>&1 | grep -E "^\*|WARNING|NOTE|ERROR|real" | tail -40

step "BiocCheck"
Rscript -e 'BiocCheck::BiocCheck("'"$TARBALL"'", `no-check-bioc-views` = FALSE)' 2>&1 | tail -80
