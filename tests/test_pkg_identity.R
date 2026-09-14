#!/usr/bin/env Rscript
## Does pkg/crossarrayEWAS still hold the same functions as bin/ewasml.R?
##
## The package is generated from the core by tools/build_pkg.py, which copies
## every function body verbatim and rewrites only the comment above it. That
## makes the package documentation and namespace safe to regenerate, but only
## as long as nobody edits R/*.R by hand. This test compares the two sources
## function by function and fails if any body has drifted.
##
## It deliberately does not install the package: sourcing the R/ files is
## enough to compare definitions, and keeps the test runnable without the
## Bioconductor toolchain.
##
## Usage:  Rscript tests/test_pkg_identity.R

suppressPackageStartupMessages({
  library(Matrix)
  library(jsonlite)
})

here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=",
        commandArgs(FALSE), value = TRUE)[1])))
root <- normalizePath(file.path(here, ".."))
core <- file.path(root, "bin", "ewasml.R")
pkg_r <- list.files(file.path(root, "pkg", "crossarrayEWAS", "R"),
                    pattern = "[.]R$", full.names = TRUE)
pkg_r <- pkg_r[!grepl("-package[.]R$", pkg_r)]

if (!length(pkg_r)) stop("no package sources found; run tools/build_pkg.py")

e_core <- new.env()
e_pkg <- new.env()
sys.source(core, envir = e_core)
for (f in pkg_r) sys.source(f, envir = e_pkg)

obj <- function(e) Filter(function(n) is.function(get(n, e)), ls(e))
n_core <- sort(obj(e_core))
n_pkg <- sort(obj(e_pkg))

fails <- 0L
report <- function(ok, msg) {
  cat(sprintf("%-5s %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

report(identical(n_core, n_pkg), sprintf(
  "same function set (core %d, package %d)%s", length(n_core), length(n_pkg),
  if (identical(n_core, n_pkg)) "" else sprintf(
    "; missing: %s; extra: %s",
    paste(setdiff(n_core, n_pkg), collapse = " "),
    paste(setdiff(n_pkg, n_core), collapse = " "))))

for (n in intersect(n_core, n_pkg))
  report(identical(deparse(get(n, e_core)), deparse(get(n, e_pkg))),
         sprintf("%s() body unchanged", n))

## Constants used by the hash-based fold assignment travel with the functions.
for (k in c("FNV_OFFSET", "FNV_PRIME"))
  report(exists(k, e_pkg, inherits = FALSE) &&
           identical(get(k, e_core), get(k, e_pkg)),
         sprintf("%s carried over", k))

## Everything the exported functions promise must actually be exported.
ns <- readLines(file.path(root, "pkg", "crossarrayEWAS", "NAMESPACE"))
exported <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
report(all(exported %in% n_pkg),
       sprintf("all %d exports exist (%s)", length(exported),
               paste(setdiff(exported, n_pkg), collapse = " ")))

cat(sprintf("\n%d checks, %d failed\n", length(n_core) + 4L, fails))
quit(status = if (fails) 1L else 0L)
