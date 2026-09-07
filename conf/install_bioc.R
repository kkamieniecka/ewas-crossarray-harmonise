#!/usr/bin/env Rscript
# Source-install the Bioconductor layer of the pipeline.
#
# Why source: Bioconda publishes no osx-arm64 builds of the bioconductor-*
# packages, and the Bioconductor project's own macOS arm64 binaries are built
# against the CRAN R framework, not this conda R, so their compiled objects
# resolve libR at a path that holds a different R version. The CRAN layer
# therefore comes from conda-forge (prebuilt) and only the Bioconductor layer
# is compiled locally.
#
# The conda env's own R library is mounted read-only, so everything installs
# into a writable workspace library instead.
# detectCores() returns NA in some sandboxes; an NA Ncpus makes
# install.packages() fail inside its own parallel-install branch.
NC <- suppressWarnings(as.integer(parallel::detectCores()))
if (!isTRUE(is.finite(NC)) || NC < 1L) NC <- 4L
NC <- max(1L, NC - 1L)
options(repos = c(CRAN = "https://cloud.r-project.org"),
        Ncpus = NC, timeout = 3600)
Sys.setenv(MAKEFLAGS = paste0("-j", NC))

LIB <- Sys.getenv("EWAS_R_LIB", file.path(getwd(), ".r-libs", "ewas"))
dir.create(LIB, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(LIB, .libPaths()))

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
msg("install library:", LIB)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = LIB)
BIOC <- as.character(BiocManager::version())
msg("Bioconductor version:", BIOC)

# preprocessCore's OpenMP threading is known to deadlock on macOS builds that
# link a different libomp than the compiler; minfi only needs it single
# threaded, so disable threading explicitly rather than debug it later.
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  msg("installing preprocessCore (threading disabled)")
  BiocManager::install("preprocessCore", type = "source", ask = FALSE,
                       update = FALSE, lib = LIB,
                       configure.args = "--disable-threading")
}

core <- c("limma", "illuminaio", "bumphunter", "minfi")
anno <- c("IlluminaHumanMethylation450kmanifest",
          "IlluminaHumanMethylationEPICmanifest",
          "IlluminaHumanMethylation450kanno.ilmn12.hg19",
          "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

for (p in c(core, anno)) {
  if (requireNamespace(p, quietly = TRUE)) { msg("have", p); next }
  msg("installing", p)
  BiocManager::install(p, type = "source", ask = FALSE, update = FALSE, lib = LIB)
  if (!requireNamespace(p, quietly = TRUE)) msg("!! FAILED", p)
}

msg("--- final check ---")
for (p in c(core, anno)) {
  ok <- requireNamespace(p, quietly = TRUE)
  v <- if (ok) as.character(utils::packageVersion(p)) else "-"
  cat(sprintf("%-46s %-6s %s\n", p, ifelse(ok, "OK", "MISSING"), v))
}
