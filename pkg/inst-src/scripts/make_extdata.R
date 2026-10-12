#!/usr/bin/env Rscript
# Build inst/extdata/GSE237561_chr21_subset.rds, the real-data subset the
# vignette analyses.
#
# Source: GEO series GSE237561 (Gillespie et al., 2024, Transl. Psychiatry 14,
# 390, doi:10.1038/s41398-024-03102-8), a longitudinal clozapine cohort of 126
# arrays from 38 subjects, 92 on HumanMethylation450 and 34 on MethylationEPIC,
# with no subject spanning the two generations. The depositors' normalised beta
# matrix is used rather than the IDATs, so this script needs no preprocessing
# decisions of its own.
#
# The subset is chromosome 21 probes present on both designs, which is what the
# harmonisation stage keeps, cut to the first 2,000 by position so the shipped
# object stays small. Nothing is filtered on effect size or significance: a
# subset chosen by outcome would make the vignette's region calls a foregone
# conclusion.
#
#   Rscript make_extdata.R [sample_sheet.csv] [out_dir]
#
# Set GSE237561_BETA to a local copy of the beta matrix to skip the download.
#
# The sample sheet is the pipeline's own (one row per array, with Sentrix_ID,
# Sentrix_Position, Subject_ID, visit, days_on_clozapine and Array_Type); it is
# written by stage 00 from the GEO series metadata. Downloads about 250 MB and
# takes a few minutes, which is why the result ships with the package.

args <- commandArgs(TRUE)
sheet_csv <- if (length(args) >= 1L) args[1] else "results/sample_sheet.csv"
out_dir <- if (length(args) >= 2L) args[2] else "inst/extdata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# attached rather than namespace-qualified: minfi resolves the annotation
# object off the search path
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

N_PROBES <- 2000L
BETA_URL <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/",
                   "GSE237561/suppl/GSE237561_ProcessedData_",
                   "LongitudinalClozapineNormalisedBetaMatrix.csv.gz")

## ---- probes on chromosome 21, in position order --------------------------
ann <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann21 <- as.data.frame(ann[ann$chr == "chr21", c("chr", "pos")])
ann21 <- ann21[order(ann21$pos), ]

## ---- read only those rows of the beta matrix -----------------------------
# the full matrix is ~800k probes by 126 arrays, so it is read in chunks
# a cached copy can be given in GSE237561_BETA, so a rebuild does not
# re-download a quarter of a gigabyte
tmp <- Sys.getenv("GSE237561_BETA", "")
if (!nzchar(tmp) || !file.exists(tmp)) {
  tmp <- tempfile(fileext = ".csv.gz")
  utils::download.file(BETA_URL, tmp, mode = "wb")
}
con <- gzfile(tmp, "rt")
header <- scan(con, what = "", nlines = 1L, sep = ",", quiet = TRUE)
samples <- header[-1]
rows <- list()
repeat {
  chunk <- try(utils::read.csv(con, header = FALSE, nrows = 20000L,
                               row.names = 1L, col.names = header),
               silent = TRUE)
  if (inherits(chunk, "try-error") || !nrow(chunk)) break
  hit <- rownames(chunk) %in% rownames(ann21)
  if (any(hit)) rows[[length(rows) + 1L]] <- chunk[hit, , drop = FALSE]
}
close(con)
beta <- as.matrix(do.call(rbind, rows))
colnames(beta) <- samples

probes <- utils::head(intersect(rownames(ann21), rownames(beta)), N_PROBES)
beta <- beta[probes, , drop = FALSE]

## ---- design, from the pipeline's sample sheet ----------------------------
ss <- utils::read.csv(sheet_csv, stringsAsFactors = FALSE)
rownames(ss) <- paste(ss$Sentrix_ID, ss$Sentrix_Position, sep = "_")
ss <- ss[colnames(beta), , drop = FALSE]
cd <- data.frame(
  subject = as.character(ss$Subject_ID),
  visit = as.character(ss$visit),
  # the exposure the published analysis used: days on treatment per 100 days,
  # so an effect reads as M-value change per 100 days
  exposure = as.numeric(ss$days_on_clozapine) / 100,
  array_type = as.character(ss$Array_Type),
  row.names = colnames(beta),
  stringsAsFactors = FALSE)
stopifnot(!anyNA(cd$exposure), nrow(cd) == ncol(beta),
          length(unique(cd$subject)) > 1L)

## ---- assemble a GenomicRatioSet ------------------------------------------
gr <- GenomicRanges::GRanges(
  ann21[probes, "chr"],
  IRanges::IRanges(ann21[probes, "pos"], width = 2L))
names(gr) <- probes
grset <- minfi::GenomicRatioSet(
  gr = gr, Beta = beta, M = NULL, CN = NULL,
  colData = S4Vectors::DataFrame(cd),
  annotation = c(array = "IlluminaHumanMethylation450k",
                 annotation = "ilmn12.hg19"),
  preprocessMethod = c(rg.norm = "GEO GSE237561 normalised beta matrix"))

out <- file.path(out_dir, "GSE237561_chr21_subset.rds")
saveRDS(grset, out, compress = "xz")
cat(sprintf("%s: %d probes x %d samples, %.2f MB\n", out, nrow(grset),
            ncol(grset), file.size(out) / 1e6))
