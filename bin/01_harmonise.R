#!/usr/bin/env Rscript
# 01_harmonise.R -- cross-array (450K + EPIC) read, QC and harmonisation.
#
# Replaces the EWASGalaxy "load + preprocess" stage, which assumed a single
# array type. Key design points:
#   * IDATs are read per array type with the array's own manifest, and the two
#     objects are merged only through minfi::combineArrays(), never by rbind-ing
#     beta matrices. This matters because bead AddressA_ID is different for
#     100% of the 454,181 CpGs shared by the GPL13534 and GPL21145 manifests
#     (see results/manifest_concordance.csv), so any address-level join yields
#     an empty intersection; the join must be on CpG name, with the manifest
#     used to resolve addresses per array.
#   * That same comparison shows the shared CpGs are chemically identical
#     across the two arrays -- Infinium design type, both probe sequences,
#     colour channel and hg19 position agree for 454,181/454,181 probes. So
#     name-level intersection does NOT mix type-I and type-II chemistry, and
#     no probe re-typing correction is required. Type-I/II scale correction
#     (BMIQ) is therefore optional here and off by default; what harmonisation
#     really costs is probe CONTENT (31,396 450K-only and 413,745 EPIC-only
#     CpGs are unusable), not chemistry.
#   * Noob dye-bias/background correction is applied PER ARRAY (it is a
#     within-array control-probe method), then arrays are combined, then
#     BMIQ puts type-I and type-II probes on a common scale.
#   * Array type is perfectly confounded with cohort and chip in this design,
#     so it is NOT removed by ComBat here. It is handled downstream as a
#     between-subject nuisance that cancels in within-subject contrasts.
#
# Usage:
#   Rscript 01_harmonise.R --sheet sample_sheet.csv --idat_dir data/idat \
#     --out_dir results/harmonised [--drop_sex TRUE] [--detp 0.01]

suppressPackageStartupMessages({
  library(optparse); library(minfi); library(data.table)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sheet",    type = "character"),
  make_option("--idat_dir", type = "character"),
  make_option("--out_dir",  type = "character", default = "results/harmonised"),
  make_option("--detp",     type = "double",    default = 0.01),
  make_option("--detp_frac",type = "double",    default = 0.05),
  make_option("--drop_sex", type = "logical",   default = TRUE),
  # off by default: every probe shared by 450K and EPIC has the same Infinium
  # design type, so harmonisation introduces no type-I/II imbalance to correct
  make_option("--bmiq",     type = "logical",   default = FALSE),
  make_option("--keep_gset",type = "logical",   default = TRUE,
              help = "store the GenomicRatioSet in harmonised.rds; required by the legacy block baseline"),
  make_option("--seed",     type = "integer",   default = 20240717)
)))
set.seed(opt$seed)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

ss <- as.data.frame(fread(opt$sheet))
stopifnot(all(c("Sample_Name","Basename","Array_Type","Subject_ID") %in% names(ss)))
ss$Basename_full <- file.path(opt$idat_dir, ss$Basename)
rownames(ss) <- ss$Sample_Name

# ---- 1. read per array type -------------------------------------------------
read_one <- function(sub, arr) {
  log("reading ", nrow(sub), " ", arr, " samples")
  rg <- read.metharray.exp(targets = data.frame(Basename = sub$Basename_full),
                           force = TRUE, verbose = FALSE)
  colnames(rg) <- sub$Sample_Name
  rg@annotation <- c(
    array      = if (arr == "450K") "IlluminaHumanMethylation450k" else "IlluminaHumanMethylationEPIC",
    annotation = if (arr == "450K") "ilmn12.hg19" else "ilm10b4.hg19")
  rg
}
rg_list <- lapply(split(ss, ss$Array_Type), function(s) read_one(s, s$Array_Type[1]))

# ---- 2. per-array QC, computed on the NATIVE probe space --------------------
# Detection p-values must be computed before intersection: a probe dropped by
# harmonisation can still flag a bad sample.
qc_rows <- list(); detp_list <- list()
for (arr in names(rg_list)) {
  rg <- rg_list[[arr]]
  dp <- detectionP(rg)
  ms <- preprocessRaw(rg)
  q  <- getQC(ms)
  qc_rows[[arr]] <- data.table(
    Sample_Name = colnames(rg), Array_Type = arr,
    mMed = as.numeric(q$mMed), uMed = as.numeric(q$uMed),
    qc_mean = (as.numeric(q$mMed) + as.numeric(q$uMed)) / 2,
    frac_failed_probes = colMeans(dp > opt$detp),
    predicted_sex = getSex(mapToGenome(ms))$predictedSex)
  detp_list[[arr]] <- dp
  rm(ms, q); gc()
}
qc <- rbindlist(qc_rows)
qc <- merge(qc, as.data.table(ss)[, .(Sample_Name, Subject_ID, sex, cohort,
                                      Sentrix_ID, visit, days_on_clozapine)],
            by = "Sample_Name")
qc[, sex_mismatch := toupper(substr(sex,1,1)) != toupper(substr(predicted_sex,1,1))]
# median-intensity threshold is the minfi convention; failed-probe fraction
# catches samples that pass on intensity but have poor bisulfite conversion.
qc[, sample_fail := qc_mean < 10.5 | frac_failed_probes > opt$detp_frac]
fwrite(qc, file.path(opt$out_dir, "qc_samples.csv"))
log("samples failing QC: ", sum(qc$sample_fail), "; sex mismatches: ", sum(qc$sex_mismatch))

keep_samp <- qc[sample_fail == FALSE, Sample_Name]
rg_list <- lapply(rg_list, function(rg) rg[, colnames(rg) %in% keep_samp])

# ---- 3. noob per array, then combine ---------------------------------------
# combineArrays() on MethylSet keeps peak memory ~half of the RGChannelSet
# route and is equivalent, because noob is a within-array correction.
ms_list <- lapply(names(rg_list), function(arr) {
  log("noob on ", arr); preprocessNoob(rg_list[[arr]])
})
names(ms_list) <- names(rg_list)
rm(rg_list); gc()

log("combineArrays -> virtual 450K probe space")
ms <- Reduce(function(a, b) combineArrays(a, b,
             outType = "IlluminaHumanMethylation450k"), ms_list)
rm(ms_list); gc()
gr <- mapToGenome(ratioConvert(ms, what = "both"))
rm(ms); gc()
log("combined: ", nrow(gr), " probes x ", ncol(gr), " samples")

# ---- 4. probe filtering ----------------------------------------------------
ann <- getAnnotation(gr)
# Each array's detection p-values live in that array's own probe space
# (485512 rows for 450K, 865859 for EPIC), so they must be intersected on
# probe name -- and put in a common row order -- before they can be bound.
detp_probes <- Reduce(intersect, lapply(detp_list, rownames))
detp <- do.call(cbind, lapply(detp_list, function(d)
  d[detp_probes, colnames(d) %in% keep_samp, drop = FALSE]))
common <- intersect(rownames(gr), rownames(detp))
detp <- detp[common, colnames(gr), drop = FALSE]

filt <- data.table(probe = rownames(gr))
filt[, fail_detp := FALSE]
filt[match(common, probe), fail_detp := rowMeans(detp > opt$detp) > 0]
filt[, is_sex := ann$chr[match(probe, rownames(ann))] %in% c("chrX", "chrY")]
# cross-reactive / multi-mapping probes (Chen 2013 for 450K, Pidsley 2016 EPIC)
xr_file <- Sys.getenv("EWAS_CROSSREACTIVE", "")
filt[, is_crossreactive := FALSE]
if (nzchar(xr_file) && file.exists(xr_file))
  filt[probe %in% fread(xr_file, header = FALSE)$V1, is_crossreactive := TRUE]
gr <- dropLociWithSnps(gr, snps = c("SBE", "CpG"), maf = 0.01)
filt[, dropped_snp := !(probe %in% rownames(gr))]
drop <- filt[fail_detp | is_crossreactive | dropped_snp |
             (opt$drop_sex & is_sex), probe]
gr <- gr[!(rownames(gr) %in% drop), ]
fwrite(filt, file.path(opt$out_dir, "probe_filter.csv"))
log("after filtering: ", nrow(gr), " probes")

# ---- 5. BMIQ: type-I / type-II scale harmonisation -------------------------
beta0 <- getBeta(gr)
beta  <- beta0
if (opt$bmiq) {
  suppressPackageStartupMessages(library(wateRmelon))
  design <- ifelse(getAnnotation(gr)$Type == "I", 1L, 2L)
  log("BMIQ on ", ncol(beta0), " samples")
  beta <- vapply(seq_len(ncol(beta0)), function(j) {
    b <- beta0[, j]; ok <- is.finite(b) & b > 0 & b < 1
    out <- b
    fit <- try(wateRmelon::BMIQ(b[ok], design.v = design[ok], plots = FALSE,
                                nfit = 50000), silent = TRUE)
    if (!inherits(fit, "try-error")) out[ok] <- fit$nbeta
    out
  }, numeric(nrow(beta0)))
  dimnames(beta) <- dimnames(beta0)
}
# M-values are the modelling scale (variance-stabilised); beta is for reporting.
bc   <- pmax(pmin(beta, 1 - 1e-6), 1e-6)
mval <- log2(bc / (1 - bc))

pheno <- as.data.table(ss)[Sample_Name %in% colnames(beta)]
pheno <- pheno[match(colnames(beta), Sample_Name)]
ann_out <- as.data.table(
  as.data.frame(getAnnotation(gr))[, c("chr","pos","strand","Type",
                                       "Relation_to_Island","UCSC_RefGene_Name")],
  keep.rownames = "probe")

# The legacy baseline (03_baseline_bumphunter.R) needs the GenomicRatioSet
# itself, because cpgCollapse()/blockFinder() only accept a minfi object. Write
# the harmonised beta back into it so both stages see the identical matrix.
gset <- NULL
if (opt$keep_gset) {
  gset <- gr
  ok <- try({ SummarizedExperiment::assay(gset, "Beta") <- beta; TRUE },
            silent = TRUE)
  if (inherits(ok, "try-error")) {
    log("!! could not write harmonised beta into the GenomicRatioSet; ",
        "storing the pre-BMIQ object instead")
    gset <- gr
  }
}

saveRDS(list(beta = beta, mval = mval, pheno = pheno, anno = ann_out,
             gset = gset),
        file.path(opt$out_dir, "harmonised.rds"))

# Plain float64 column-major dump + JSON header, so the Python region-finding
# modules read the same matrix without an R<->Python bridge dependency.
export_f64 <- function(mat, stem) {
  con <- file(paste0(stem, ".f64"), "wb")
  writeBin(as.vector(mat), con, size = 8); close(con)
  writeLines(jsonlite::toJSON(list(nrow = nrow(mat), ncol = ncol(mat),
    order = "F", rownames = rownames(mat), colnames = colnames(mat)),
    auto_unbox = TRUE), paste0(stem, "_dims.json"))
}
suppressPackageStartupMessages(library(jsonlite))
export_f64(mval, file.path(opt$out_dir, "mval"))
fwrite(ann_out, file.path(opt$out_dir, "probe_annotation.csv"))
fwrite(pheno,   file.path(opt$out_dir, "pheno_used.csv"))
log("wrote harmonised.rds: ", nrow(beta), " probes x ", ncol(beta), " samples")
