#!/usr/bin/env Rscript
## 01_probe_sets.R -- native 450K, native EPIC and shared probe sets under one
## filter rule, for the probe-spacing baseline of Figure 1B.
##
## The filters are the two the task specifies and the two that can be applied to
## a manifest without samples: autosomes only, and loci dropped by
## dropLociWithSnps(snps = c("SBE", "CpG"), maf = 0.01), which is the call
## bin/01_harmonise.R makes on the shared panel. Detection-p filtering is
## sample-dependent and cannot be reproduced for a native set without IDATs;
## cross-reactive filtering was NOT active in the recorded run (probe_filter.csv
## has is_crossreactive = 0 for all 452,567 probes, because EWAS_CROSSREACTIVE
## was unset), so it is not applied here either.
##
## Usage: Rscript fig1/01_probe_sets.R <out-dir>

suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})

out <- commandArgs(trailingOnly = TRUE)[1]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

AUTOSOMES <- paste0("chr", 1:22)

probe_set <- function(anno_pkg, label) {
  ann <- minfi::getAnnotation(anno_pkg)
  n0 <- nrow(ann)
  keep_auto <- ann$chr %in% AUTOSOMES
  snp_cpg <- !is.na(ann$CpG_maf) & ann$CpG_maf >= 0.01
  snp_sbe <- !is.na(ann$SBE_maf) & ann$SBE_maf >= 0.01
  keep <- keep_auto & !(snp_cpg | snp_sbe)
  cat(sprintf("%-14s manifest %7d | autosomal %7d | SNP at CpG %6d | SNP at SBE %6d | kept %7d\n",
              label, n0, sum(keep_auto), sum(snp_cpg & keep_auto),
              sum(snp_sbe & keep_auto), sum(keep)))
  data.frame(probe = rownames(ann)[keep],
             chr = sub("^chr", "", ann$chr[keep]),
             pos = ann$pos[keep],
             stringsAsFactors = FALSE)
}

a450 <- probe_set(IlluminaHumanMethylation450kanno.ilmn12.hg19, "450K native")
aepic <- probe_set(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, "EPIC native")
shared <- a450[a450$probe %in% aepic$probe, ]
cat(sprintf("%-14s %7d probes present on both manifests after the same filters\n",
            "shared", nrow(shared)))

write.csv(a450, file.path(out, "probes_450k_native.csv"), row.names = FALSE)
write.csv(aepic, file.path(out, "probes_epic_native.csv"), row.names = FALSE)
write.csv(shared, file.path(out, "probes_shared_manifest.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_R.txt"))
