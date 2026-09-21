#!/usr/bin/env Rscript
## Stage 07 -- what genotype information a 450K/EPIC methylation experiment
## actually carries, applied to the CYP pharmacogenes.
##
## Scope, stated up front because the name of the task invites a stronger
## reading than the assay supports: Infinium methylation arrays are not
## genotyping arrays. No CYP star allele (CYP2D6 *2/*4/*10, CYP2C19 *2/*17,
## CYP1A2 *1F, ...) is interrogated by any probe on either array, and CYP2D6
## structural variation -- deletions, duplications, hybrid genes, which drive
## much of its phenotype -- is not assessed here at all. Star-allele
## assignment needs targeted genotyping, an imputed GWAS panel or sequencing,
## and then a caller such as PharmCAT, Aldy or StellarPGx.
##
## What the arrays do carry, and what this stage calls:
##   1. The 'rs' identity probes (65 on 450K, 59 on EPIC), designed to be
##      genotyped rather than methylation-read. They are genome-wide sample
##      fingerprints, not pharmacogenes; used here to verify that the
##      longitudinal samples of a subject are the same person, including
##      across the 450K/EPIC boundary -- the assumption every within-subject
##      effect in this pipeline rests on.
##   2. cg/ch probes whose measured beta is driven by a SNP under the probe
##      body, the single-base extension site or the interrogated CpG. These
##      give a gapped, near-trimodal beta distribution that tracks genotype
##      rather than methylation. minfi::gaphunter is the standard detector;
##      the 450K annotation's Probe_rs / SBE_rs / CpG_rs columns say which
##      hits have a documented SNP to explain them. Restricted to probes
##      annotated to CYP genes, this is the set of CYP-locus variants the
##      experiment can genotype, and nothing beyond it.
##
## Inputs: the stage 01 harmonised object (betas for 2) and the IDATs (for the
## rs probes, which stage 01 filters out). Outputs to --out-dir:
##   pgx_rs_genotypes.csv          rs probe x sample calls (AA/AB/BB)
##   pgx_identity_concordance.csv  all sample pairs, within/between subject
##   pgx_cyp_variant_probes.csv    CYP probes with genotype-like beta gaps
##   pgx_cyp_genotypes.csv         per-sample classes at those probes
##   pgx_cyp_beta.csv              their betas, for plotting
##   pgx_summary.json              counts and parameters

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(minfi)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

opt <- parse_args(OptionParser(option_list = list(
  make_option("--in-dir",   type = "character", default = "results/harmonised"),
  make_option("--sheet",    type = "character"),
  make_option("--idat-dir", type = "character"),
  make_option("--out-dir",  type = "character", default = "results/pgx"),
  make_option("--genes",    type = "character", default = "(^|;)CYP",
              help = "regex on UCSC_RefGene_Name selecting the gene family"),
  make_option("--subject",  type = "character", default = "Subject_ID"),
  make_option("--array-col",type = "character", default = "Array_Type"),
  make_option("--het-lo",   type = "double", default = 0.25,
              help = "beta below this is the low homozygous class"),
  make_option("--het-hi",   type = "double", default = 0.75),
  make_option("--gap-threshold", type = "double", default = 0.15,
              help = "minimum beta gap for gaphunter"),
  make_option("--gap-outcutoff", type = "double", default = 0.01,
              help = "gaphunter outlier cutoff, as a fraction of samples")
)), convert_hyphens_to_underscores = TRUE)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------- 1. rs probes
## getSnpBeta() returns the rs probes' allelic ratio. Stage 01 discards them
## with the rest of the non-CpG panel, so the IDATs have to be re-read.
ss <- read.csv(opt$sheet, stringsAsFactors = FALSE)
ss$Basename <- file.path(opt$idat_dir, basename(ss$Basename))
arrays <- split(ss, ss[[opt$array_col]])
log_msg("reading rs probes from ", nrow(ss), " samples on ", length(arrays), " array types")

snp_list <- lapply(names(arrays), function(a) {
  s <- arrays[[a]]
  rg <- read.metharray.exp(targets = s, force = TRUE)
  colnames(rg) <- s$Sample_Name
  getSnpBeta(rg)
})
names(snp_list) <- names(arrays)
common_rs <- Reduce(intersect, lapply(snp_list, rownames))
log_msg("rs probes: ", paste(sprintf("%s=%d", names(arrays), sapply(snp_list, nrow)), collapse = ", "),
        "; shared by both arrays: ", length(common_rs))
snp_beta <- do.call(cbind, lapply(snp_list, function(m) m[common_rs, , drop = FALSE]))
snp_beta <- snp_beta[, ss$Sample_Name, drop = FALSE]

call_geno <- function(b, lo, hi)
  ifelse(is.na(b), NA_character_, ifelse(b < lo, "AA", ifelse(b > hi, "BB", "AB")))
rs_geno <- call_geno(snp_beta, opt$het_lo, opt$het_hi)

## Identity: fraction of rs probes carrying the same call in a pair of samples.
agree <- function(i, j) mean(rs_geno[, i] == rs_geno[, j], na.rm = TRUE)
subj <- ss[[opt$subject]]; arr <- ss[[opt$array_col]]
pairs <- t(combn(ncol(rs_geno), 2))
conc <- data.frame(
  sample_a = colnames(rs_geno)[pairs[, 1]], sample_b = colnames(rs_geno)[pairs[, 2]],
  subject_a = subj[pairs[, 1]], subject_b = subj[pairs[, 2]],
  array_a = arr[pairs[, 1]], array_b = arr[pairs[, 2]],
  same_subject = subj[pairs[, 1]] == subj[pairs[, 2]],
  cross_array = arr[pairs[, 1]] != arr[pairs[, 2]],
  concordance = mapply(agree, pairs[, 1], pairs[, 2]),
  stringsAsFactors = FALSE)
log_msg("identity concordance: within-subject median ",
        sprintf("%.3f", median(conc$concordance[conc$same_subject])),
        ", between-subject median ",
        sprintf("%.3f", median(conc$concordance[!conc$same_subject])))

## --------------------------------------------- 2. genotype-driven CYP probes
## Stage 01 saves a list (beta, mval, pheno, anno, gset); its annotation is the
## harmonised one and carries no SNP columns, so those are read from the 450K
## annotation package -- the space the two arrays were combined into.
h <- readRDS(file.path(opt$in_dir, "harmonised.rds"))
ann <- as.data.frame(h$anno)
stopifnot(identical(rownames(h$beta), ann$probe))
sel <- grepl(opt$genes, ann$UCSC_RefGene_Name)
log_msg("probes annotated to ", opt$genes, ": ", sum(sel))
bet <- h$beta[sel, , drop = FALSE][, ss$Sample_Name, drop = FALSE]
snp_ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

gh <- gaphunter(bet, threshold = opt$gap_threshold, keepOutliers = FALSE,
                outCutoff = opt$gap_outcutoff, verbose = FALSE)
hit <- rownames(gh$proberesults)
log_msg("gapped (genotype-like) probes among them: ", length(hit))

a  <- ann[match(hit, ann$probe), , drop = FALSE]
sa <- snp_ann[hit, , drop = FALSE]
snp_cols <- intersect(c("Probe_rs","Probe_maf","CpG_rs","CpG_maf","SBE_rs","SBE_maf"), names(sa))
documented <- !is.na(sa$Probe_rs) | !is.na(sa$CpG_rs) | !is.na(sa$SBE_rs)
vp <- data.frame(probe = hit, chr = a$chr, pos = a$pos,
                 gene = a$UCSC_RefGene_Name,
                 island = a$Relation_to_Island,
                 n_groups = gh$proberesults[hit, "Groups"],
                 sa[, snp_cols, drop = FALSE],
                 documented_snp = documented,
                 stringsAsFactors = FALSE)
vp <- vp[order(!vp$documented_snp, vp$chr, vp$pos), ]

## Per-sample calls at the gapped probes. The labels are the beta classes the
## gap defines, not reference/alternate alleles: the array reads an allelic
## ratio whose direction depends on probe design, so these are genotype
## CLASSES, and which class is the minor homozygote follows from its frequency.
cls <- t(apply(bet[hit, , drop = FALSE], 1, function(b)
  ifelse(b < opt$het_lo, "class_low", ifelse(b > opt$het_hi, "class_high", "class_mid"))))
colnames(cls) <- colnames(bet)

write.csv(data.frame(probe = rownames(rs_geno), rs_geno, check.names = FALSE),
          file.path(opt$out_dir, "pgx_rs_genotypes.csv"), row.names = FALSE)
write.csv(conc, file.path(opt$out_dir, "pgx_identity_concordance.csv"), row.names = FALSE)
write.csv(vp, file.path(opt$out_dir, "pgx_cyp_variant_probes.csv"), row.names = FALSE)
write.csv(data.frame(probe = rownames(cls), cls, check.names = FALSE),
          file.path(opt$out_dir, "pgx_cyp_genotypes.csv"), row.names = FALSE)
bh <- as.data.frame(bet[hit, , drop = FALSE])
write.csv(data.frame(probe = rownames(bh), bh, check.names = FALSE),
          file.path(opt$out_dir, "pgx_cyp_beta.csv"), row.names = FALSE)

summary_json <- list(
  n_samples = ncol(rs_geno),
  rs_probes_shared = length(common_rs),
  rs_probes_per_array = as.list(sapply(snp_list, nrow)),
  identity = list(
    within_subject_median = median(conc$concordance[conc$same_subject]),
    within_subject_min = min(conc$concordance[conc$same_subject]),
    between_subject_median = median(conc$concordance[!conc$same_subject]),
    between_subject_max = max(conc$concordance[!conc$same_subject]),
    cross_array_within_subject_median =
      if (any(conc$same_subject & conc$cross_array))
        median(conc$concordance[conc$same_subject & conc$cross_array]) else NA,
    n_pairs_same_subject = sum(conc$same_subject),
    n_pairs_cross_array_same_subject = sum(conc$same_subject & conc$cross_array)),
  gene_pattern = opt$genes,
  probes_in_genes = sum(sel),
  gapped_probes = length(hit),
  gapped_with_documented_snp = sum(documented),
  genes_with_gapped_probes = sort(unique(unlist(strsplit(vp$gene, ";")))),
  star_alleles_callable = FALSE,
  args = opt[setdiff(names(opt), "help")])
write_json(summary_json, file.path(opt$out_dir, "pgx_summary.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 6)
log_msg("done")
