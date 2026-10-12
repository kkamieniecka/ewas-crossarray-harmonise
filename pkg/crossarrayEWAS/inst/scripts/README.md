# Scripts

`make_extdata.R` regenerates `inst/extdata/GSE237561_chr21_subset.rds`, the 450K/EPIC subset the vignette analyses, from the public GEO series GSE237561 and the harmonised probe map. It downloads ~250 MB and takes a few minutes, which is why the result ships with the package.

# Pipeline drivers

The command-line stage drivers (`04_dmr_ml.R`, `05_blocks_hsmm.R`) live in the pipeline repository, which generates this package. They do what `ewas_regions()` and `ewas_blocks()` do, with file handling, run records and provenance around it, and they source the pipeline's own copy of the numerical core rather than loading this package.
