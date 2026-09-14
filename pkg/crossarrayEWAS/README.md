# crossarrayEWAS

Design-aware region and block detection for longitudinal EWAS that combines
Illumina 450K and EPIC samples. This is the numerical core of the
[ewas-crossarray-harmonise](https://github.com/kkamieniecka/ewas-crossarray-harmonise)
pipeline, packaged so the Galaxy wrappers and the Nextflow stages can depend on
a version rather than vendoring the script.

## Installation

Once accepted into Bioconductor:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("crossarrayEWAS")
```

Development version:

```r
BiocManager::install("kkamieniecka/ewas-crossarray-harmonise",
                     subdir = "pkg/crossarrayEWAS")
```

## Scope

`fit_within()` estimates the exposure effect within subject for every probe,
absorbing subject, array type and chip. `comethylation_clusters()` groups
adjacent probes that co-vary within subject. `distance_penalty()` and
`tv_denoise()` produce a piecewise-constant effect track whose jumps are the
region boundaries, which `call_regions()` reports with a precision-weighted
effect and standard error. `fit_block_hsmm()` and `call_blocks()` call
long-range blocks from a distance-aware three-state hidden semi-Markov model.
`stability_selection()` and `within_subject_permutation()` resample subjects,
not samples.

Status: pre-submission skeleton. See `NEWS.md` for what is still missing.
