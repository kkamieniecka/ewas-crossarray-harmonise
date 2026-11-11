# Regenerating the harmonised M-value matrix from the GSE237561 IDATs

Nothing has been started. This is what it would take, with the recorded
parameters, and where each figure comes from.

## What is missing and what survives

`harmonised.rds` is in no archive. What survives is enough to verify a rebuild
rather than to avoid one: `results/harmonised/mval_dims.json` carries the
422,898 row names and the 126 column names in order, so a regenerated matrix can
be checked row-for-row and column-for-column against the one the published
results were computed on. `probe_filter.csv` carries the per-probe filter flags,
and `qc_samples.csv` the per-sample QC values, so the filtering decisions can be
compared too.

The IDATs themselves are also gone. `data/sample_sheets/sample_sheet_GSE237561.csv`
holds the per-sample supplementary-file URLs, and `tools/fetch_idats.py` fetches
from those URLs and verifies each file as a complete gzip stream, which exists
because a truncated IDAT only surfaces later as an obscure reader error.

## The recorded parameters

From `bin/01_harmonise.R` and confirmed against `results/GSE237561/logs/harmonise.log`:

| Parameter | Value | Note |
|---|---|---|
| `--detp` | 0.01 | detection-p threshold |
| `--detp_frac` | 0.05 | sample fails if more than this fraction of probes fail detection |
| `--drop_sex` | TRUE | chrX and chrY dropped |
| `--bmiq` | FALSE | off by default: every shared probe has the same Infinium design type on both arrays |
| `--keep_gset` | TRUE | the GenomicRatioSet is stored, which the blockFinder baseline needs |
| `--seed` | 20240717 | |
| Normalisation | `preprocessNoob` per array type | within-array dye-bias and background correction, applied before `combineArrays` |
| Merge | `combineArrays(outType = "IlluminaHumanMethylation450k")` | on MethylSets, never by binding beta matrices |
| SNP filter | `dropLociWithSnps(snps = c("SBE", "CpG"), maf = 0.01)` | 16,098 probes |
| Cross-reactive filter | not applied | the list is read only when `EWAS_CROSSREACTIVE` names a file, and it was unset |
| Sample QC | `qc_mean < 10.5` or failed-probe fraction above `--detp_frac` | 0 samples failed, 0 sex mismatches |

**One thing to check against the thesis before rebuilding.** The run normalised
with noob, not with quantile normalisation. If the Methods text says "quantile
normalisation within array type", it does not describe this run, and a rebuild
under the recorded parameters will not produce what that sentence claims.

## What it needs

- **Packages**: `minfi`, both manifests and both annotation packages
  (`IlluminaHumanMethylation450kmanifest`, `...EPICmanifest`,
  `...450kanno.ilmn12.hg19`, `...EPICanno.ilm10b4.hg19`), plus `data.table` and
  `optparse`. All are in the restored `.r-libs/ewas` bundle, which also needs
  `matrixStats` in the environment for `MatrixGenerics` to load. `conf/` carries
  the conda specs and the source-install script.
- **Memory**: the recorded run's garbage-collector report shows a high-water
  mark of 10,350.6 MB of vector memory. On a 16 GB machine that leaves little
  headroom; it will run, but not alongside anything else.
- **Disk**: the IDATs plus the stored GenomicRatioSet. The 252 IDAT files are
  a per-sample download whose exact size I have not measured — a HEAD request
  against the URLs in the sample sheet gives the bytes, and I can run that in
  under a minute when you want it. `harmonised.rds` at 422,898 x 126 with Beta,
  M and CN retained is of order 1 GB in memory; compressed on disk it will be
  several hundred megabytes.

## Time

| Step | Measured or estimated | Source |
|---|---|---|
| IDAT download, 252 files | not measured | depends on the GEO link; `tools/fetch_idats.py` verifies each file, so a failed fetch costs a retry rather than a silent corruption |
| Stage 01, read to `harmonised.rds` | **5 min 34 s** | `harmonise.log`, 16:54:22 to 16:59:56 |
| of which: reading 126 arrays | 3 min 9 s | same log |
| of which: noob on both arrays | 36 s | same log |
| of which: combine, filter, write | 1 min 15 s | same log |

So the matrix itself is quick once the IDATs are local; the download dominates.
If the whole chain is wanted rather than the matrix alone, the recorded
downstream runtimes are 338 s for the bumphunter and blockFinder baseline and
9,399 s for the distance-aware stage at 200 permutations and 50 bootstrap
subsamples, which is where the hours are.

## Suggested order

1. `tools/fetch_idats.py` with the committed sample sheet, and let it verify.
2. Stage 01 with the parameters above, into a scratch directory.
3. Check the regenerated matrix against `mval_dims.json`: the row names and
   column names must match in order, and `probe_filter.csv` and
   `qc_samples.csv` must reproduce. If they do, every published downstream
   number is reproducible from this matrix; if they do not, stop and find out
   why before rerunning anything downstream.
