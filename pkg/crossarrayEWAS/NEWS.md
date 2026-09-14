# crossarrayEWAS 0.99.0

* First version, extracted from the `ewas-crossarray-harmonise` pipeline
  (`bin/ewasml.R`). Function bodies are unchanged from the pipeline source;
  the package adds documentation, a namespace and unit tests.
* Not yet present, and required before submission: entry points accepting
  `SummarizedExperiment` / `GenomicRatioSet` and returning `GRanges`, a
  `BiocStyle` vignette, and the cross-array harmonisation stage (still
  straight-line script code in `scripts/01_harmonise.R`).
