# Pipeline drivers

The command-line stage drivers (`04_dmr_ml.R`, `05_blocks_hsmm.R`) live in the pipeline repository and are not yet installed here. They call only exported functions, and prefer the installed package over sourcing the pipeline's `ewasml.R` when one is present.
