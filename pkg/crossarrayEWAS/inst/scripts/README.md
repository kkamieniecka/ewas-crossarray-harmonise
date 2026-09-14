# Pipeline drivers

The command-line stage drivers (`04_dmr_ml.R`, `05_blocks_hsmm.R`) live in the pipeline repository and are not yet installed here. They call the exported functions plus one internal helper (`pinv()`), which has to be replaced by `MASS::ginv()` or dropped before they can run against the installed package.
