# Manuscript figures

Scripts that build the figures for the applications note from the recorded
GSE237561 run. They read result tables and manifests only; they never rerun a
pipeline stage and they never generate data.

| Script | Produces |
|---|---|
| `01_probe_sets.R` | `probes_450k_native.csv`, `probes_epic_native.csv`, `probes_shared_manifest.csv` and `sessionInfo_R.txt` — the three probe sets of the spacing baseline, filtered identically |
| `02_figure1.py` | `docs/figures/figure1.pdf`, `.tiff` and `.png`, plus `figure1_stats.json` |

## Running them

```sh
Rscript figures/01_probe_sets.R figures           # needs minfi and both annotation packages
python  figures/02_figure1.py  figures            # writes the figure beside its inputs
```

`01_probe_sets.R` writes three large CSVs of probe coordinates, which are
regenerable from the annotation packages in under a minute and are therefore
ignored by git rather than committed. The summary it feeds,
`results/GSE237561/spacing_baseline.csv`, is committed.

`02_figure1.py` expects its inputs in the directory given as the argument:
`spacing_baseline.csv`, `gaps.npz` and the four panel tables, which are the
committed copies under `results/GSE237561/` with the names recorded in
`figure1_provenance.json`. Copy them in, or point the script at a working
directory holding both.

## Probe sets and filters

All three sets carry the two filters that can be applied to a manifest without
samples: autosomes only, and loci dropped by
`dropLociWithSnps(snps = c("SBE", "CpG"), maf = 0.01)`, which is the call
`bin/01_harmonise.R` makes. Detection-p filtering is sample-dependent and is
not applied to a native set, so the native rows differ from the analysed panel
by that filter as well as by probe content. `docs/figure1-report.md` records the
counts and the one place where a published figure and the panel it describes
come from different probe sets.
