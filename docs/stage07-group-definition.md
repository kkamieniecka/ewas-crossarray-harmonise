# Stage 07 re-run with an explicit group definition

## What changed in the pipeline

`bin/07_pgx_genotypes.R` now states, computes and records which partition of the
samples a gapped probe is being interpreted on.

1. New option `--group-definition {gap,class}`, default `gap`. The `gap`
   partition is the one gaphunter itself returns, taken from `gh$sampleresults`
   when present and otherwise reconstructed by cutting the ordered beta values
   wherever consecutive values differ by more than `--gap-threshold`. The
   `class` partition is the `--het-lo`/`--het-hi` boundary pair used to call
   identity-probe genotypes; it is retained only for comparison with earlier
   runs, and the header comment says why it is not a group definition: a fixed
   boundary applied to a CpG probe can cut through a unimodal distribution and
   manufacture two classes where the data contain no gap.
2. Both partitions are written out. `pgx_cyp_groups.csv` is new and carries gap
   group membership per sample; `pgx_cyp_genotypes.csv` keeps the class labels
   and is now commented as such.
3. `pgx_cyp_probe_interpretation.csv` is computed on the selected partition and
   gains the confounding and frequency columns: group size in samples and in
   subjects, the per-subject visit counts inside and outside the group, the
   number of subjects straddling it, the number of Sentrix chips the group
   occupies, the top-chip share, a permutation P value for group membership
   against chip (`--n-perm`, default 5,000, seeded by `--seed`), a Fisher exact
   P value against array platform, one-way variance explained in beta by chip,
   subject and platform, and, where the manifest documents a variant, the
   expected number of carrier subjects under Hardy-Weinberg computed on
   independent subjects with an exact binomial P value. The class split is kept
   as a labelled column so the two partitions can be compared directly.
4. The run summary records `group_definition`, `gap_threshold`, `het_lo` and
   `het_hi`, so a reader of the outputs never has to infer which partition was
   used.

## What could and could not be re-executed

The harmonised object and the IDATs are not in this workspace, so the detection
step itself, gap hunting across the 668 CYP-annotated probes, was not
re-executed and the set of five flagged probes is carried over from the
original run. Everything downstream of detection was recomputed from the
deposited per-probe beta matrix (`pgx_cyp_beta.csv`) and the sample sheet.

The patched script parses cleanly, and its three new helpers were extracted and
run against the deposited betas. They reproduce the independent Python
computation exactly for group size, subject composition, straddling counts,
chip counts and all three variance-explained terms. The permutation P value is
the one quantity that moves between runs: at 2,000 permutations cg14711976
returned 0.057 and at 5,000 it returned 0.048, so the test was settled at
50,000.

## Results under the explicit definition

| Probe | Gene | Group | Subjects | Visits in group | Chips | Chip P (50,000 perm) | Platform P |
|---|---|---|---|---|---|---|---|
| cg08561071 | CYP2F1 | 5 samples | 2 | 0044: 2 of 2; 0048: 3 of 4 | 2 of 14 | 0.029 | 0.32 |
| cg07512206 | CYP4F8 | 4 samples | 1 | 0010: 4 of 4 | 1 of 14 | 0.0033 | 0.57 |
| cg14711976 | CYP4F12 | 4 samples | 1 | 0042: 4 of 4 | 2 of 14 | 0.057 | 0.57 |
| cg09199330 | CYP27C1 | 3 samples | 1 | CLZ1022: 3 of 3 | 1 of 14 | 0.00002 | 0.018 |
| cg10985347 | CYP27C1 | 3 samples | 1 | CLZ1022: 3 of 3 | 1 of 14 | 0.00004 | 0.018 |

Two figures move against what was reported before this re-run: cg07512206 from
0.0024 to 0.0033, and cg14711976 from 0.048 to 0.057, so the chip association at
that probe does not reach significance on its own. The conclusion rests on the
set rather than on any single probe and is unaffected, but the manuscript now
says so explicitly rather than listing five significant associations.

## Figures

Both figures that encoded the class split were redrawn on the gap partition.

- `fig1_pgx_identity_and_gaps.png`. Panel (a) is unchanged in substance: all
  7,875 identity-probe concordances by pair type, with the 0.339 separation
  shaded. Panel (b) plots every sample's beta at the five probes with the gap
  boundary drawn as a dashed rule and the gap group picked out, replacing the
  version in which the drawn boundary was the 0.25/0.75 call threshold.
- `fig3_pgx_within_subject_invariance_v2.png`. Per-subject beta ranges across
  visits, five panels, with the gap boundary dashed and each contributing
  subject labelled with how many of its visits fall inside the group. The
  earlier version drew the class boundary, which is why a third of the subjects
  appeared to straddle it at cg07512206.

Both are embedded in the revised manuscript in place of the originals, width
preserved and height set from the new aspect ratio.

## Still outstanding

Re-running the detection step needs the harmonised object, which means
re-downloading the 1.5 GB IDAT archive from GEO and re-running stage 01. That is
worth doing once before submission so the flagged set, and not only its
interpretation, comes from the current code. A genome-wide gap-hunting run would
also supply the background rate against which five flagged probes out of 668 can
be judged; no run has produced one yet.
