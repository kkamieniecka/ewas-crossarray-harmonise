"""Validation of the genotype-content audit's gap arm by simulation.

The identity probes validate genotype handling, but not gap detection, which is
a different procedure on a different probe class. This script simulates probes
on the cohort's own design (38 subjects, their observed visit counts, their
observed Sentrix chip assignment) and measures, per scenario:

  flagged          a gap wider than --threshold separates a group of >= 2 samples
  invariant        no subject has visits on both sides of the gap
  freq_consistent  the number of subjects in the group agrees with Hardy-Weinberg
                   at the known allele frequency (exact binomial, two-sided);
                   the expectation matches the class the group represents - the
                   rare homozygote when three groups are occupied, carriers when
                   two
  chip_separable   group membership is not associated with chip (permutation)
  all_four         all of the above, i.e. the probe would be reported as
                   genotype-compatible

The counterfactual layout, in which subjects are spread across chips and nothing
else changes, isolates design from method.

    python validate_audit.py sample_sheet.csv out.csv
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

THRESHOLD = 0.15
CLASS_MEANS = np.array([0.15, 0.50, 0.85])   # AA, AB, BB at a SNP-driven probe
NOISE_SD = 0.03                              # matched to observed within-subject SD
N_REP = 150
N_PERM = 1000


def chi2_stat(minor, codes, n_chip):
    a = np.zeros(n_chip); b = np.zeros(n_chip)
    np.add.at(a, codes[minor], 1); np.add.at(b, codes[~minor], 1)
    tot = a + b
    exp = np.outer(np.array([a.sum(), b.sum()]), tot) / len(codes)
    obs = np.vstack([a, b])
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(exp > 0, (obs - exp) ** 2 / exp, 0.0).sum()


def chip_p(minor, codes, n_chip, n_perm, rng):
    obs = chi2_stat(minor, codes, n_chip)
    idx = np.arange(len(minor)); hits = 0
    for _ in range(n_perm):
        rng.shuffle(idx)
        hits += chi2_stat(minor[idx], codes, n_chip) >= obs - 1e-12
    return (hits + 1) / (n_perm + 1)


def audit(beta, maf, subj, codes, n_chip, n_subj, rng, n_perm=N_PERM):
    o = np.argsort(beta); s = beta[o]
    if not (np.diff(s) > THRESHOLD).any():
        return dict(flagged=False)
    lab = np.zeros(len(beta), int)
    lab[o] = np.cumsum(np.r_[0, np.diff(s) > THRESHOLD])
    vc = pd.Series(lab).value_counts()
    minor = lab == vc.idxmin()
    if minor.sum() < 2:
        return dict(flagged=False)
    strad = any(0 < minor[subj == x].sum() < (subj == x).sum()
                for x in np.unique(subj[minor]))
    k = len(np.unique(subj[minor]))
    p_exp = maf ** 2 if len(vc) >= 3 else 2 * maf - maf ** 2
    p_freq = min(1.0, 2 * min(stats.binom.cdf(k, n_subj, p_exp),
                              stats.binom.sf(k - 1, n_subj, p_exp)))
    p_chip = chip_p(minor, codes, n_chip, n_perm, rng)
    return dict(flagged=True, invariant=not strad, freq_ok=p_freq >= 0.05,
                chip_ok=p_chip >= 0.05,
                called=bool((not strad) and p_freq >= 0.05 and p_chip >= 0.05))


def summarise(res, layout, scenario):
    f = [r for r in res if r.get("flagged")]
    pct = lambda k: round(100 * np.mean([r[k] for r in f]), 1) if f else np.nan
    return dict(layout=layout, scenario=scenario, n=len(res),
                flagged=round(100 * len(f) / len(res), 1),
                invariant=pct("invariant"), freq_consistent=pct("freq_ok"),
                chip_separable=pct("chip_ok"),
                all_four=round(100 * np.mean([r.get("called", False) for r in res]), 1))


def main(sheet_csv, out_csv, seed=23):
    ss = pd.read_csv(sheet_csv)
    subj = ss.Subject_ID.astype(str).values
    codes, chips = pd.factorize(ss.Sentrix_ID.values)
    n_chip, subs = len(chips), pd.unique(subj)
    n_subj, n = len(subs), len(subj)
    rng = np.random.default_rng(seed)
    spread = rng.permutation(np.repeat(np.arange(n_chip), int(np.ceil(n / n_chip)))[:n])

    def variant(maf):
        g = {s: v for s, v in zip(subs, rng.choice(
            [0, 1, 2], n_subj, p=[(1 - maf) ** 2, 2 * maf * (1 - maf), maf ** 2]))}
        return np.clip(CLASS_MEANS[np.array([g[x] for x in subj])]
                       + rng.normal(0, NOISE_SD, n), 0.01, 0.99)

    rows = []
    for layout, cc in [("observed", codes), ("subjects spread across chips", spread)]:
        for maf in (0.10, 0.20, 0.30):
            rows.append(summarise([audit(variant(maf), maf, subj, cc, n_chip, n_subj, rng)
                                   for _ in range(N_REP)],
                                  layout, "true variant, MAF %.2f" % maf))
    nulls = [("no signal", lambda: rng.normal(0.5, 0.05, n)),
             ("chip artefact, one chip shifted 0.3",
              lambda: rng.normal(0.5, 0.05, n) + 0.3 * (codes == rng.integers(n_chip))),
             ("subject artefact, one subject shifted 0.3",
              lambda: rng.normal(0.5, 0.05, n) + 0.3 * (subj == rng.choice(subs)))]
    for name, gen in nulls:
        rows.append(summarise([audit(np.clip(gen(), 0.01, 0.99), 0.05, subj, codes,
                                     n_chip, n_subj, rng) for _ in range(N_REP)],
                              "observed", name))
    out = pd.DataFrame(rows)
    out.to_csv(out_csv, index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
