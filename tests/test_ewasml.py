"""
Numerical checks for ewasml. Each test states the property being verified and
the tolerance, so a failure localises to one estimator.

Run:  python ewas-harmonise/tests/test_ewasml.py
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
import ewasml as E

rng = np.random.default_rng(0)
FAIL = []


def check(name, cond, detail=""):
    print(f"{'PASS' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not cond:
        FAIL.append(name)


# --- 1. within-subject estimator recovers a known effect under a large,
#        perfectly confounded array/cohort offset -----------------------------
n_sub, n_vis, n_probe = 30, 4, 400
subject = np.repeat(np.arange(n_sub), n_vis)
time = np.tile(np.arange(n_vis), n_sub).astype(float)
array = (subject >= n_sub // 2).astype(float)          # array nested in subject
true_beta = np.zeros(n_probe)
true_beta[:40] = 0.30                                   # 40 true probes
subj_off = rng.normal(0, 2.0, n_sub)[subject]
arr_off = 1.5 * array                                   # huge array shift
M = (true_beta[:, None] * time[None, :]
     + subj_off[None, :] + arr_off[None, :]
     + rng.normal(0, 0.25, (n_probe, len(subject))))
fit = E.fit_within(M, time, subject, shrink_var=False)
bias = np.abs(fit.beta[:40].mean() - 0.30)
null_bias = np.abs(fit.beta[40:].mean())
check("fit_within recovers effect despite confounded array offset",
      bias < 0.02 and null_bias < 0.02,
      f"mean beta(true)={fit.beta[:40].mean():.4f} (truth 0.30), "
      f"mean beta(null)={fit.beta[40:].mean():.4f}")
check("fit_within z separates true from null probes",
      np.abs(fit.z[:40]).mean() > 5 * np.abs(fit.z[40:]).mean(),
      f"|z| true={np.abs(fit.z[:40]).mean():.2f} null={np.abs(fit.z[40:]).mean():.2f}")

# time-invariant covariate must be annihilated, not break the fit
age = rng.normal(40, 10, n_sub)[subject]
fit2 = E.fit_within(M, time, subject, covars=age[:, None], shrink_var=False)
check("time-invariant covariate is dropped by the within transform",
      np.allclose(fit2.beta, fit.beta, atol=1e-10))

# --- 2. variance moderation shrinks towards the trend ----------------------
fit3 = E.fit_within(M, time, subject, shrink_var=True)
spread_raw = np.std(fit.se); spread_mod = np.std(fit3.se)
check("variance moderation reduces SE dispersion",
      spread_mod < spread_raw, f"sd(se) {spread_raw:.4f} -> {spread_mod:.4f}")

# --- 3. TV denoising recovers a piecewise-constant track -------------------
n = 600
pos = np.sort(rng.integers(0, 400_000, n)).astype(float)
truth = np.zeros(n)
truth[150:230] = 0.4
truth[400:460] = -0.35
se = rng.uniform(0.05, 0.20, n)
y = truth + rng.normal(0, se)
wt = 1.0 / se ** 2
lam = E.distance_penalty(pos, lam0=1.0, decay_bp=2000.0, w=wt)
theta = E.tv_denoise(y, wt, lam)
mse_raw = np.mean((y - truth) ** 2); mse_tv = np.mean((theta - truth) ** 2)
check("tv_denoise reduces MSE against the true step signal",
      mse_tv < 0.5 * mse_raw, f"MSE {mse_raw:.4f} -> {mse_tv:.4f}")
check("tv_denoise output is piecewise constant (few distinct levels)",
      len(np.unique(np.round(theta, 6))) < n / 3,
      f"{len(np.unique(np.round(theta,6)))} levels for n={n}")
theta_flat = E.tv_denoise(y, 1.0 / se ** 2, np.full(n - 1, 1e6))
w = 1.0 / se ** 2
check("tv_denoise -> weighted mean as penalty grows",
      np.allclose(theta_flat, (w * y).sum() / w.sum(), atol=1e-3))

# --- 4. region calling finds the two simulated regions ---------------------
chrom = np.array(["chr1"] * n)
regs = E.call_regions(chrom, pos, y, se, theta, min_probes=5, min_effect=0.10)
hits = [r for r in regs.table if abs(r["effect_M"]) > 0.15]
check("call_regions recovers both simulated regions", len(hits) >= 2,
      f"{len(hits)} regions with |effect|>0.15")

# --- 5. within-subject permutation preserves subject structure -------------
perm = E.within_subject_permutation(time, subject, rng)
same_multiset = all(
    np.array_equal(np.sort(perm[subject == s]), np.sort(time[subject == s]))
    for s in np.unique(subject))
check("within_subject_permutation keeps each subject's exposure multiset",
      same_multiset)
check("within_subject_permutation actually permutes", not np.array_equal(perm, time))
zmax = []
for b in range(60):
    p = E.within_subject_permutation(time, subject, rng)
    f = E.fit_within(M, p, subject, shrink_var=False)
    zmax.append(np.abs(f.z).max())
obs = np.abs(fit.z).max()
check("observed max|z| exceeds the within-subject permutation null",
      obs > np.quantile(zmax, 0.95),
      f"obs={obs:.1f}, null 95th pct={np.quantile(zmax,0.95):.1f}")

# --- 6. block HSMM recovers simulated blocks -------------------------------
nk = 500
kpos = np.sort(rng.integers(0, 60_000_000, nk)).astype(float)
kchr = np.array(["chr1"] * nk)
bstate = np.zeros(nk)
bstate[120:200] = 1.0     # hyper block
bstate[330:400] = -1.0    # hypo block
kse = rng.uniform(0.02, 0.08, nk)
kyy = 0.25 * bstate + rng.normal(0, kse)
hs = E.fit_block_hsmm(kyy, kse, kpos, kchr, length_scale=2_000_000.0, seed=1)
decoded = hs.posterior.argmax(axis=1) - 1
acc_in = np.mean(decoded[120:200] == 1), np.mean(decoded[330:400] == -1)
fp = np.mean(decoded[np.abs(bstate) < 0.5] != 0)
check("block HSMM recovers hyper and hypo blocks",
      acc_in[0] > 0.7 and acc_in[1] > 0.7,
      f"hyper recall={acc_in[0]:.2f}, hypo recall={acc_in[1]:.2f}")
check("block HSMM keeps neutral background mostly neutral", fp < 0.15,
      f"false-positive cluster rate={fp:.3f}")
blocks = E.call_blocks(kchr, kpos, kpos, hs.posterior, min_post=0.8, min_clusters=3)
check("call_blocks returns both blocks", len(blocks) >= 2, f"{len(blocks)} blocks")
check("HSMM state means are ordered hypo<neutral<hyper",
      hs.mu[0] < hs.mu[1] <= hs.mu[2], f"mu={np.round(hs.mu,3)}")

# --- 7. co-methylation clustering ------------------------------------------
npc = 300
cpos = np.arange(npc) * 120.0
cchr = np.array(["chr1"] * npc)
shared = rng.normal(0, 1, (3, 40))
resid = rng.normal(0, 1, (npc, 40))
resid[50:56] = shared[0] + rng.normal(0, 0.2, (6, 40))    # co-methylated run
resid[150:158] = shared[1] + rng.normal(0, 0.2, (8, 40))
cl = E.comethylation_clusters(cchr, cpos, resid, max_gap=1000, rho_min=0.5)
c1 = set(cl[50:56]); c2 = set(cl[150:158])
check("co-methylated runs form single clusters",
      len(c1) == 1 and -1 not in c1 and len(c2) == 1 and -1 not in c2,
      f"cluster ids {c1}, {c2}")
check("uncorrelated neighbours are not clustered",
      np.mean(cl[:40] == -1) > 0.8, f"{np.mean(cl[:40]==-1):.2f} singletons")

print("\n" + ("ALL PASS" if not FAIL else f"{len(FAIL)} FAILED: {FAIL}"))
sys.exit(1 if FAIL else 0)
