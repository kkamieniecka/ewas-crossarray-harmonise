#!/usr/bin/env python3
"""
Stage-level check: bin/05_blocks_hsmm.R must reproduce bin/05_blocks_hsmm.py.

tests/test_equivalence.R proves the numerical core agrees function by function.
This runs the two stage DRIVERS end to end on the same generated input and
compares every output file. It exists because both defects the port exposed
lived in driver code, out of reach of the function-level fixtures:

  * direction was labelled by assuming the middle HMM state is the zero-effect
    one, so when all cluster effects fall on one side of zero a real effect
    state was reported as "no change"
  * with an array arm skipped or no legacy units surviving, non-finite values
    reached the parameters record as a bare NaN (Python) or the string "NA"
    (R), so neither file parsed strictly and the two disagreed

Baum-Welch from fixed starting values uses no RNG, so this stage is held to
exact agreement -- including the called blocks -- not agreement in
distribution as stage 04 requires.

Usage
  python tests/gen_equivalence_fixtures.py --out-dir fx --stage-dir stage
  python tests/test_stage05_equivalence.py stage
"""
from __future__ import annotations

import json, os, subprocess, sys, tempfile
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "bin")
FAILED: list[str] = []


def check(name, ok, detail=""):
    print(f"{'ok  ' if ok else 'FAIL'}  {name}{'' if not detail else '  -- ' + detail}")
    if not ok:
        FAILED.append(name)


def run(cmd, out_dir, stage_dir, covars):
    subprocess.run(cmd + ["--in-dir", stage_dir, "--out-dir", out_dir,
                          "--probe-map", os.path.join(stage_dir,
                                                      "crossarray_probe_map.csv.gz"),
                          "--exposure-scale", "1", "--covars", covars,
                          # the fixture's dense probe runs collapse into too few
                          # clusters at the default 1500 bp gap; a tighter gap
                          # splits them, which is a fixture choice and not a
                          # change to the stage defaults
                          "--max-gap", "100", "--rho-min", "0.05",
                          "--min-clusters", "3", "--fixed-collapse"],
                   check=True, capture_output=True, text=True)


def strict_json(path):
    """Reject NaN/Infinity tokens: they are not JSON and jsonlite cannot read
    them, so a driver that emits one has produced an unusable record."""
    def boom(tok):
        raise ValueError(f"non-standard JSON token {tok!r} in {path}")
    with open(path) as fh:
        return json.load(fh, parse_constant=boom)


def compare(py_dir, r_dir, tag):
    bl_p = pd.read_csv(os.path.join(py_dir, "blocks_hsmm.csv"))
    bl_r = pd.read_csv(os.path.join(r_dir, "blocks_hsmm.csv"))
    check(f"[{tag}] block table has the same shape and columns",
          bl_p.shape == bl_r.shape and list(bl_p.columns) == list(bl_r.columns),
          f"{bl_p.shape} vs {bl_r.shape}")
    if bl_p.shape == bl_r.shape and list(bl_p.columns) == list(bl_r.columns):
        for c in bl_p.columns:
            if bl_p[c].dtype.kind == "f":
                d = float(np.nanmax(np.abs(bl_p[c].values - bl_r[c].values)))
                check(f"[{tag}] blocks.{c} agrees", d < 1e-9, f"max|diff|={d:.2e}")
            else:
                same = bool((bl_p[c].values == bl_r[c].values).all())
                check(f"[{tag}] blocks.{c} identical", same)

    cl_p = pd.read_csv(os.path.join(py_dir, "openSea_cluster_effects.csv.gz"))
    cl_r = pd.read_csv(os.path.join(r_dir, "openSea_cluster_effects.csv.gz"))
    check(f"[{tag}] cluster table has the same shape and columns",
          cl_p.shape == cl_r.shape and list(cl_p.columns) == list(cl_r.columns),
          f"{cl_p.shape} vs {cl_r.shape}")
    if cl_p.shape == cl_r.shape and list(cl_p.columns) == list(cl_r.columns):
        worst, where = 0.0, ""
        for c in cl_p.columns:
            if cl_p[c].dtype.kind == "f":
                d = float(np.nanmax(np.abs(cl_p[c].values - cl_r[c].values)))
                if d > worst:
                    worst, where = d, c
            else:
                check(f"[{tag}] clusters.{c} identical",
                      bool((cl_p[c].values == cl_r[c].values).all()))
        check(f"[{tag}] every cluster-level numeric column agrees",
              worst < 1e-8, f"worst {where} max|diff|={worst:.2e}")

    jp = strict_json(os.path.join(py_dir, "hsmm_params.json"))
    jr = strict_json(os.path.join(r_dir, "hsmm_params.json"))
    for k in ("n_iter", "neutral_state", "n_clusters", "n_openSea_probes",
              "n_blocks", "state_labels", "legacy_fixed_collapse"):
        check(f"[{tag}] params.{k} agrees", jp[k] == jr[k], f"{jp[k]!r} vs {jr[k]!r}")
    for k in ("loglik", "length_scale", "cross_array_r",
              "cross_array_sign_concordance"):
        a, b = jp[k], jr[k]
        if a is None or b is None:
            check(f"[{tag}] params.{k} null in both", a is None and b is None,
                  f"{a!r} vs {b!r}")
        else:
            d = abs(a - b) / max(1.0, abs(a))
            check(f"[{tag}] params.{k} agrees", d < 1e-9, f"rel diff {d:.2e}")
    for k in ("mu", "sigma", "pi"):
        d = float(np.max(np.abs(np.array(jp[k], float) - np.array(jr[k], float))))
        check(f"[{tag}] params.{k} agrees", d < 1e-9, f"max|diff|={d:.2e}")
    return jp


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        raise SystemExit(__doc__.strip().splitlines()[-1])
    stage_dir = os.path.abspath(argv[0])

    with tempfile.TemporaryDirectory() as tmp:
        # one covariate: both array arms keep residual degrees of freedom, so
        # the per-array comparison runs
        for impl, cmd in (("py", [sys.executable, os.path.join(BIN, "05_blocks_hsmm.py")]),
                          ("r", ["Rscript", os.path.join(BIN, "05_blocks_hsmm.R")])):
            run(cmd, os.path.join(tmp, f"cross_{impl}"), stage_dir, "cov1")
        jp = compare(os.path.join(tmp, "cross_py"), os.path.join(tmp, "cross_r"),
                     "cross-array")
        check("[cross-array] the per-array comparison actually ran",
              jp["cross_array_r"] is not None)

        # three covariates: the 450K arm of this design has no residual degrees
        # of freedom, so the identification guard fires and the cross-array
        # fields are undefined -- the case that used to emit invalid JSON
        for impl, cmd in (("py", [sys.executable, os.path.join(BIN, "05_blocks_hsmm.py")]),
                          ("r", ["Rscript", os.path.join(BIN, "05_blocks_hsmm.R")])):
            run(cmd, os.path.join(tmp, f"skip_{impl}"), stage_dir, "cov1,cov2,cov3")
        jp = compare(os.path.join(tmp, "skip_py"), os.path.join(tmp, "skip_r"),
                     "arm skipped")
        check("[arm skipped] cross-array fields are null, not NaN",
              jp["cross_array_r"] is None
              and jp["cross_array_sign_concordance"] is None)

    print("\n" + ("ALL STAGE-05 CHECKS PASSED" if not FAILED
                  else f"{len(FAILED)} FAILED: " + "; ".join(FAILED)))
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main())
