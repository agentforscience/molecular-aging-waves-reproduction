"""Validate src/fastloess.LoessEngine against two references.

Reference 1 (exactness): `code/loess_py/rloess.py`, the pure-Python
reimplementation of R's loess that was itself validated against R output during
the resource-gathering phase. The vectorised engine should agree with it to
floating-point precision, since it computes the same linear functional.

Reference 2 (fidelity): `datasets/ipop_loess_R_csv/`, the LOESS-interpolated
matrices produced by the original authors' R pipeline and shipped in their
analysis repo. Agreement here is close but not bit-exact -- R's default
`surface="interpolate"` uses a kd-tree approximation and CV span ties may break
differently.

Run:  python src/validate_fastloess.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "code", "loess_py"))

from fastloess import DEFAULT_GRID, LoessEngine  # noqa: E402
import rloess  # noqa: E402

LAYER = "plasma_proteomics"
N_EXACT = 40      # variables used for the exact cross-check (reference is slow)
SEED = 42


def load_layer(layer: str):
    expr = pd.read_csv(f"{ROOT}/datasets/ipop_cross_section/{layer}__expression.csv", index_col=0)
    info = pd.read_csv(f"{ROOT}/datasets/ipop_cross_section/{layer}__sample_info.csv")
    info = info.set_index("subject_id").reindex(expr.columns)
    return expr, info["adjusted_age"].to_numpy(float)


def main() -> int:
    rng = np.random.default_rng(SEED)
    expr, ages = load_layer(LAYER)
    Y = expr.to_numpy(float)
    print(f"[data] {LAYER}: {Y.shape[0]} variables x {Y.shape[1]} subjects, "
          f"age {ages.min():.1f}-{ages.max():.1f}")

    # ---------------------------------------------------------------- build
    t0 = time.time()
    eng = LoessEngine(ages, grid=DEFAULT_GRID)
    t_build = time.time() - t0
    print(f"[build] operators for {len(eng.spans)} spans in {t_build:.1f}s")

    t0 = time.time()
    Yhat = eng.smooth(Y)
    t_fast = time.time() - t0
    print(f"[fast ] smoothed {Y.shape[0]} variables in {t_fast:.2f}s "
          f"({1000 * t_fast / Y.shape[0]:.2f} ms/variable)")

    # -------------------------------------------- reference 1: exact agreement
    idx = rng.choice(Y.shape[0], N_EXACT, replace=False)
    t0 = time.time()
    ref, _, ref_spans = rloess.run_loess(Y[idx], ages, grid=DEFAULT_GRID)
    t_ref = time.time() - t0
    print(f"[ref  ] reference implementation: {N_EXACT} variables in {t_ref:.1f}s "
          f"({1000 * t_ref / N_EXACT:.0f} ms/variable) -> speedup "
          f"{(t_ref / N_EXACT) / (t_fast / Y.shape[0]):.0f}x")

    fast_spans, _ = eng.select_spans(Y)
    span_match = float(np.mean(fast_spans[idx] == ref_spans))
    absdiff = np.abs(Yhat[idx] - ref)
    scale = np.maximum(np.abs(ref).mean(axis=1, keepdims=True), 1e-12)
    reldiff = absdiff / scale
    print(f"[chk 1] span agreement {span_match:.3f}; "
          f"max |diff| {absdiff.max():.3e}; median rel diff {np.median(reldiff):.3e}")
    assert span_match == 1.0, "vectorised CV selected different spans than the reference"
    assert absdiff.max() < 1e-8, "vectorised smoother disagrees with the reference fit"

    # ------------------------------------- reference 2: published R LOESS output
    rpath = (f"{ROOT}/datasets/ipop_loess_R_csv/"
             f"{LAYER}_object_cross_section_loess__expression.csv")
    rdf = pd.read_csv(rpath, index_col=0)
    common = expr.index.intersection(rdf.index)
    A = pd.DataFrame(Yhat, index=expr.index).loc[common].to_numpy()
    B = rdf.loc[common].to_numpy(float)
    assert A.shape == B.shape, f"grid mismatch: {A.shape} vs {B.shape}"
    corr = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(A, B)])
    relerr = np.abs(A - B) / np.maximum(np.abs(B).mean(axis=1, keepdims=True), 1e-12)
    print(f"[chk 2] vs published R output on {len(common)} variables: "
          f"median corr {np.median(corr):.4f}, median rel err {np.median(relerr):.2e}")
    assert np.median(corr) > 0.99, "poor agreement with the authors' published R LOESS"

    # ---------------------------------------- permutation invariance sanity check
    # Permuting age labels must not change the *set* of fitted design points, so
    # the operator matrices are reusable. Verify against an explicit re-fit.
    perm = rng.permutation(len(ages))
    eng_perm = LoessEngine(ages[perm], grid=DEFAULT_GRID)
    for s in eng.spans:
        assert np.allclose(eng.L[s], eng_perm.L[s]), "operator depends on label order"
    print("[chk 3] smoother operators invariant to age relabelling: OK")

    print("\nALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
