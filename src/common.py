"""Shared I/O, configuration and logging for all experiments."""
from __future__ import annotations

import json
import logging
import os
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CS = os.path.join(ROOT, "datasets", "ipop_cross_section")
RESULTS = os.path.join(ROOT, "results")
TABLES = os.path.join(RESULTS, "tables")
FIGURES = os.path.join(ROOT, "figures")
LOGS = os.path.join(ROOT, "logs")
for _d in (RESULTS, TABLES, FIGURES, LOGS):
    os.makedirs(_d, exist_ok=True)

SEED = 42

# The 10 iPOP omics layers, ordered so that the two that drive the published
# crests (transcriptome, proteomics) come first.
LAYERS = [
    "plasma_transcriptome",
    "plasma_proteomics",
    "plasma_metabolomics_metabolite",
    "plasma_lipidomics",
    "plasma_cytokine",
    "clinical_test",
    "gut_microbiome",
    "oral_microbiome",
    "skin_microbiome",
    "nasal_microbiome",
]

# Layers small enough to run every ablation and permutation study on in full.
FAST_LAYERS = ["plasma_proteomics", "plasma_cytokine", "clinical_test",
               "plasma_metabolomics_metabolite", "plasma_lipidomics"]


def get_logger(name: str) -> logging.Logger:
    log = logging.getLogger(name)
    if log.handlers:
        return log
    log.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%H:%M:%S")
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    log.addHandler(sh)
    fh = logging.FileHandler(os.path.join(LOGS, f"{name}.log"), mode="w")
    fh.setFormatter(fmt)
    log.addHandler(fh)
    return log


def set_seed(seed: int = SEED) -> np.random.Generator:
    """Seed global RNG state and return a fresh generator for local use."""
    import random
    random.seed(seed)
    np.random.seed(seed)
    return np.random.default_rng(seed)


def load_layer(layer: str, log_transform: str = "auto"):
    """Load an iPOP cross-sectional layer.

    Returns
    -------
    Y : (p, n) ndarray, variables x subjects
    ages : (n,) subject ages
    info : DataFrame of subject covariates, aligned to the columns of Y
    var_ids : Index of variable names

    Microbiome layers are relative-abundance counts with extreme skew and many
    structural zeros; `log_transform="auto"` applies log1p to those layers only.
    This matches the original pipeline's treatment (the ipop_aging code
    log-transforms microbiome data before smoothing) and is recorded in the
    per-layer metadata written by each experiment.
    """
    expr = pd.read_csv(os.path.join(CS, f"{layer}__expression.csv"), index_col=0)
    info = pd.read_csv(os.path.join(CS, f"{layer}__sample_info.csv"))
    info = info.set_index("subject_id").reindex(expr.columns)
    ages = info["adjusted_age"].to_numpy(float)

    Y = expr.to_numpy(float)
    transformed = False
    if log_transform == "auto" and "microbiome" in layer:
        Y = np.log1p(Y)
        transformed = True
    elif log_transform is True:
        Y = np.log1p(Y)
        transformed = True

    # Drop variables that are constant or all-missing: no test is defined for them.
    keep = np.isfinite(Y).sum(axis=1) >= 10
    with np.errstate(invalid="ignore"):
        keep &= np.nanstd(Y, axis=1) > 0
    Y, var_ids = Y[keep], expr.index[keep]

    info = info.copy()
    info.attrs["log1p"] = transformed
    return Y, ages, info, var_ids


def complete_rows(Y: np.ndarray):
    """Mask of variables with no missing values.

    The least-squares machinery in `corrected.py` works on whole matrices and
    propagates NaN, so variables with any missing observation must be excluded
    explicitly rather than silently producing NaN statistics. Only
    `clinical_test` is affected among the iPOP layers (2.9% missing).
    """
    return np.isfinite(np.asarray(Y, float)).all(axis=1)


def zscore_rows(Y: np.ndarray) -> np.ndarray:
    """Standardise each variable to mean 0 / sd 1 across subjects.

    The published pipeline z-scores before smoothing so that clusters and
    trajectories are comparable across molecules. Rank-based tests are invariant
    to this, but effect sizes and spline fits are not, so it is applied
    explicitly and identically in every arm.
    """
    mu = np.nanmean(Y, axis=1, keepdims=True)
    sd = np.nanstd(Y, axis=1, keepdims=True)
    sd = np.where(sd > 0, sd, 1.0)
    return (Y - mu) / sd


def regress_out(Y: np.ndarray, info: pd.DataFrame, cols=("Gender", "Ethnicity", "IRIS")):
    """Residualise each variable on the available confounders.

    Shen et al. regress out BMI, sex, IRIS status and ethnicity and analyse the
    residuals. BMI is not present in the recovered sample tables, so this uses
    the three that are. Used only in the confounder-adjustment sensitivity
    analysis, never silently in the main arms.
    """
    present = [c for c in cols if c in info.columns and info[c].notna().any()]
    if not present:
        return Y, []
    X = pd.get_dummies(info[present].astype(str), drop_first=True).to_numpy(float)
    X = np.column_stack([np.ones(X.shape[0]), X])
    R = np.empty_like(Y)
    for i in range(Y.shape[0]):
        y = Y[i]
        m = np.isfinite(y)
        beta, *_ = np.linalg.lstsq(X[m], y[m], rcond=None)
        R[i] = y - X @ beta
        R[i, ~m] = np.nan
    return R, present


def save_table(df: pd.DataFrame, name: str, log=None) -> str:
    path = os.path.join(TABLES, name)
    df.to_csv(path, index=False)
    if log:
        log.info(f"wrote {path}  ({len(df)} rows)")
    return path


def save_json(obj, name: str, log=None) -> str:
    path = os.path.join(RESULTS, name)

    def _default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o) if np.isfinite(o) else None
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        raise TypeError(type(o))

    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=_default)
    if log:
        log.info(f"wrote {path}")
    return path


def env_report() -> dict:
    import scipy
    import sklearn
    import statsmodels
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "statsmodels": statsmodels.__version__,
        "scikit-learn": sklearn.__version__,
        "seed": SEED,
    }
