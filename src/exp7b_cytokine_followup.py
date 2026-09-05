"""E7b -- Follow-up on the single positive finding: the iPOP cytokine "transition".

E5 found 35/66 cytokine analytes with a significant breakpoint at age 41, and E7
showed that a rank-inverse-normal transform (RINT) reduces but does not abolish
it: 4-7 analytes survive with breakpoints at 41-44, i.e. within 2 years of the
published transition age of 44. Taken at face value that would be the one place
in this project where a corrected analysis of iPOP localises change near a
published crest -- so it is exactly the finding that must be attacked hardest.

E7's influence analysis pointed at a single subject aged 25.9 whose removal left
2.8% of the median sup-F statistic. This script completes the argument by
re-running the full corrected battery on nested subject subsets and reporting
what survives. The RINT is retained throughout, so any result here is already
robust to the *magnitude* of an outlier and can only be driven by its rank
position.

Outputs
-------
results/tables/e7b_cytokine_sensitivity.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import env_report, get_logger, load_layer, save_json, save_table, set_seed
from corrected import (linear_and_spline_tests, perm_pvalues,
                       permutation_null_supF, segmented_supF)
from fastdeswan import bh

LOG = get_logger("exp7b_cytokine_followup")

PSI_GRID = np.arange(35.0, 70.5, 1.0)
N_PERM = 300
ALPHA = 0.05


def rint(Y: np.ndarray) -> np.ndarray:
    from scipy import stats
    n = Y.shape[1]
    r = stats.rankdata(Y, axis=1, method="average")
    return stats.norm.ppf((r - 0.375) / (n + 0.25))


def battery(Yt, ages, rng, label):
    t = linear_and_spline_tests(Yt, ages)
    supF, psi, _ = segmented_supF(Yt, ages, PSI_GRID)
    null = permutation_null_supF(Yt, ages, PSI_GRID, N_PERM, rng)
    q_nl = bh(perm_pvalues(t["F_nonlinear"], null["F_nonlinear"]))
    q_seg = bh(perm_pvalues(supF, null["supF"]))
    sig = q_seg < ALPHA
    row = dict(subset=label, n_subjects=int(ages.size), n_analytes=int(Yt.shape[0]),
               age_min=float(ages.min()), age_max=float(ages.max()),
               n_nonlinear=int((q_nl < ALPHA).sum()), n_segmented=int(sig.sum()),
               breakpoints=";".join(f"{v:.0f}" for v in np.sort(psi[sig])) if sig.sum() else "",
               n_breakpoints_within2_of_44=int(np.sum(np.abs(psi[sig] - 44) <= 2)) if sig.sum() else 0)
    LOG.info(f"{label:44s} n={ages.size:3d}  nonlinear {row['n_nonlinear']:2d}/66  "
             f"segmented {row['n_segmented']:2d}/66  {row['breakpoints']}")
    return row


def main() -> None:
    rng = set_seed()
    Y, ages, info, _ = load_layer("plasma_cytokine")
    order = np.argsort(ages)
    LOG.info(f"cytokine ages, five youngest subjects: {ages[order][:5].round(1)}")

    rows = []
    subsets = [
        ("all subjects (RINT)", np.ones(ages.size, bool)),
        ("drop the single youngest subject (age 25.9)", ages > 27.0),
        ("drop all subjects under 30", ages >= 30.0),
        ("drop all subjects under 35", ages >= 35.0),
        ("drop the single oldest subject", ages < ages.max()),
    ]
    for label, keep in subsets:
        rows.append(battery(rint(Y[:, keep]), ages[keep], rng, label))

    # A matched control: drop one randomly chosen subject, repeatedly. If losing
    # ANY single subject destroyed the finding, that would be a small-sample
    # fragility statement; if only the youngest one does, the finding IS that
    # subject.
    n_seg = []
    for b in range(20):
        keep = np.ones(ages.size, bool)
        keep[rng.integers(0, ages.size)] = False
        r = battery(rint(Y[:, keep]), ages[keep], rng, f"drop one random subject (rep {b+1})")
        n_seg.append(r["n_segmented"])
        rows.append(r)
    LOG.info(f"dropping a RANDOM single subject: median {np.median(n_seg):.0f} analytes "
             f"still significant (range {min(n_seg)}-{max(n_seg)}), versus "
             f"{rows[1]['n_segmented']} when the youngest subject specifically is dropped")

    save_table(pd.DataFrame(rows), "e7b_cytokine_sensitivity.csv", LOG)
    save_json({"n_perm": N_PERM, "alpha": ALPHA, "psi_grid": PSI_GRID.tolist(),
               "env": env_report()}, "e7b_meta.json", LOG)


if __name__ == "__main__":
    main()
