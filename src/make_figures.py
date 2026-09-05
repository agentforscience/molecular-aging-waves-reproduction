"""Generate all figures from the saved result tables.

Design rules followed (see the dataviz reference):
  * categorical hues assigned in a FIXED order, never cycled, and never more
    than three in one panel -- more categories become small multiples;
  * NO dual-axis charts. Where a DE-SWAN curve must be read against its
    per-window group sizes, the group sizes go in a stacked panel sharing the
    x axis, not on a second y scale;
  * one hue light-to-dark for sequential encodings; recessive grid and axes;
  * a legend whenever two or more series are present; text in ink, not in the
    series colour;
  * thin marks, labelled axes with units, and a caption-ready title.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES, TABLES, get_logger

LOG = get_logger("make_figures")

# Validated categorical palette, light mode, fixed order.
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"     # blue, orange, aqua
C4, C5 = "#eda100", "#4a3aa7"                      # yellow, violet
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#b8b7b2"
GRID = dict(color="#e6e6e3", linewidth=0.7)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 8.5, "axes.titlesize": 9.5, "axes.labelsize": 8.5,
    "axes.edgecolor": "#c9c8c4", "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 7.5,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


def _grid(ax, axis="y"):
    ax.grid(axis=axis, **GRID)
    ax.set_axisbelow(True)


def _t(name):
    p = os.path.join(TABLES, name)
    return pd.read_csv(p) if os.path.exists(p) else None


def _save(fig, name):
    path = os.path.join(FIGURES, name)
    fig.savefig(path)
    plt.close(fig)
    LOG.info(f"wrote {path}")


# ---------------------------------------------------------------- Figure 1
def fig1_reproduction():
    """All 10 iPOP layers: published pipeline vs the same test without LOESS."""
    d = _t("e1_curves.csv")
    if d is None:
        return
    layers = [l for l in d.layer.unique()]
    fig, axes = plt.subplots(2, 5, figsize=(13.5, 5.8), sharex=True,
                             gridspec_kw=dict(hspace=.52, wspace=.30))
    for ax, layer in zip(axes.ravel(), layers):
        s = d[d.layer == layer]
        a = s[s.arm == "loess_deswan"].sort_values("midpoint")
        b = s[s.arm == "deswan_only"].sort_values("midpoint")
        p = a.n_significant.max() / max(a.frac_significant.max(), 1e-12)
        ax.plot(a.midpoint, 100 * a.frac_significant, color=C1, lw=1.8,
                label="LOESS + DE-SWAN (published)")
        ax.plot(b.midpoint, 100 * b.frac_significant, color=C2, lw=1.8,
                label="DE-SWAN only (no LOESS)")
        ax.axvline(44, color=MUTED, lw=0.8, ls=":")
        ax.axvline(60, color=MUTED, lw=0.8, ls=":")
        ax.set_ylim(-3, 100)
        ax.set_title(f"{layer.replace('plasma_','').replace('_',' ')}\n"
                     f"{int(round(p)):,} variables", fontsize=8)
        _grid(ax)
    for ax in axes[1]:
        ax.set_xlabel("window centre (years)")
    for ax in axes[:, 0]:
        ax.set_ylabel("% of variables significant\n(BH q < 0.05)")
    axes[0, 0].legend(loc="lower left", bbox_to_anchor=(0, 1.28), ncol=1)
    fig.suptitle("The published pipeline calls most molecules significant at every age; "
                 "removing LOESS leaves almost nothing\n"
                 "Dotted lines mark the published transition ages 44 and 60",
                 y=1.10, fontsize=10.5, ha="center")
    _save(fig, "fig1_reproduction_all_layers.png")


# ---------------------------------------------------------------- Figure 2
def fig2_decomposition():
    """Smoothing vs densification, and the instability of the crest age."""
    comp = _t("e2_components.csv")
    one = _t("e2_onefactor.csv")
    if comp is None:
        return
    order = ["raw", "interp_only", "loess_at_ages", "loess_grid"]
    labels = ["raw\n(neither)", "interpolate\n(densify only)",
              "LOESS at ages\n(smooth only)", "LOESS on grid\n(PUBLISHED)"]
    layers = comp.layer.unique()

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2),
                             gridspec_kw=dict(width_ratios=[1.15, 1], wspace=.34))

    # (a) peak % significant per representation -- sequential single hue
    ax = axes[0]
    x = np.arange(len(order))
    w = 0.78 / len(layers)
    shades = plt.cm.Blues(np.linspace(0.40, 0.92, len(layers)))
    for i, (layer, col) in enumerate(zip(layers, shades)):
        s = comp[comp.layer == layer].set_index("representation").reindex(order)
        ax.bar(x + i * w - 0.39 + w / 2, 100 * s.peak_frac, w * 0.9, color=col,
               label=layer.replace("plasma_", "").replace("_", " "))
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("peak % of variables significant")
    ax.set_title("(a) Smoothing, not densification, manufactures the signal")
    ax.legend(ncol=2, loc="upper left", fontsize=6.8)
    _grid(ax)

    # (b) crest age across every ablation -- one row per layer, dots per setting
    ax = axes[1]
    if one is not None:
        sub = one[one.arm == "loess_grid"].dropna(subset=["crest_age"])
        ly = list(sub.layer.unique())
        for i, layer in enumerate(ly):
            v = sub[sub.layer == layer].crest_age.to_numpy(float)
            ax.scatter(v, np.full(v.size, i) + np.random.default_rng(i).normal(0, .06, v.size),
                       s=22, color=C1, alpha=.55, edgecolor="white", linewidth=.5)
            ax.plot([v.min(), v.max()], [i, i], color=C1, lw=1.0, alpha=.35, zorder=0)
        ax.set_yticks(range(len(ly)))
        ax.set_yticklabels([l.replace("plasma_", "").replace("metabolomics_metabolite",
                                                             "metabolites").replace("_", " ")
                            for l in ly], fontsize=7.5)
        ax.set_ylim(-0.7, len(ly) - 0.3)
        ax.axvline(44, color=C2, lw=1.2, ls="--")
        ax.axvline(60, color=C2, lw=1.2, ls="--")
        ax.text(44, -0.62, " published 44", color=C2, fontsize=7, va="bottom")
        ax.text(60, -0.62, " 60", color=C2, fontsize=7, va="bottom")
        ax.set_xlabel("crest age (years)")
        ax.set_title("(b) The crest age moves with arbitrary analysis choices\n"
                     "each dot = one setting of test / window scheme / span / adjustment")
        _grid(ax, axis="x")
    _save(fig, "fig2_decomposition.png")


# ---------------------------------------------------------------- Figure 3
def fig3_permutation():
    """Observed curve against the valid and the invalid permutation null."""
    d = _t("e3_null_curves.csv")
    pv = _t("e3_pvalues.csv")
    if d is None:
        return
    layers = list(d.layer.unique())
    ncol = min(3, len(layers))
    nrow = int(np.ceil(len(layers) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 3.3 * nrow), squeeze=False)
    for ax, layer in zip(axes.ravel(), layers):
        s = d[d.layer == layer].sort_values("midpoint")
        p = s.n_variables.iloc[0]
        ax.fill_between(s.midpoint, 100 * s.null_valid_q025 / p, 100 * s.null_valid_q975 / p,
                        color=C3, alpha=.22, lw=0,
                        label="valid null: shuffle ages BEFORE LOESS (95%)")
        ax.plot(s.midpoint, 100 * s.null_valid_mean / p, color=C3, lw=1.4, ls="--")
        ax.plot(s.midpoint, 100 * s.null_invalid_mean / p, color=C2, lw=1.6,
                label="invalid null: shuffle AFTER LOESS (Shen et al.)")
        ax.plot(s.midpoint, 100 * s.observed / p, color=C1, lw=2.0, label="observed data")
        ax.axvline(44, color=MUTED, lw=0.8, ls=":")
        ax.axvline(60, color=MUTED, lw=0.8, ls=":")
        ax.set_ylim(-3, 100)
        row = pv[pv.layer == layer].iloc[0] if pv is not None and (pv.layer == layer).any() else None
        ttl = layer.replace("plasma_", "").replace("_", " ")
        if row is not None:
            ttl += (f"\np = {row.p_max_valid:.2f} (valid)   "
                    f"p = {row.p_max_invalid:.3f} (invalid)")
        ax.set_title(ttl, fontsize=8.5)
        _grid(ax)
    for ax in axes.ravel()[len(layers):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("window centre (years)")
    for ax in axes[:, 0]:
        ax.set_ylabel("% significant (BH q < 0.05)")
    axes[0, 0].legend(loc="lower left", bbox_to_anchor=(0, 1.32), fontsize=7)
    fig.suptitle("The observed curve is indistinguishable from the VALID permutation null, "
                 "and far above the INVALID one\n"
                 "The published significance rests entirely on the ordering of the shuffle",
                 y=1.03, fontsize=10.5)
    fig.tight_layout()
    _save(fig, "fig3_permutation_nulls.png")


# ---------------------------------------------------------------- Figure 4
def fig4_nullsim():
    """Crests from pure noise; curve shape set by the age distribution; DE-SWAN
    on strictly linear data."""
    # Panel (a) uses the widest admissible range of window centres (36-65), not
    # the published 40-65 plotting range: within the narrower range the curve's
    # maximum falls on the left boundary and would read as "a crest at 40".
    a = _t("e4a_null_curves_extended.csv")
    if a is None:
        a = _t("e4a_null_curves.csv")
    b = _t("e4b_age_dists.csv")
    dlin = _t("e4d_linear.csv")
    if a is None:
        return
    fig = plt.figure(figsize=(13.5, 6.9))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=.62, wspace=.32)

    # (a) mean null curve at the real iPOP ages
    ax = fig.add_subplot(gs[0, :2])
    ax.fill_between(a.midpoint, 100 * a.q025 / a.n_molecules, 100 * a.q975 / a.n_molecules,
                    color=C1, alpha=.20, lw=0, label="95% across simulations")
    ax.plot(a.midpoint, 100 * a.mean_significant / a.n_molecules, color=C1, lw=2.2,
            label="mean over simulations")
    ax.axvspan(40, a.midpoint.max(), color=MUTED, alpha=.13, lw=0,
               label="range plotted by Shen et al. (40–65)")
    ax.axvline(44, color=C2, lw=1.2, ls="--")
    ax.axvline(60, color=C2, lw=1.2, ls="--")
    lo, hi = ax.get_ylim()
    ax.text(44, hi, " 44 ", color=C2, fontsize=8, va="top")
    ax.text(60, hi, " 60 ", color=C2, fontsize=8, va="top")
    ax.set_xlabel("window centre (years)")
    ax.set_ylabel("% significant (BH q < 0.05)")
    ax.set_title("(a) 1,000 molecules of PURE NOISE, no age relationship at all, at the "
                 "real iPOP ages.\nThe curve has an interior local maximum at 61, and "
                 "inside the published plotting\nwindow its maximum sits at the left "
                 "boundary — reading as a crest in the low 40s.", fontsize=8.8)
    ax.legend(loc="lower left", fontsize=7)
    _grid(ax)

    # (b) small multiples: curve shape by age distribution (no colour cycling)
    if b is not None:
        dists = list(b.age_dist.unique())
        for i, dist in enumerate(dists):
            ax = fig.add_subplot(gs[0, 2 + i % 2] if i < 2 else gs[1, 2 + (i - 2) % 2])
            s = b[b.age_dist == dist].sort_values("midpoint")
            ax.plot(s.midpoint, 100 * s.mean_significant / s.n_molecules, color=C1, lw=1.8)
            ax.fill_between(s.midpoint, 100 * s.q025 / s.n_molecules,
                            100 * s.q975 / s.n_molecules, color=C1, alpha=.18, lw=0)
            ax.axvline(44, color=MUTED, lw=.8, ls=":")
            ax.axvline(60, color=MUTED, lw=.8, ls=":")
            ax.set_title(f"{'(b) ' if i == 0 else ''}null data, {dist} ages", fontsize=8.5)
            ax.set_xlabel("window centre (years)")
            if i % 2 == 0:
                ax.set_ylabel("% significant")
            _grid(ax)

    # (c) DE-SWAN on strictly LINEAR data: counts and group sizes in STACKED
    # panels sharing x (never a second y axis)
    if dlin is not None:
        inner = gs[1, :2].subgridspec(2, 1, height_ratios=[2, 1], hspace=.12)
        axc = fig.add_subplot(inner[0])
        axn = fig.add_subplot(inner[1], sharex=axc)
        # Scenarios differ by an order of magnitude in absolute count, and the
        # claim here is about SHAPE, so each curve is scaled to its own peak.
        scen = ["homoskedastic", "variance_increasing", "outlier_cluster_55"]
        for name, col in zip(scen, (C1, C2, C3)):
            s = dlin[dlin.scenario == name].sort_values("midpoint")
            if s.empty:
                continue
            y = s.mean_significant.to_numpy(float)
            axc.plot(s.midpoint, 100 * y / max(y.max(), 1e-9), color=col, lw=1.8,
                     label=f"{name.replace('_', ' ')} (peak {y.max():.0f}/1000)")
        s0 = dlin[dlin.scenario == "homoskedastic"].sort_values("midpoint")
        hm = 2 / (1 / s0.n_young.clip(lower=1) + 1 / s0.n_old.clip(lower=1))
        axn.bar(s0.midpoint, hm, width=.8, color=MUTED)
        axn.set_ylabel("harmonic\nmean n", fontsize=7.5)
        axn.set_xlabel("window centre (years)")
        axc.set_ylabel("% of that scenario's\npeak count")
        axc.tick_params(labelbottom=False)
        axc.set_title("(c) DE-SWAN alone (no LOESS) on STRICTLY LINEAR data:\n"
                      "the 'wave' tracks per-window sample size, not biology", fontsize=9)
        axc.legend(loc="upper left", fontsize=6.8)
        _grid(axc)
        _grid(axn)
    _save(fig, "fig4_null_simulations.png")


# ---------------------------------------------------------------- Figure 5
def _short(name):
    return (name.replace("plasma_", "")
                .replace("metabolomics_metabolite", "metabolites")
                .replace("_", " "))


def fig5_corrected():
    """iPOP under corrected inference: what survives, how much of it is one
    subject, and how uncertain the published crest actually is."""
    bt = _t("e5_crest_bootstrap.csv")
    summ = _t("e5_tests_summary.csv")
    infl = _t("e7_influence.csv")
    cyt = _t("e7b_cytokine_sensitivity.csv")
    if summ is None:
        return
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2),
                             gridspec_kw=dict(wspace=.36))

    # (a) counts surviving each of the three separated claims
    ax = axes[0]
    lay = summ.layer.tolist()
    x = np.arange(len(lay))
    for i, (col, key, lab) in enumerate([
            (C1, "frac_linear_bh", "any age association"),
            (C2, "frac_nonlinear_bh", "nonlinear"),
            (C3, "frac_segmented_bh", "discrete transition")]):
        ax.bar(x + (i - 1) * .27, 100 * summ[key].fillna(0), .25, color=col, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels([_short(l) for l in lay], fontsize=7, rotation=32, ha="right")
    ax.set_ylabel("% of variables surviving\n(permutation-calibrated, BH q < 0.05)")
    ax.set_title("(a) iPOP under corrected inference", fontsize=9)
    ax.legend(fontsize=7)
    _grid(ax)

    # (b) the one positive finding, and the single subject it rests on
    ax = axes[1]
    if infl is not None and (infl.layer == "plasma_cytokine").any():
        s = infl[infl.layer == "plasma_cytokine"]
        ax.scatter(s.dropped_subject_age, 100 * s.supF_retention, s=26, color=C1,
                   alpha=.65, edgecolor="white", linewidth=.5,
                   label="drop one subject, refit")
        w = s.loc[s.supF_retention.idxmin()]
        ax.scatter([w.dropped_subject_age], [100 * w.supF_retention], s=90, color=C2,
                   zorder=4, label=f"the subject aged {w.dropped_subject_age:.0f}")
        ax.annotate(f"{100*w.supF_retention:.0f}% of the statistic remains",
                    xy=(w.dropped_subject_age, 100 * w.supF_retention),
                    xytext=(w.dropped_subject_age + 6, 100 * w.supF_retention + 28),
                    fontsize=7, color=INK2,
                    arrowprops=dict(arrowstyle="->", color=MUTED, lw=.9))
        ax.axhline(100, color=MUTED, lw=.9, ls=":")
        ax.set_xlabel("age of the subject removed (years)")
        ax.set_ylabel("median sup-F retained (%)")
        ax.set_title("(b) The one positive finding — 35/66 cytokines\n"
                     "'transitioning' at 41 — is a single subject", fontsize=9)
        ax.legend(fontsize=7, loc="lower right")
        _grid(ax)

    # (c) bootstrap CI of the published crest age
    ax = axes[2]
    if bt is not None:
        y = np.arange(len(bt))
        ax.hlines(y, bt.ci_lo, bt.ci_hi, color=C1, lw=5, alpha=.30,
                  label="bootstrap 95% CI")
        ax.scatter(bt.crest_point_estimate, y, s=34, color=C1, zorder=3,
                   label="published-style point estimate")
        ax.axvline(44, color=C2, lw=1.2, ls="--")
        ax.axvline(60, color=C2, lw=1.2, ls="--")
        ax.set_yticks(y)
        ax.set_yticklabels([_short(l) for l in bt.layer], fontsize=7.5)
        ax.set_ylim(-0.8, len(bt) - 0.2)
        ax.set_xlabel("crest age (years)")
        ax.set_title("(c) Bootstrap 95% CI of the published crest", fontsize=9)
        ax.legend(fontsize=7, loc="lower right")
        _grid(ax, axis="x")
    _save(fig, "fig5_corrected_ipop.png")


def fig5b_powerequalised():
    """Power-equalised DE-SWAN curves against their permutation nulls."""
    pe = _t("e5_powerequal_deswan.csv")
    if pe is None:
        return
    layers = list(pe.layer.unique())
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 6.0), squeeze=False,
                             gridspec_kw=dict(hspace=.62, wspace=.30))
    for ax, layer in zip(axes.ravel(), layers):
        s = pe[pe.layer == layer].sort_values("midpoint")
        inner = ax.get_subplotspec().subgridspec(2, 1, height_ratios=[2.2, 1], hspace=.12)
        ax.remove()
        axc = fig.add_subplot(inner[0])
        axn = fig.add_subplot(inner[1], sharex=axc)
        axc.fill_between(s.midpoint, 0, s.null_q975, color=C3, alpha=.22, lw=0,
                         label="97.5th pct of permutation null")
        axc.plot(s.midpoint, s.n_significant, color=C1, lw=2.0, label="observed")
        axc.axvline(44, color=MUTED, lw=.8, ls=":")
        axc.axvline(60, color=MUTED, lw=.8, ls=":")
        axc.set_title(_short(layer), fontsize=8.8)
        axc.tick_params(labelbottom=False)
        axc.set_ylabel("variables\nsignificant", fontsize=7.5)
        # Group sizes are constant by construction; the window WIDTH varies instead.
        axn.bar(s.midpoint, s.window_width_years, width=.8, color=MUTED)
        axn.set_ylabel("window\nwidth (yr)", fontsize=7)
        axn.set_xlabel("window centre (years)", fontsize=8)
        _grid(axc)
        _grid(axn)
        if layer == layers[0]:
            axc.legend(fontsize=6.8, loc="upper right")
    fig.suptitle("Power-equalised DE-SWAN: group sizes held constant at (k, k) by construction, "
                 "so window WIDTH varies instead (lower panels).\n"
                 "Counts are zero at almost every centre in every layer; the single cytokine "
                 "analyte at 47 gives family-wise p = 0.25.\n"
                 "Centres below ~45 are undefined because fewer than k subjects lie below them — "
                 "power equalisation costs edge coverage.",
                 y=1.05, fontsize=9.5)
    _save(fig, "fig5b_power_equalised.png")


# ---------------------------------------------------------------- Figure 6
def fig6_gse():
    """GSE40279 (n=656): does a well-powered cohort localise change at 44/60?"""
    cur = _t("e6_curves.csv")
    bpd = _t("e6_breakpoint_dist.csv")
    if cur is None:
        return
    fig = plt.figure(figsize=(13.5, 6.4))
    gs = fig.add_gridspec(2, 3, hspace=.55, wspace=.28)

    sel = "top_variance"
    # (a)-(c): the three DE-SWAN arms, adjusted vs unadjusted, with the group
    # sizes in a stacked sub-panel (never a second y axis)
    arms = [("A_loess_deswan", "(a) LOESS + DE-SWAN (published pipeline)"),
            ("B_deswan_only", "(b) DE-SWAN only, fixed-width windows"),
            ("C_power_equalised", "(c) DE-SWAN, power-equalised windows")]
    for i, (arm, title) in enumerate(arms):
        inner = gs[0, i].subgridspec(2, 1, height_ratios=[2.4, 1], hspace=.10)
        axc = fig.add_subplot(inner[0])
        axn = fig.add_subplot(inner[1], sharex=axc)
        for adj, col, lab in ((False, C1, "unadjusted"), (True, C2, "batch-adjusted")):
            s = cur[(cur.selection == sel) & (cur.arm == arm) & (cur.adjusted == adj)]
            s = s.sort_values("midpoint")
            if s.empty:
                continue
            axc.plot(s.midpoint, 100 * s.frac_significant, color=col, lw=1.8, label=lab)
            if s.null_q975.notna().any():
                axc.plot(s.midpoint, 100 * s.null_q975 / (s.n_significant.max() /
                         max(s.frac_significant.max(), 1e-9)), color=col, lw=.9, ls=":")
        s0 = cur[(cur.selection == sel) & (cur.arm == arm) & (~cur.adjusted)].sort_values("midpoint")
        if not s0.empty:
            hm = 2 / (1 / s0.n_young.clip(lower=1) + 1 / s0.n_old.clip(lower=1))
            axn.bar(s0.midpoint, hm, width=.8, color=MUTED)
        axn.set_ylabel("harmonic\nmean n", fontsize=7)
        axn.set_xlabel("window centre (years)")
        axc.tick_params(labelbottom=False)
        axc.axvline(44, color=MUTED, lw=.8, ls=":")
        axc.axvline(60, color=MUTED, lw=.8, ls=":")
        axc.set_title(title, fontsize=8.8)
        if i == 0:
            axc.set_ylabel("% CpGs significant")
            axc.legend(fontsize=7, loc="lower left")
        _grid(axc)
        _grid(axn)

    # (d) breakpoint distributions, observed vs null
    ax = fig.add_subplot(gs[1, :2])
    if bpd is not None:
        for adj, col, lab in ((False, C1, "observed, unadjusted"),
                              (True, C2, "observed, batch-adjusted")):
            s = bpd[(bpd.selection == sel) & (bpd.adjusted == adj)].sort_values("breakpoint")
            if s.empty or s.observed_frac.isna().all():
                continue
            ax.plot(s.breakpoint, 100 * s.observed_frac, color=col, lw=1.8, label=lab)
        s = bpd[(bpd.selection == sel) & (~bpd.adjusted)].sort_values("breakpoint")
        ax.plot(s.breakpoint, 100 * s.null_frac, color=C3, lw=1.4, ls="--",
                label="permutation null")
        ax.axvline(44, color=MUTED, lw=.9, ls=":")
        ax.axvline(60, color=MUTED, lw=.9, ls=":")
        ax.set_xlabel("estimated breakpoint age (years)")
        ax.set_ylabel("% of CpGs with a\nsignificant transition here")
        ax.set_title("(d) Where do transitions actually localise in a well-powered cohort?",
                     fontsize=9)
        ax.legend(fontsize=7)
        _grid(ax)

    # (e) selection sensitivity: top-variance vs random CpGs
    ax = fig.add_subplot(gs[1, 2])
    for s_, col, lab in (("top_variance", C1, "top-variance CpGs"),
                         ("random", C2, "random CpGs")):
        s = cur[(cur.selection == s_) & (cur.arm == "C_power_equalised") &
                (cur.adjusted)].sort_values("midpoint")
        if s.empty:
            continue
        ax.plot(s.midpoint, 100 * s.frac_significant, color=col, lw=1.8, label=lab)
    ax.axvline(44, color=MUTED, lw=.8, ls=":")
    ax.axvline(60, color=MUTED, lw=.8, ls=":")
    ax.set_xlabel("window centre (years)")
    ax.set_ylabel("% CpGs significant")
    ax.set_title("(e) CpG selection sensitivity\n(power-equalised, adjusted)", fontsize=8.8)
    ax.legend(fontsize=7)
    _grid(ax)

    fig.suptitle("GSE40279 whole-blood methylation, n = 656, ages 19–101",
                 y=.99, fontsize=10.5)
    _save(fig, "fig6_gse40279.png")


# ---------------------------------------------------------------- Figure 7
def fig7_headline():
    """One-panel headline: achieved false discovery proportion vs nominal."""
    pv = _t("e3_pvalues.csv")
    if pv is None:
        return
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    y = np.arange(len(pv))
    ax.barh(y, 100 * pv.achieved_fdp_mean, .62, color=C1, label="achieved (measured on null data)")
    ax.axvline(100 * pv.nominal_fdr.iloc[0], color=C2, lw=1.8, ls="--",
               label="nominal FDR (5%)")
    for i, v in enumerate(pv.achieved_fdp_mean):
        ax.text(100 * v + 1.2, i, f"{100*v:.0f}%", va="center", fontsize=7.5, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels([l.replace("plasma_", "").replace("_", " ") for l in pv.layer], fontsize=7.5)
    ax.set_xlabel("false discovery proportion (%)")
    ax.set_title("The published pipeline's achieved error rate is ~13× its nominal rate\n"
                 "measured on data where every null hypothesis is true by construction",
                 fontsize=9.5)
    ax.legend(fontsize=7.5, loc="lower right")
    _grid(ax, axis="x")
    _save(fig, "fig7_achieved_fdr.png")


def main():
    for fn in (fig1_reproduction, fig2_decomposition, fig3_permutation, fig4_nullsim,
               fig5_corrected, fig5b_powerequalised, fig6_gse, fig7_headline):
        try:
            fn()
        except Exception as exc:
            LOG.error(f"{fn.__name__} FAILED: {exc!r}")
            import traceback
            LOG.error(traceback.format_exc())


if __name__ == "__main__":
    main()
