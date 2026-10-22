"""
The six standard review figures, plus the per-claim panel.

Each figure answers one review question and is built on the corresponding module
of the data-visualization skill:

  F1 design map          -> abundance_and_proportion  (how many? of what unit?)
  F2 raw behind claim    -> distributions             (how does it vary?)
  F3 assumption Q-Q      -> distributions             (what shape is it?)
  F4 claim ledger        -> relationships             (reported vs recomputed)
  F5 null distribution   -> distributions             (where does the observed fall?)
  F6 multiplicity        -> patterns_across_features  (pattern across many tests)
  per-claim panel        -> small multiple of F2/F5 plus the three routes

Shared conventions, from the data-visualization universal principles: axis
titles name the measurement and its units; every panel states its n; a test name
travels with any p-value; verdict is encoded with a colourblind-safe palette and
always also written as text, so colour is never the only channel.

Every function returns a matplotlib Figure and, when `path` is given, saves it.
A figure that cannot be built from the available inputs returns
`not_assessable_figure` -- the fixed figure spine means a missing diagnostic
stays visible instead of silently disappearing.
"""

import os
import sys
from typing import Dict, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Verdict vocabulary shared by the ledger, the report template and the agent.
VERDICT_ORDER = (
    "supported",
    "supported-with-caveat",
    "assumption-dependent",
    "underpowered-resolution",
    "not-reproducible",
    "contradicted",
    "not-assessable",
)

# Okabe-Ito, so figures stay readable in greyscale and to colourblind readers.
VERDICT_COLORS = {
    "supported": "#009E73",
    "supported-with-caveat": "#56B4E9",
    "assumption-dependent": "#E69F00",
    "underpowered-resolution": "#CC79A7",
    "not-reproducible": "#D55E00",
    "contradicted": "#000000",
    "not-assessable": "#999999",
}

FIG_DPI = 300
_FALLBACK_COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]


def _dataviz():
    """
    Import the data-visualization skill's utilities if reachable.

    Returns a dict of its modules, or None. The figures fall back to local
    equivalents when it is absent, so this module works standalone; when both
    skills are present the palettes and annotation helpers are shared.
    """
    skill_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidate = os.path.join(os.path.dirname(skill_root), "data-visualization")
    if os.path.isdir(candidate) and candidate not in sys.path:
        sys.path.insert(0, candidate)
    try:
        from utilities import annotations, color_scales, data_inspection
        return {"annotations": annotations, "color_scales": color_scales,
                "data_inspection": data_inspection}
    except Exception:
        return None


def _palette(n: int):
    dv = _dataviz()
    if dv is not None:
        try:
            return dv["color_scales"].categorical_palette(n, palette="colorblind")
        except Exception:
            pass
    return (_FALLBACK_COLORS * (n // len(_FALLBACK_COLORS) + 1))[:n]


def _finish(fig, path: Optional[str]):
    fig.tight_layout()
    if path:
        parent = os.path.dirname(os.path.abspath(path))
        os.makedirs(parent, exist_ok=True)
        fig.savefig(path, dpi=FIG_DPI, bbox_inches="tight")
    return fig


def not_assessable_figure(title: str, reason: str, path: Optional[str] = None):
    """
    Placeholder for a figure slot the paper or data cannot support.

    Keeping the slot and naming the reason is the point: an absent diagnostic
    and an impossible one look identical in a report that simply omits it.
    """
    fig, ax = plt.subplots(figsize=(6, 2.4))
    ax.axis("off")
    ax.text(0.5, 0.70, title, ha="center", va="center", fontsize=11, fontweight="bold",
            transform=ax.transAxes)
    ax.text(0.5, 0.36, f"Not assessable - {reason}", ha="center", va="center", fontsize=9,
            color="#555555", wrap=True, transform=ax.transAxes)
    ax.add_patch(plt.Rectangle((0.02, 0.05), 0.96, 0.9, fill=False, linestyle="--",
                               edgecolor="#999999", transform=ax.transAxes))
    return _finish(fig, path)


def fig_design_map(data: pd.DataFrame, unit_col: str, group_col: str, value_col: str,
                   independence: Optional[dict] = None, path: Optional[str] = None,
                   unit_label: str = "experimental unit"):
    """
    F1 - What is n, and of what?

    Left: observations per experimental unit, so a reader sees at once whether
    rows are independent units or repeated measures within units.
    Right: units per group, which is the n a test should be using.

    Annotated with the ICC and effective n from
    assumption_checks.design_independence: the gap between the row count and the
    effective n is the finding.
    """
    from .assumption_checks import design_independence

    if independence is None:
        independence = design_independence(data, unit_col, group_col, value_col)

    per_unit = data.groupby(unit_col).size().sort_values(ascending=False)
    unit_group = data.groupby(unit_col)[group_col].first()
    groups = list(pd.unique(data[group_col]))
    cmap = dict(zip(groups, _palette(len(groups))))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [2, 1]})

    ax = axes[0]
    ax.bar(range(len(per_unit)), per_unit.to_numpy(),
           color=[cmap[unit_group[u]] for u in per_unit.index],
           edgecolor="white", linewidth=0.4)
    ax.set_xlabel(f"{unit_label} (one bar per unit, {independence['n_units']} units)")
    ax.set_ylabel(f"observations of {value_col} per unit (count)")
    ax.set_title("Observations nested within units", fontsize=10, loc="left")
    ax.axhline(1, color="#444444", linestyle=":", linewidth=1)
    eff = independence["effective_n"]
    ax.text(0.99, 0.95,
            f"rows:         {independence['n_observations']}\n"
            f"units:        {independence['n_units']}\n"
            f"ICC:          {independence['icc']:.2f}\n"
            f"effective n:  {eff:.1f}" if np.isfinite(eff) else
            f"rows: {independence['n_observations']}\nunits: {independence['n_units']}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="#cccccc"))
    handles = [plt.Line2D([], [], marker="s", linestyle="", color=cmap[g], label=str(g))
               for g in groups]
    # Below the axis: the bars fill the panel, so an in-axes legend collides
    # with them and with the statistics box.
    ax.legend(handles=handles, title=group_col, fontsize=8, title_fontsize=8,
              frameon=False, loc="upper left", bbox_to_anchor=(0, -0.16),
              ncol=min(len(groups), 4))

    ax = axes[1]
    upg = pd.Series(independence["units_per_group"])
    ax.bar([str(i) for i in upg.index], upg.to_numpy(),
           color=[cmap.get(g, "#999999") for g in upg.index], edgecolor="white")
    for i, v in enumerate(upg.to_numpy()):
        ax.text(i, v, f"n={v}", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel(f"independent units (count)")
    ax.set_xlabel(group_col)
    ax.set_title("n available to a test", fontsize=10, loc="left")
    ax.set_ylim(0, float(max(upg.to_numpy())) * 1.3)

    fig.suptitle("F1 - Design map: what the sample actually contains", fontsize=12,
                 x=0.01, ha="left")
    fig.text(0.01, -0.05, independence["note"], fontsize=8, color="#444444", ha="left",
             wrap=True)
    return _finish(fig, path)


def fig_raw_behind_claim(groups: Dict[str, Sequence[float]], value_label: str,
                         path: Optional[str] = None, test_name: Optional[str] = None,
                         p_value: Optional[float] = None,
                         reported_error: Optional[str] = None):
    """
    F2 - The raw data behind a headline claim.

    Left: every observation with a box summary, so identical means over
    different distributions cannot hide.
    Right: the mean with SD and the mean with SEM drawn side by side. Putting
    both error definitions on one axis makes the sqrt(n) difference visible,
    which is the clearest way to show what a bar-and-SEM plot was hiding.

    `reported_error` is the paper's own stated error definition, printed on the
    panel so reader and reviewer compare like with like.
    """
    names = list(groups)
    arrays = [np.asarray(groups[k], float) for k in names]
    arrays = [a[np.isfinite(a)] for a in arrays]
    colors = _palette(len(names))

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4), sharey=True,
                             gridspec_kw={"width_ratios": [1.4, 1]})
    rng = np.random.default_rng(0)

    ax = axes[0]
    bp = ax.boxplot(arrays, positions=range(len(names)), widths=0.5, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.25); patch.set_edgecolor(c)
    for i, (a, c) in enumerate(zip(arrays, colors)):
        ax.scatter(i + rng.uniform(-0.12, 0.12, a.size), a, s=18, color=c,
                   edgecolor="white", linewidth=0.4, zorder=3)
    ylo = ax.get_ylim()[0]
    for i, a in enumerate(arrays):
        ax.text(i, ylo, f"n={a.size}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names)
    ax.set_ylabel(value_label)
    # The test result goes on a second title line rather than a separate text
    # call at y=1.02, which lands on top of the title.
    title = "Every observation, with the box summary"
    if p_value is not None:
        title += f"\np = {p_value:.3g}" + (f" ({test_name})" if test_name else "")
    ax.set_title(title, fontsize=10, loc="left")

    ax = axes[1]
    for i, (a, c) in enumerate(zip(arrays, colors)):
        m, sd = float(np.mean(a)), float(np.std(a, ddof=1))
        sem = sd / np.sqrt(a.size)
        ax.errorbar(i - 0.12, m, yerr=sd, fmt="o", color=c, capsize=5, markersize=6,
                    label="mean +/- SD" if i == 0 else None)
        ax.errorbar(i + 0.12, m, yerr=sem, fmt="s", color=c, capsize=5, markersize=5,
                    markerfacecolor="white", label="mean +/- SEM" if i == 0 else None)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names)
    ax.set_title("SD describes spread; SEM describes the mean", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False, loc="best")
    note = (f"paper states its error bars are: {reported_error}" if reported_error else
            "paper does not state what its error bars represent")
    ax.text(0.5, -0.17, note, transform=ax.transAxes, ha="center", va="top", fontsize=8,
            color="#444444")

    fig.suptitle("F2 - Raw data behind the claim", fontsize=12, x=0.01, ha="left")
    return _finish(fig, path)


def fig_assumption_qq(groups: Dict[str, Sequence[float]], path: Optional[str] = None,
                      paired_differences: Optional[Sequence[float]] = None,
                      value_label: str = "observed value"):
    """
    F3 - Distribution shape, as a Q-Q plot per group.

    For paired designs pass `paired_differences`; the assumption lives on the
    differences, not on either arm.

    Each panel carries n and the shape statistics, and the figure carries the
    interpretation note stating that the diagnostic is descriptive -- never a
    pass/fail stamp.
    """
    from .assumption_checks import distribution_shape, qq_points

    series = dict(groups)
    if paired_differences is not None:
        series = {"paired differences": paired_differences, **series}
    n_panels = len(series)
    fig, axes = plt.subplots(1, n_panels, figsize=(3.5 * n_panels, 3.7), squeeze=False)

    for ax, (name, vals) in zip(axes[0], series.items()):
        q = qq_points(vals)
        shape = distribution_shape(vals, name=name)
        if len(q):
            ax.scatter(q["theoretical"], q["sample"], s=18, color="#0072B2",
                       edgecolor="white", linewidth=0.4)
            if np.isfinite(q.attrs["slope"]):
                xs = np.array([q["theoretical"].min(), q["theoretical"].max()])
                ax.plot(xs, q.attrs["intercept"] + q.attrs["slope"] * xs, color="#D55E00",
                        linewidth=1.2, label="normal reference")
                ax.legend(fontsize=7, frameon=False, loc="upper left")
        ax.set_xlabel("theoretical normal quantile (z)")
        ax.set_ylabel(value_label if ax is axes[0][0] else "")
        ax.set_title(name, fontsize=10, loc="left")
        box = f"n={shape['n']}\nskew={shape['skewness']:.2f}"
        if np.isfinite(shape["excess_kurtosis"]):
            box += f"\nex-kurt={shape['excess_kurtosis']:.2f}"
        if np.isfinite(shape["shapiro_p"]):
            box += f"\nShapiro p={shape['shapiro_p']:.3f}"
        ax.text(0.98, 0.02, box, transform=ax.transAxes, ha="right", va="bottom",
                fontsize=7.5, family="monospace",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="#cccccc"))

    first = distribution_shape(list(series.values())[0], name="first")
    fig.suptitle("F3 - Distribution shape (descriptive, not a pass/fail gate)",
                 fontsize=12, x=0.01, ha="left")
    fig.text(0.01, -0.07, first["interpretation_note"], fontsize=8, color="#444444",
             ha="left", wrap=True)
    return _finish(fig, path)


def fig_claim_ledger(ledger: pd.DataFrame, path: Optional[str] = None,
                     effect_label: str = "effect (paper's units)"):
    """
    F4 - Reported against recomputed, one row per claim.

    Expects columns claim_id, claim_short, verdict, and where available
    effect_reported, effect_recomputed, ci_low, ci_high, p_reported,
    p_recomputed.

    Left: a forest plot of the recomputed effect and interval with the reported
    point overlaid as an open marker, so a discrepancy reads as a gap between
    two markers on one row.
    Right: reported against recomputed p-values on log axes with the 0.05 lines
    drawn, so decision-changing discrepancies land in the off-diagonal
    quadrants.
    """
    led = ledger.copy()
    led["verdict"] = pd.Categorical(led["verdict"].astype(str), categories=VERDICT_ORDER,
                                    ordered=True)
    led = led.sort_values(["verdict", "claim_id"]).reset_index(drop=True)
    ypos = np.arange(len(led))[::-1]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, max(3.4, 0.55 * len(led) + 2.0)),
                             gridspec_kw={"width_ratios": [1.6, 1]})

    ax = axes[0]
    any_effect = False
    for y, (_, row) in zip(ypos, led.iterrows()):
        color = VERDICT_COLORS.get(str(row["verdict"]), "#999999")
        lo = float(row.get("ci_low", np.nan)); hi = float(row.get("ci_high", np.nan))
        rec = float(row.get("effect_recomputed", np.nan))
        rep = float(row.get("effect_reported", np.nan))
        if np.isfinite(lo) and np.isfinite(hi):
            ax.plot([lo, hi], [y, y], color=color, linewidth=2, solid_capstyle="round")
            any_effect = True
        if np.isfinite(rec):
            ax.plot(rec, y, "o", color=color, markersize=7, zorder=3); any_effect = True
        if np.isfinite(rep):
            ax.plot(rep, y, "o", markerfacecolor="white", markeredgecolor=color,
                    markersize=9, markeredgewidth=1.4, zorder=2); any_effect = True
    ax.axvline(0, color="#444444", linestyle=":", linewidth=1)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{r.claim_id}  {str(r.claim_short)[:46]}" for r in led.itertuples()],
                       fontsize=8)
    ax.set_xlabel(effect_label)
    ax.set_title("Recomputed effect (filled) vs reported (open), with 95% CI",
                 fontsize=10, loc="left")
    if not any_effect:
        ax.text(0.5, 0.5, "no claim carried a numeric effect size", ha="center", va="center",
                fontsize=9, color="#555555", transform=ax.transAxes)

    ax = axes[1]
    pr = pd.to_numeric(led.get("p_reported", pd.Series(dtype=float)), errors="coerce")
    pc = pd.to_numeric(led.get("p_recomputed", pd.Series(dtype=float)), errors="coerce")
    mask = np.isfinite(pr) & np.isfinite(pc)
    if bool(mask.any()):
        # Range the axes to the p-values actually present, with a decade of
        # padding: a fixed 1e-12 floor leaves most of the panel empty.
        smallest = float(min(pr[mask].min(), pc[mask].min()))
        floor = max(10 ** np.floor(np.log10(smallest) - 1), 1e-12) if smallest > 0 else 1e-12
        for (_, row), ok in zip(led.iterrows(), mask):
            if not ok:
                continue
            color = VERDICT_COLORS.get(str(row["verdict"]), "#999999")
            xv = max(float(row["p_reported"]), floor)
            yv = max(float(row["p_recomputed"]), floor)
            ax.scatter(xv, yv, s=46, color=color, edgecolor="white", linewidth=0.5, zorder=3)
            ax.annotate(str(row["claim_id"]), (xv, yv), textcoords="offset points",
                        xytext=(5, 3), fontsize=7)
        lims = [floor, 1]
        ax.plot(lims, lims, color="#888888", linewidth=1, linestyle="--",
                label="exact agreement")
        ax.axvline(0.05, color="#D55E00", linewidth=1, linestyle=":")
        ax.axhline(0.05, color="#D55E00", linewidth=1, linestyle=":")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.legend(fontsize=7, frameon=False, loc="lower right")
        # The explanation goes in the figure footer, not inside the axes: claim
        # labels can land anywhere in the panel and an in-axes note clips them.
        fig.text(0.01, -0.03,
                 "Right panel: points off the diagonal that cross 0.05 are decision-changing.",
                 fontsize=8, color="#444444", ha="left")
    else:
        ax.text(0.5, 0.5, "no claim had both a reported and a\nrecomputed p-value",
                ha="center", va="center", fontsize=9, color="#555555", transform=ax.transAxes)
    ax.set_xlabel("p reported in the paper")
    ax.set_ylabel("p recomputed from the data")
    ax.set_title("Reported vs recomputed p-values", fontsize=10, loc="left")

    present = [v for v in VERDICT_ORDER if (led["verdict"].astype(str) == v).any()]
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=VERDICT_COLORS[v], label=v)
               for v in present]
    fig.legend(handles=handles, loc="lower center", ncol=min(4, max(len(handles), 1)),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("F4 - Claim ledger", fontsize=12, x=0.01, ha="left")
    return _finish(fig, path)


def fig_null_distribution(perm_result: dict, boot: Optional[dict] = None,
                          path: Optional[str] = None,
                          statistic_label: Optional[str] = None):
    """
    F5 - Where the observed statistic falls in its own null distribution.

    Left: the permutation null with the observed value marked, annotated with
    the p-value, the resolution limit and the design's combinatorial floor. When
    the floor binds, the panel says so -- a p-value that could not have been
    smaller is a property of the design, not of the effect.
    Right: the bootstrap distribution of the estimate with its interval, which
    answers the separate question of how precise the estimate is.
    """
    stat_label = statistic_label or str(perm_result.get("statistic_name", "test statistic"))
    n_panels = 2 if boot else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(5.6 * n_panels, 3.9), squeeze=False)

    ax = axes[0][0]
    null = np.asarray(perm_result["null"], float)
    ax.hist(null, bins=int(min(60, max(15, np.sqrt(null.size)))), color="#CCCCCC",
            edgecolor="white", linewidth=0.3)
    ax.axvline(perm_result["observed"], color="#D55E00", linewidth=2,
               label=f"observed = {perm_result['observed']:.3g}")
    ax.set_xlabel(f"{stat_label} under the null (label permutation)")
    ax.set_ylabel(f"permutations (count of {perm_result['n_perm']:,})")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    floor_text = (f"\nfloor:      p >= {perm_result['min_achievable_p']:.4f}"
                  if perm_result.get("floor_binds") else "")
    ax.text(0.98, 0.97,
            f"p =         {perm_result['p_value']:.4g}\n"
            f"resolution: {perm_result['p_resolution']:.1e}{floor_text}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, family="monospace",
            bbox=dict(boxstyle="round",
                      facecolor="#FFF1E0" if perm_result.get("floor_binds") else "white",
                      alpha=0.9, edgecolor="#cccccc"))
    ax.set_title("Permutation null", fontsize=10, loc="left")

    if boot:
        ax = axes[0][1]
        dist = np.asarray(boot["distribution"], float)
        ax.hist(dist, bins=int(min(60, max(15, np.sqrt(dist.size)))), color="#CCE3F0",
                edgecolor="white", linewidth=0.3)
        ax.axvspan(boot["ci_low"], boot["ci_high"], color="#0072B2", alpha=0.12,
                   label=f"{int(boot['level'] * 100)}% {boot['method']} CI")
        ax.axvline(boot["point_estimate"], color="#0072B2", linewidth=2,
                   label=f"estimate = {boot['point_estimate']:.3g}")
        ax.set_xlabel(f"bootstrap replicates of {boot['statistic']}")
        ax.set_ylabel(f"resamples (count of {boot['n_boot']:,})")
        ax.legend(fontsize=8, frameon=False, loc="upper left")
        ax.set_title("Bootstrap distribution of the estimate", fontsize=10, loc="left")

    fig.suptitle("F5 - Null distribution and estimate precision", fontsize=12, x=0.01,
                 ha="left")
    fig.text(0.01, -0.06, perm_result.get("note", ""), fontsize=8, color="#444444",
             ha="left", wrap=True)
    return _finish(fig, path)


def fig_multiplicity(pvals: Sequence[float], alpha: float = 0.05,
                     path: Optional[str] = None, labels: Optional[Sequence[str]] = None):
    """
    F6 - The pattern across every test the paper ran.

    Left: p-value histogram against the uniform expectation, with the
    just-below-threshold window highlighted.
    Right: how many results survive no correction, BH and Bonferroni, with the
    number expected by chance alone drawn as a reference line.
    """
    from .multiplicity import discovery_recount, just_below_rate, p_value_histogram

    p = np.asarray(pvals, float); p = p[np.isfinite(p)]
    if p.size == 0:
        return not_assessable_figure("F6 - Multiplicity", "no p-values available", path)

    hist = p_value_histogram(p, bins=20)
    rc = discovery_recount(p, alpha=alpha, labels=labels)
    jb = just_below_rate(p)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4), gridspec_kw={"width_ratios": [1.5, 1]})

    ax = axes[0]
    width = float(hist["bin_high"].iloc[0] - hist["bin_low"].iloc[0])
    bar_colors = ["#CC79A7" if (lo >= 0.04 and hi <= 0.05) else "#999999"
                  for lo, hi in zip(hist["bin_low"], hist["bin_high"])]
    ax.bar(hist["bin_low"], hist["count"], width=width, align="edge",
           color=bar_colors, edgecolor="white", linewidth=0.4)
    ax.axhline(float(hist["expected_if_uniform"].iloc[0]), color="#0072B2", linestyle="--",
               linewidth=1.2, label="expected if all nulls true (uniform)")
    ax.axvline(alpha, color="#D55E00", linestyle=":", linewidth=1.2, label=f"alpha = {alpha}")
    ax.set_xlabel("p-value")
    ax.set_ylabel(f"tests (count of {p.size})")
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.set_title("Distribution of every reported p-value", fontsize=10, loc="left")
    ax.text(0.02, 0.97,
            (f"[0.045, 0.05): {jb['n_in_suspect_window']}\n"
             f"[0.040, 0.045): {jb['n_in_reference_window']}\n"
             f"binomial p = {jb['binomial_p_excess']:.3f}")
            if jb["interpretable"] else
            (f"{p.size} p-values: too few for the\njust-below-threshold diagnostic"),
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5, family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="#cccccc"))

    ax = axes[1]
    summary = rc["summary"]
    bars = ax.bar(summary["correction"], summary["n_significant"],
                  color=["#999999", "#56B4E9", "#0072B2"], edgecolor="white")
    for b, v in zip(bars, summary["n_significant"]):
        ax.text(b.get_x() + b.get_width() / 2, v, str(int(v)), ha="center", va="bottom",
                fontsize=9)
    ax.axhline(rc["expected_false_positives_uncorrected"], color="#D55E00", linestyle="--",
               linewidth=1.2,
               label=f"expected by chance = {rc['expected_false_positives_uncorrected']:.1f}")
    ax.set_ylabel(f"results called significant (of {rc['m_tests']} tests)")
    ax.set_xlabel("multiple-testing correction")
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ax.set_title("Hits surviving each correction", fontsize=10, loc="left")
    ax.tick_params(axis="x", labelrotation=12)

    fig.suptitle("F6 - Multiplicity across the whole analysis", fontsize=12, x=0.01, ha="left")
    fig.text(0.01, -0.06, rc["note"], fontsize=8, color="#444444", ha="left", wrap=True)
    return _finish(fig, path)


def fig_claim_panel(claim_id: str, groups: Dict[str, Sequence[float]], value_label: str,
                    perm_result: dict, boot: Optional[dict] = None,
                    routes: Optional[pd.DataFrame] = None,
                    verdict: str = "not-assessable", path: Optional[str] = None):
    """
    Per-claim small multiple: raw data, null distribution, and the three routes.

    One of these is produced for each headline claim and appended to the report
    as the per-claim supplement. The three-routes panel is the point of the
    figure: it shows whether the claim survives the choice of method or depends
    on it.
    """
    n_panels = 2 + (1 if routes is not None and len(routes) else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(4.3 * n_panels, 4.0), squeeze=False)
    rng = np.random.default_rng(0)
    names = list(groups)
    arrays = [np.asarray(groups[k], float) for k in names]
    arrays = [a[np.isfinite(a)] for a in arrays]
    colors = _palette(len(names))

    ax = axes[0][0]
    for i, (a, c) in enumerate(zip(arrays, colors)):
        ax.scatter(i + rng.uniform(-0.1, 0.1, a.size), a, s=20, color=c,
                   edgecolor="white", linewidth=0.4)
        m, sd = float(np.mean(a)), float(np.std(a, ddof=1))
        ax.errorbar(i, m, yerr=sd, fmt="_", color="black", capsize=6, markersize=18,
                    linewidth=1.4)
    ylo = ax.get_ylim()[0]
    for i, a in enumerate(arrays):
        ax.text(i, ylo, f"n={a.size}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names)
    ax.set_ylabel(value_label)
    ax.set_title("raw observations (bar = mean +/- SD)", fontsize=9, loc="left")

    ax = axes[0][1]
    null = np.asarray(perm_result["null"], float)
    ax.hist(null, bins=40, color="#CCCCCC", edgecolor="white", linewidth=0.3)
    ax.axvline(perm_result["observed"], color="#D55E00", linewidth=2)
    ax.set_xlabel(f"{perm_result.get('statistic_name', 'statistic')} under the null")
    ax.set_ylabel(f"permutations ({perm_result['n_perm']:,})")
    ax.set_title(f"permutation p = {perm_result['p_value']:.4g}", fontsize=9, loc="left")
    if boot:
        ax.text(0.98, 0.97,
                f"estimate {boot['point_estimate']:.3g}\n"
                f"{int(boot['level'] * 100)}% CI [{boot['ci_low']:.3g}, {boot['ci_high']:.3g}]",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5, family="monospace",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="#cccccc"))

    if n_panels == 3:
        ax = axes[0][2]
        r = routes.reset_index(drop=True)
        ypos = np.arange(len(r))[::-1]
        ax.scatter(r["p_value"].clip(lower=1e-6), ypos,
                   color=["#009E73" if s else "#D55E00" for s in r["agrees_at_0.05"]],
                   s=55, zorder=3)
        for y, row in zip(ypos, r.itertuples()):
            if row.min_achievable_p > 0.05:
                ax.plot([row.min_achievable_p, 1], [y, y], color="#CC79A7", linewidth=6,
                        alpha=0.25, solid_capstyle="butt")
        ax.axvline(0.05, color="#D55E00", linestyle=":", linewidth=1.2)
        ax.set_xscale("log"); ax.set_xlim(1e-6, 1)
        ax.set_yticks(ypos); ax.set_yticklabels(r["test"], fontsize=8)
        ax.set_xlabel("p-value (log scale)")
        ax.set_title("three routes to the same claim", fontsize=9, loc="left")
        ax.text(0.02, 0.02, "shaded = p unreachable\nfor this design",
                transform=ax.transAxes, fontsize=7, color="#8A5577", va="bottom")

    color = VERDICT_COLORS.get(verdict, "#999999")
    fig.suptitle(f"{claim_id} - verdict: {verdict}", fontsize=11, x=0.01, ha="left",
                 color=color)
    return _finish(fig, path)
