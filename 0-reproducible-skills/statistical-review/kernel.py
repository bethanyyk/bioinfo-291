"""
Statistical review: high-level workflow.

Load this module to get the whole review pipeline as a handful of calls. The
per-check functions live in `utilities/`; everything here is orchestration.

    paper  = load_paper("paper.pdf")
    data   = load_data("data/raw/figure2_source.csv")
    audit  = audit_reported_statistics(paper["inventory"])
    claim  = review_claim(data, value_col="volume", group_col="arm",
                          unit_col="mouse", claim_id="C1",
                          claim_short="treated tumours are smaller",
                          p_reported=0.033)
    figs   = build_review_figures(data, "mouse", "arm", "volume",
                                  claims=[claim], audit=audit,
                                  outdir="figures/")

`review_claim` proposes a verdict from the evidence, with its reasons. The
proposal is a starting point for the reviewer, not the finding: a verdict is
only final once a human (or the reviewing agent) has read the methods section
and confirmed that the design supports the claim at all.
"""

import os
import sys
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd

_SKILL_ROOT = os.path.dirname(os.path.abspath(__file__))
if _SKILL_ROOT not in sys.path:
    sys.path.insert(0, _SKILL_ROOT)

from utilities import assumption_checks as _ac
from utilities import claim_extraction as _ce
from utilities import multiplicity as _mult
from utilities import recompute as _rc
from utilities import resampling as _rs
from utilities import review_figures as _rf

# The verdict vocabulary. Ordered from strongest support to strongest problem;
# `propose_verdict` returns the first that applies.
VERDICT_ORDER = _rf.VERDICT_ORDER

# The fixed figure spine. Every review attempts all six; a slot that cannot be
# built is rendered as an explicit "not assessable" panel naming the reason.
FIGURE_SPINE = (
    ("F1_design_map", "What is n, and of what unit?"),
    ("F2_raw_behind_claim", "The raw observations behind the headline claim"),
    ("F3_assumption_shape", "Distribution shape, descriptively"),
    ("F4_claim_ledger", "Reported against recomputed, per claim"),
    ("F5_null_distribution", "Where the observed statistic falls in its null"),
    ("F6_multiplicity", "The pattern across every test in the paper"),
)

DEFAULT_ALPHA = 0.05


def load_paper(pdf_path: str) -> dict:
    """
    Read the paper PDF and build the claim inventory.

    Returns {'path', 'pages', 'inventory', 'text_quality', 'summary'}.

    If `text_quality['usable_text_layer']` is False the PDF has no usable text
    layer (scanned image, or extraction failed). That is a reportable finding --
    "claims could not be extracted automatically" -- not an absence of claims.
    """
    pages = _ce.read_pdf_text(pdf_path)
    inventory = _ce.build_claim_inventory(pages)
    quality = _ce.pdf_text_quality(pages)
    return {
        "path": pdf_path,
        "pages": pages,
        "inventory": inventory,
        "text_quality": quality,
        "summary": _ce.inventory_summary(inventory),
    }


def load_data(source: Union[str, pd.DataFrame], **read_kwargs) -> pd.DataFrame:
    """
    Load the paper's raw data from csv/tsv/xlsx/parquet, or pass a DataFrame through.

    Does no cleaning. Inspect what you loaded before testing anything -- the
    lecture's first instruction is to look at the data.
    """
    if isinstance(source, pd.DataFrame):
        return source
    ext = os.path.splitext(str(source))[1].lower()
    if ext in (".csv", ".txt"):
        return pd.read_csv(source, **read_kwargs)
    if ext in (".tsv", ".tab"):
        return pd.read_csv(source, sep="\t", **read_kwargs)
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(source, **read_kwargs)
    if ext == ".parquet":
        return pd.read_parquet(source, **read_kwargs)
    raise ValueError(f"unsupported data format: {ext}")


def describe_data(data: pd.DataFrame, value_col: str, group_col: str,
                  unit_col: Optional[str] = None) -> dict:
    """
    Look at the data before running numbers: shape per group, spread, and -- if
    `unit_col` is given -- whether rows are independent units or clustered
    within them.
    """
    groups = {str(g): d[value_col].to_numpy() for g, d in data.groupby(group_col)}
    out = {
        "n_rows": int(len(data)),
        "groups": {k: int(np.isfinite(v).sum()) for k, v in groups.items()},
        "shape_table": _ac.shape_table(groups),
        "missing_values": int(data[value_col].isna().sum()),
    }
    if unit_col is not None:
        out["independence"] = _ac.design_independence(data, unit_col, group_col, value_col)
    if len(groups) == 2:
        a, b = list(groups.values())
        out["variance"] = _ac.variance_comparison(a, b)
        out["sample_size_notes"] = _ac.small_sample_notes(a.size, b.size)
    return out


def audit_reported_statistics(inventory: dict, alpha: float = DEFAULT_ALPHA) -> dict:
    """
    Check the paper's own arithmetic, using only the PDF.

    Runs statcheck over every (statistic, p-value) pair, recounts the reported
    p-values under BH and Bonferroni, and measures clustering just below the
    threshold.

    Returns {'statcheck', 'n_checked', 'n_inconsistent', 'n_decision_errors',
    'n_not_checkable', 'multiplicity', 'just_below', 'pvalues'}.
    """
    tests = inventory.get("tests", pd.DataFrame())
    checked = _rc.statcheck_table(tests, alpha=alpha) if len(tests) else tests

    pv = inventory.get("p_values", pd.DataFrame())
    pvals = pv["p_reported"].astype(float).to_numpy() if len(pv) else np.array([])
    pvals = pvals[np.isfinite(pvals)]

    n_checked = int(checked["consistent"].notna().sum()) if len(checked) else 0
    return {
        "statcheck": checked,
        "n_checked": n_checked,
        "n_inconsistent": int((checked["consistent"] == False).sum()) if len(checked) else 0,  # noqa: E712
        "n_decision_errors": int((checked["decision_error"] == True).sum()) if len(checked) else 0,  # noqa: E712
        "n_not_checkable": int(len(checked) - n_checked) if len(checked) else 0,
        "multiplicity": _mult.discovery_recount(pvals, alpha=alpha) if pvals.size else None,
        "just_below": _mult.just_below_rate(pvals) if pvals.size else None,
        "pvalues": pvals,
    }


def propose_verdict(evidence: dict, alpha: float = DEFAULT_ALPHA) -> dict:
    """
    Map assembled evidence to a proposed verdict, with reasons.

    Cascade, first match wins:

      not-assessable          -- nothing in the data addresses the claim
      contradicted            -- the data point the other way
      underpowered-resolution -- the design's p-value floor is above alpha
      assumption-dependent    -- the parametric, rank and permutation routes
                                 disagree at alpha
      not-reproducible        -- the routes are unanimous and all contradict the
                                 reported conclusion
      supported-with-caveat   -- the conclusion stands but something about the
                                 design or the reporting qualifies it
      supported               -- reproduces and survives every route

    Two ordering decisions matter:

    * The resolution floor is checked before any discrepancy. Where a design
      cannot produce a p-value below alpha at all, that fact *explains* the
      discrepancy rather than merely co-occurring with it.
    * Route disagreement is checked before the discrepancy. When the routes
      disagree, at least one supports the paper's own decision, so the claim
      depends on the modelling choice rather than failing to reproduce.

    In both cases the discrepancy is still recorded in `reasons`.

    `evidence` keys consumed: routes (DataFrame), permutation (dict),
    p_reported, p_recomputed, effect_reported, effect_recomputed, independence,
    error_bar, survives_correction.

    This is deterministic bookkeeping over the computed diagnostics. It cannot
    see the methods section, so it never returns a verdict that depends on
    design information the reviewer has to read for themselves.
    """
    reasons: List[str] = []
    routes = evidence.get("routes")
    perm = evidence.get("permutation") or {}
    p_rep = float(evidence.get("p_reported", np.nan) or np.nan)
    p_rec = float(evidence.get("p_recomputed", np.nan) or np.nan)
    eff_rep = float(evidence.get("effect_reported", np.nan) or np.nan)
    eff_rec = float(evidence.get("effect_recomputed", np.nan) or np.nan)

    if routes is None and not perm:
        return {"verdict": "not-assessable",
                "reasons": ["no recomputation was possible from the supplied data"]}

    if np.isfinite(eff_rep) and np.isfinite(eff_rec) and eff_rep * eff_rec < 0:
        reasons.append(
            f"recomputed effect ({eff_rec:.3g}) has the opposite sign to the reported "
            f"effect ({eff_rep:.3g})"
        )
        return {"verdict": "contradicted", "reasons": reasons}

    crosses_alpha = bool(
        np.isfinite(p_rep) and np.isfinite(p_rec) and (p_rep < alpha) != (p_rec < alpha)
    )

    if perm.get("floor_binds"):
        reasons.append(
            f"the design's permutation floor is p >= {perm['min_achievable_p']:.4f}, so no "
            f"result from it can reach alpha = {alpha}"
        )
        if crosses_alpha:
            reasons.append(
                f"the reported p = {p_rep:.3g} is below that floor; recomputing gives "
                f"p = {p_rec:.3g}"
            )
        return {"verdict": "underpowered-resolution", "reasons": reasons}

    # Route disagreement is checked before the reported-vs-recomputed discrepancy.
    # When the routes disagree, at least one of them supports the paper's own
    # decision, so the claim depends on the modelling choice rather than failing
    # to reproduce; calling that 'not-reproducible' would misattribute it.
    # 'not-reproducible' is reserved for the case where the routes are unanimous
    # and all of them contradict the reported conclusion.
    if routes is not None and len(routes) and routes["agrees_at_0.05"].nunique() > 1:
        disagreeing = routes.loc[routes["agrees_at_0.05"] != routes["agrees_at_0.05"].mode()[0],
                                 "test"].tolist()
        reasons.append(
            "the routes disagree at alpha: " + "; ".join(
                f"{r.test} p = {r.p_value:.4g}" for r in routes.itertuples()
            ) + f" (minority: {', '.join(disagreeing)})"
        )
        if crosses_alpha:
            reasons.append(
                f"the reported p = {p_rep:.3g} agrees with some routes but not with the "
                f"permutation reference p = {p_rec:.3g}"
            )
        return {"verdict": "assumption-dependent", "reasons": reasons}

    if crosses_alpha:
        reasons.append(
            f"reported p = {p_rep:.3g} and recomputed p = {p_rec:.3g} fall on opposite sides "
            f"of alpha = {alpha}; every route agrees with the recomputed value"
        )
        return {"verdict": "not-reproducible", "reasons": reasons}

    indep = evidence.get("independence") or {}
    if indep.get("pseudoreplication_risk"):
        reasons.append(indep["note"])
    err = evidence.get("error_bar") or {}
    if err.get("resolved") and err.get("best_match") == "SEM":
        reasons.append(
            "the plotted error bars match SEM, which describes uncertainty in the mean rather "
            f"than the spread of observations (factor {err['sd_to_sem_factor']:.2f})"
        )
    survives = evidence.get("survives_correction")
    if survives is not None and survives == "uncorrected":
        reasons.append("the result does not survive multiple-testing correction")

    if reasons:
        return {"verdict": "supported-with-caveat", "reasons": reasons}

    agree = "all routes agree at alpha" if routes is not None else "permutation agrees"
    return {"verdict": "supported",
            "reasons": [f"{agree}; recomputation reproduces the reported result"]}


def review_claim(data: pd.DataFrame, value_col: str, group_col: str,
                 claim_id: str = "C1", claim_short: str = "",
                 unit_col: Optional[str] = None, paired: bool = False,
                 pair_col: Optional[str] = None,
                 p_reported: Optional[float] = None,
                 effect_reported: Optional[float] = None,
                 reported_error: Optional[float] = None,
                 statistic: str = "mean_diff",
                 n_perm: int = _rs.DEFAULT_N_PERM, n_boot: int = _rs.DEFAULT_N_BOOT,
                 alpha: float = DEFAULT_ALPHA, random_state: int = 0,
                 figure_path: Optional[str] = None) -> dict:
    """
    Full review of one two-group claim against the raw data.

    Runs, in order: group description, assumption diagnostics, the three routes
    (parametric / rank / permutation), a permutation null, a bootstrap interval
    for the difference, and -- where `reported_error` is given -- the error-bar
    identity check. Then proposes a verdict.

    For paired designs set `paired=True` and give `pair_col`, the column
    identifying the subject each pair belongs to; the differences are formed
    per subject and the permutation sign-flips within pair.

    Returns a dict with every diagnostic, the proposed verdict and its reasons,
    plus a one-row `ledger_row` DataFrame ready to concatenate into the claim
    ledger for figure F4.
    """
    levels = list(pd.unique(data[group_col].dropna()))
    if len(levels) != 2:
        raise ValueError(f"{group_col} must have exactly 2 levels for this review, got {levels}")

    if paired:
        if pair_col is None:
            raise ValueError("paired review needs pair_col")
        wide = data.pivot_table(index=pair_col, columns=group_col, values=value_col,
                               aggfunc="mean").dropna()
        a = wide[levels[0]].to_numpy(float)
        b = wide[levels[1]].to_numpy(float)
        groups = {str(levels[0]): a, str(levels[1]): b}
        diffs = a - b
        shape = _ac.paired_difference_check(a, b)
        perm = _rs.permutation_test(a, b, paired=True, n_perm=n_perm, random_state=random_state)
        boot = _rs.bootstrap_ci(diffs, np.mean, n_boot=n_boot, random_state=random_state)
        routes = _rs.three_routes(a, b, paired=True, n_perm=n_perm, random_state=random_state)
        observed_effect = float(np.mean(diffs))
        variance = None
    else:
        a = data.loc[data[group_col] == levels[0], value_col].to_numpy(float)
        b = data.loc[data[group_col] == levels[1], value_col].to_numpy(float)
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        groups = {str(levels[0]): a, str(levels[1]): b}
        shape = _ac.shape_table(groups)
        perm = _rs.permutation_test(a, b, statistic=statistic, n_perm=n_perm,
                                    random_state=random_state)
        boot = _bootstrap_difference(a, b, n_boot=n_boot, level=0.95, random_state=random_state)
        routes = _rs.three_routes(a, b, n_perm=n_perm, random_state=random_state)
        observed_effect = float(np.mean(a) - np.mean(b))
        variance = _ac.variance_comparison(a, b)
        diffs = None

    independence = (_ac.design_independence(data, unit_col, group_col, value_col)
                    if unit_col is not None else None)

    error_bar = None
    if reported_error is not None:
        pooled_sd = float(np.std(np.concatenate([a, b]), ddof=1))
        error_bar = _rc.error_bar_identity(reported_error, pooled_sd, int(a.size + b.size))

    # The permutation route is the reference p-value: it assumes the least.
    p_recomputed = float(perm["p_value"])
    evidence = {
        "routes": routes, "permutation": perm,
        "p_reported": p_reported, "p_recomputed": p_recomputed,
        "effect_reported": effect_reported, "effect_recomputed": observed_effect,
        "independence": independence, "error_bar": error_bar,
    }
    verdict = propose_verdict(evidence, alpha=alpha)

    figure = None
    if figure_path:
        figure = _rf.fig_claim_panel(claim_id, groups, value_col, perm, boot, routes,
                                     verdict=verdict["verdict"], path=figure_path)

    ledger_row = pd.DataFrame([{
        "claim_id": claim_id,
        "claim_short": claim_short or f"{value_col} differs between {levels[0]} and {levels[1]}",
        "verdict": verdict["verdict"],
        "effect_reported": effect_reported if effect_reported is not None else np.nan,
        "effect_recomputed": observed_effect,
        "ci_low": boot.get("ci_low", np.nan),
        "ci_high": boot.get("ci_high", np.nan),
        "p_reported": p_reported if p_reported is not None else np.nan,
        "p_recomputed": p_recomputed,
        "reasons": " | ".join(verdict["reasons"]),
    }])

    return {
        "claim_id": claim_id, "claim_short": claim_short,
        "groups": groups, "levels": [str(x) for x in levels],
        "paired": paired, "differences": diffs,
        "shape": shape, "variance": variance, "independence": independence,
        "error_bar": error_bar,
        "routes": routes, "permutation": perm, "bootstrap": boot,
        "effect_recomputed": observed_effect, "p_recomputed": p_recomputed,
        "verdict": verdict["verdict"], "reasons": verdict["reasons"],
        "ledger_row": ledger_row, "figure": figure, "figure_path": figure_path,
    }


def _bootstrap_difference(a, b, n_boot: int = 10000, level: float = 0.95,
                          random_state: int = 0) -> dict:
    """
    Bootstrap interval for a difference in means between two independent groups.

    Resamples each group separately, which preserves the two-sample structure;
    `resampling.bootstrap_ci` handles the one-sample case.
    """
    rng = np.random.default_rng(random_state)
    a = np.asarray(a, float); b = np.asarray(b, float)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        draws[i] = (rng.choice(a, a.size, replace=True).mean()
                    - rng.choice(b, b.size, replace=True).mean())
    point = float(a.mean() - b.mean())
    lo, hi = np.percentile(draws, [100 * (1 - level) / 2, 100 * (1 + level) / 2])
    return {
        "n": int(a.size + b.size), "statistic": "difference in means",
        "point_estimate": point, "ci_low": float(lo), "ci_high": float(hi),
        "level": level, "method": "percentile", "n_boot": int(n_boot),
        "bootstrap_se": float(np.std(draws, ddof=1)),
        "bootstrap_bias": float(draws.mean() - point),
        "distribution": draws,
        "note": ("tail quantiles are estimated from few resamples; raise n_boot above 2000 "
                 "before quoting interval endpoints." if n_boot < 2000 else ""),
    }


def build_review_figures(data: Optional[pd.DataFrame], unit_col: Optional[str],
                         group_col: Optional[str], value_col: Optional[str],
                         claims: Sequence[dict], audit: Optional[dict] = None,
                         outdir: str = "figures", value_label: Optional[str] = None,
                         reported_error_label: Optional[str] = None) -> Dict[str, str]:
    """
    Render the fixed six-figure spine plus one panel per claim.

    Returns {figure_key: path}. Any slot that cannot be built from the supplied
    inputs is still written, as a panel naming the reason -- the spine is fixed
    so that a missing diagnostic stays visible in the report.
    """
    os.makedirs(outdir, exist_ok=True)
    paths: Dict[str, str] = {}
    label = value_label or (value_col or "measured value")

    def _p(key):
        return os.path.join(outdir, f"{key}.png")

    # F1
    key = "F1_design_map"
    if data is not None and unit_col and group_col and value_col:
        _rf.fig_design_map(data, unit_col, group_col, value_col, path=_p(key))
    else:
        _rf.not_assessable_figure(
            "F1 - Design map",
            "the data do not identify the experimental unit; n cannot be verified", _p(key))
    paths[key] = _p(key)

    # F2 and F5 use the first claim as the headline
    head = claims[0] if claims else None
    key = "F2_raw_behind_claim"
    if head:
        _rf.fig_raw_behind_claim(head["groups"], label, path=_p(key),
                                 test_name="permutation (mean difference)",
                                 p_value=head["p_recomputed"],
                                 reported_error=reported_error_label)
    else:
        _rf.not_assessable_figure("F2 - Raw data behind the claim",
                                  "no claim could be matched to the supplied data", _p(key))
    paths[key] = _p(key)

    key = "F3_assumption_shape"
    if head:
        _rf.fig_assumption_qq(head["groups"], path=_p(key), value_label=label,
                              paired_differences=head.get("differences"))
    else:
        _rf.not_assessable_figure("F3 - Distribution shape",
                                  "no claim could be matched to the supplied data", _p(key))
    paths[key] = _p(key)

    key = "F4_claim_ledger"
    if claims:
        ledger = pd.concat([c["ledger_row"] for c in claims], ignore_index=True)
        _rf.fig_claim_ledger(ledger, path=_p(key), effect_label=f"difference in {label}")
    else:
        _rf.not_assessable_figure("F4 - Claim ledger", "no claims were reviewed", _p(key))
    paths[key] = _p(key)

    key = "F5_null_distribution"
    if head:
        _rf.fig_null_distribution(head["permutation"], head["bootstrap"], path=_p(key))
    else:
        _rf.not_assessable_figure("F5 - Null distribution",
                                  "no claim could be matched to the supplied data", _p(key))
    paths[key] = _p(key)

    key = "F6_multiplicity"
    pvals = audit["pvalues"] if audit is not None else np.array([])
    if len(pvals):
        _rf.fig_multiplicity(pvals, path=_p(key))
    else:
        _rf.not_assessable_figure(
            "F6 - Multiplicity",
            "no p-values could be extracted from the PDF, so the set of tests run cannot "
            "be assessed", _p(key))
    paths[key] = _p(key)

    for claim in claims:
        if claim.get("figure_path"):
            paths[f"panel_{claim['claim_id']}"] = claim["figure_path"]
    return paths


def review_summary_table(claims: Sequence[dict]) -> pd.DataFrame:
    """The claim ledger as one table, ready for the report body."""
    if not claims:
        return pd.DataFrame(columns=["claim_id", "claim_short", "verdict"])
    return pd.concat([c["ledger_row"] for c in claims], ignore_index=True)


def verdict_counts(claims: Sequence[dict]) -> pd.DataFrame:
    """Count claims by verdict, in the canonical verdict order, for the summary."""
    led = review_summary_table(claims)
    counts = led["verdict"].value_counts().reindex(VERDICT_ORDER, fill_value=0)
    return counts.rename_axis("verdict").reset_index(name="n_claims")
