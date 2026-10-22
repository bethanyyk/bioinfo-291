"""
Re-derive the paper's reported numbers.

Three independent checks, cheapest first:

1. `statcheck_table` -- does the reported p-value follow from the reported test
   statistic and degrees of freedom? Needs only the PDF.
2. `grim_check` / `error_bar_identity` -- are the reported summaries even
   attainable for the stated n? Needs only the PDF (plus the raw SD).
3. `recompute_from_data` -- does the test, run on the authors' own raw data,
   return what the paper reports? Needs the data.

An inconsistency here is an arithmetic or reporting fact, not a judgement about
the science. Report it as such: state the reported value, the recomputed value,
and the inputs used, so the authors can point at the discrepancy if the inputs
were wrong.
"""

import numpy as np
import pandas as pd
from scipy import stats

DEFAULT_TAIL = 2


def recompute_p(stat_type: str, stat_value: float, df1: float = np.nan,
                df2: float = np.nan, tail: int = DEFAULT_TAIL) -> float:
    """
    Recompute a p-value from a test statistic and its degrees of freedom.

    stat_type: 't', 'F', 'r', 'chi2' or 'z'. F and chi2 are upper-tail by
    convention, so `tail` is ignored for those. Returns nan when required
    degrees of freedom are missing.
    """
    s = abs(float(stat_value))
    if stat_type == "t":
        if not np.isfinite(df1):
            return np.nan
        p = stats.t.sf(s, df1)
    elif stat_type == "F":
        if not (np.isfinite(df1) and np.isfinite(df2)):
            return np.nan
        return float(stats.f.sf(s, df1, df2))
    elif stat_type == "chi2":
        if not np.isfinite(df1):
            return np.nan
        return float(stats.chi2.sf(s, df1))
    elif stat_type == "r":
        if not np.isfinite(df1) or s >= 1:
            return np.nan
        t = s * np.sqrt(df1) / np.sqrt(max(1 - s ** 2, 1e-300))
        p = stats.t.sf(t, df1)
    elif stat_type == "z":
        p = stats.norm.sf(s)
    else:
        return np.nan
    return float(min(p * tail, 1.0))


def _rounding_window(value: float, text: str = None) -> float:
    """Half-width of the rounding interval implied by how a number was printed."""
    if not np.isfinite(value):
        return np.nan
    s = text if text is not None else repr(float(value))
    if "e" in s or "E" in s:
        return abs(value) * 0.05
    decimals = len(s.split(".")[1]) if "." in s else 0
    return 0.5 * 10 ** (-decimals)


def statcheck_table(tests: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """
    Check reported p-values against their reported statistics (statcheck-style).

    Input: the `tests` table from claim_extraction.build_claim_inventory().

    Adds:
      p_recomputed   -- p implied by stat_value and df
      consistent     -- reported p agrees with the recomputed p within the
                        rounding window implied by the reported precision; for
                        '<' and '>' relations, whether the bound actually holds
      decision_error -- reported and recomputed p fall on opposite sides of
                        alpha; this is the subset that can change a conclusion
      note           -- why a row could not be checked

    Rows that cannot be checked (missing df, unsupported statistic) come back
    with consistent = NA. Count them in the report; never silently drop them.
    """
    cols = ["p_recomputed", "consistent", "decision_error", "note"]
    if not len(tests):
        out = tests.copy()
        for c in cols:
            out[c] = pd.Series(dtype="object")
        return out

    out = tests.copy()
    recomputed, consistent, decision_error, notes = [], [], [], []

    for _, row in out.iterrows():
        tail = int(row.get("tail_hint") or DEFAULT_TAIL)
        p_rec = recompute_p(row["stat_type"], row["stat_value"],
                            row.get("df1", np.nan), row.get("df2", np.nan), tail)
        recomputed.append(p_rec)
        p_rep = float(row.get("p_reported", np.nan))
        rel = (row.get("p_relation") or "=")

        if not np.isfinite(p_rec):
            consistent.append(pd.NA); decision_error.append(pd.NA)
            notes.append("not checkable: missing df or unsupported statistic")
            continue
        if not np.isfinite(p_rep):
            consistent.append(pd.NA); decision_error.append(pd.NA)
            notes.append("not checkable: no p-value reported alongside statistic")
            continue

        window = _rounding_window(p_rep)
        if rel in ("<", "<="):
            ok = p_rec < p_rep + window
        elif rel in (">", ">="):
            ok = p_rec > p_rep - window
        else:
            ok = abs(p_rec - p_rep) <= window

        consistent.append(bool(ok))
        decision_error.append(bool((p_rep < alpha) != (p_rec < alpha)) if not ok else False)
        notes.append("" if ok else "reported p does not follow from reported statistic and df")

    out["p_recomputed"] = recomputed
    out["consistent"] = pd.array(consistent, dtype="boolean")
    out["decision_error"] = pd.array(decision_error, dtype="boolean")
    out["note"] = notes
    return out


def grim_check(mean: float, n: int, decimals: int, n_items: int = 1) -> dict:
    """
    GRIM: is a reported mean attainable from n integer-valued observations?

    For data that are counts, or sums of Likert-type items, the mean must be a
    multiple of 1/(n * n_items). A reported mean that is not is an arithmetic
    impossibility -- a typo, a different n than stated, or an unreported
    exclusion.

    Only valid for genuinely integer-valued measurements, and only informative
    while n is small relative to the reported precision. Returns
    {'applicable': False, 'reason': ...} otherwise.
    """
    if n <= 0 or n_items <= 0:
        return {"applicable": False, "reason": "n and n_items must be positive"}
    if n * n_items > 2 * 10 ** decimals:
        return {"applicable": False,
                "reason": f"n*items={n * n_items} too large for {decimals} decimals to constrain"}
    granularity = 1.0 / (n * n_items)
    nearest = round(mean / granularity) * granularity
    tol = 0.5 * 10 ** (-decimals) + 1e-12
    return {
        "applicable": True,
        "mean_reported": float(mean),
        "nearest_attainable": float(nearest),
        "granularity": float(granularity),
        "consistent": bool(abs(mean - nearest) <= tol),
    }


def error_bar_identity(reported_error: float, sd: float, n: int, level: float = 0.95) -> dict:
    """
    Decide which quantity a reported +/- value actually is.

    Given the raw-data SD and n, compute SD, SEM and the half-width of a t-based
    CI, and report which the paper's number matches. Papers that plot SEM while
    describing "variability" are making a confidence statement about the mean,
    not a description of the spread; the two differ by a factor of sqrt(n),
    which is the most common error-bar error there is.
    """
    if n <= 1 or not np.isfinite(sd) or sd <= 0:
        return {"resolved": False, "reason": "need n > 1 and a positive finite SD"}
    sem = sd / np.sqrt(n)
    half_ci = float(stats.t.ppf(0.5 + level / 2, n - 1) * sem)
    candidates = {"SD": float(sd), "SEM": float(sem), f"CI{int(level * 100)}_half_width": half_ci}
    ratios = {k: float(reported_error / v) for k, v in candidates.items()}
    best = min(ratios, key=lambda k: abs(np.log(max(abs(ratios[k]), 1e-300))))
    return {
        "resolved": True,
        "reported_error": float(reported_error),
        "candidates": candidates,
        "ratio_to_candidates": ratios,
        "best_match": best,
        "sd_to_sem_factor": float(np.sqrt(n)),
        "unambiguous": bool(abs(np.log(max(abs(ratios[best]), 1e-300))) < np.log(1.1)),
    }


def effect_from_summary(m1: float, s1: float, n1: int, m2: float, s2: float, n2: int,
                        level: float = 0.95) -> dict:
    """
    Recover a two-group effect size and interval from published summaries.

    Use when the paper reports mean +/- SD and n but no effect size. Returns
    Welch's t and df, the mean difference with its CI, Cohen's d and Hedges' g
    (the small-sample-corrected standardized difference).
    """
    se = float(np.sqrt(s1 ** 2 / n1 + s2 ** 2 / n2))
    diff = float(m1 - m2)
    df = se ** 4 / ((s1 ** 2 / n1) ** 2 / (n1 - 1) + (s2 ** 2 / n2) ** 2 / (n2 - 1))
    tcrit = float(stats.t.ppf(0.5 + level / 2, df))
    s_pool = float(np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2)))
    d = diff / s_pool
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    return {
        "mean_difference": diff,
        "se_difference": se,
        "welch_t": diff / se,
        "welch_df": float(df),
        "p_two_sided": float(2 * stats.t.sf(abs(diff / se), df)),
        "ci_low": diff - tcrit * se,
        "ci_high": diff + tcrit * se,
        "cohens_d": d,
        "hedges_g": d * J,
        "level": level,
    }


def recompute_from_data(x, y=None, paired: bool = False, test: str = "welch_t",
                        alternative: str = "two-sided") -> dict:
    """
    Re-run a named test on the authors' raw data.

    test: 'welch_t', 'student_t', 'paired_t', 'one_sample_t', 'mann_whitney',
    'wilcoxon_signed_rank'.

    This answers "is the number reproducible from the deposited data?". It does
    not answer "was this the right test", which the assumption checks and the
    three-routes comparison address.
    """
    x = np.asarray(x, dtype=float)
    y_arr = None if y is None else np.asarray(y, dtype=float)

    if paired or test in ("paired_t", "wilcoxon_signed_rank"):
        if y_arr is None:
            raise ValueError(f"{test} needs paired y")
        mask = np.isfinite(x) & np.isfinite(y_arr)
        x, y_arr = x[mask], y_arr[mask]
    else:
        x = x[np.isfinite(x)]
        if y_arr is not None:
            y_arr = y_arr[np.isfinite(y_arr)]

    if test == "welch_t":
        res = stats.ttest_ind(x, y_arr, equal_var=False, alternative=alternative)
    elif test == "student_t":
        res = stats.ttest_ind(x, y_arr, equal_var=True, alternative=alternative)
    elif test == "paired_t":
        res = stats.ttest_rel(x, y_arr, alternative=alternative)
    elif test == "one_sample_t":
        res = stats.ttest_1samp(x, 0.0, alternative=alternative)
    elif test == "mann_whitney":
        res = stats.mannwhitneyu(x, y_arr, alternative=alternative)
    elif test == "wilcoxon_signed_rank":
        res = stats.wilcoxon(x, y_arr, alternative=alternative)
    else:
        raise ValueError(f"unknown test: {test}")

    return {
        "test": test,
        "statistic": float(res.statistic),
        "p_value": float(res.pvalue),
        "n1": int(x.size),
        "n2": int(y_arr.size) if y_arr is not None else None,
        "alternative": alternative,
        "paired": bool(paired or test in ("paired_t", "wilcoxon_signed_rank")),
    }


def reproduction_verdict(p_reported: float, p_recomputed: float, alpha: float = 0.05,
                         rel_tol: float = 0.15) -> str:
    """
    Classify how a recomputed p-value relates to the reported one.

    Returns 'reproduced', 'minor-discrepancy' (same side of alpha, differs by
    more than rel_tol), 'decision-changing' (opposite sides of alpha), or
    'not-checkable'.

    Only 'decision-changing' warrants a headline finding; the others belong in
    the claim ledger with both numbers shown.
    """
    if not (np.isfinite(p_reported) and np.isfinite(p_recomputed)):
        return "not-checkable"
    if (p_reported < alpha) != (p_recomputed < alpha):
        return "decision-changing"
    denom = max(abs(p_reported), 1e-300)
    if abs(p_recomputed - p_reported) / denom > rel_tol:
        return "minor-discrepancy"
    return "reproduced"
