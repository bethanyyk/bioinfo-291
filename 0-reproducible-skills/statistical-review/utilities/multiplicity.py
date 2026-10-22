"""
Multiple testing and p-hacking diagnostics.

Two questions, given the set of p-values a paper reports (or a set recomputed
from its data):

1. **How many of these survive correction?** Bonferroni controls the chance of
   *any* false positive and does not look at the data. Benjamini-Hochberg
   controls the expected *fraction* of false positives among the calls and
   adapts to how many hits there are. They answer different questions;
   `discovery_recount` shows both, so the report can state what changes.

2. **Does the distribution of p-values look steered?** Under a true null,
   p-values are uniform; under a real effect they pile up near zero. Neither
   produces a spike immediately below 0.05. `just_below_rate` measures that
   spike. It is a property of a *set* of p-values and is never evidence about
   any single result.
"""

from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats

SUSPECT_WINDOW = (0.045, 0.05)
REFERENCE_WINDOW = (0.04, 0.045)


def bonferroni(pvals: Sequence[float], alpha: float = 0.05) -> pd.DataFrame:
    """Bonferroni-adjusted p-values. Data-independent, controls FWER."""
    p = np.asarray(pvals, float)
    adj = np.minimum(p * p.size, 1.0)
    return pd.DataFrame({"p_raw": p, "p_adjusted": adj, "reject": adj < alpha})


def benjamini_hochberg(pvals: Sequence[float], alpha: float = 0.05) -> pd.DataFrame:
    """
    Benjamini-Hochberg adjusted p-values (FDR).

    Data-dependent: the effective threshold rises with the number of small
    p-values, which is why BH recovers hits Bonferroni discards. The price is
    that up to alpha of the calls are expected to be false.
    """
    p = np.asarray(pvals, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adj_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(adj_sorted, 1.0)
    return pd.DataFrame({"p_raw": p, "p_adjusted": adj, "reject": adj < alpha})


def discovery_recount(pvals: Sequence[float], alpha: float = 0.05,
                      labels: Optional[Sequence[str]] = None) -> dict:
    """
    Recount discoveries under no correction, Bonferroni and BH.

    Returns {'table': per-test DataFrame, 'summary': counts DataFrame,
    'm_tests', 'alpha', 'expected_false_positives_uncorrected', 'note'}.

    The expected-false-positive figure is the lecture's point made arithmetic:
    at alpha = 0.05 over m tests you expect 0.05*m significant results even if
    every null is true. Put that number next to the paper's hit count.
    """
    p = np.asarray(pvals, float)
    p = p[np.isfinite(p)]
    m = p.size
    if m == 0:
        return {"table": pd.DataFrame(), "summary": pd.DataFrame(), "m_tests": 0,
                "alpha": alpha, "expected_false_positives_uncorrected": 0.0,
                "note": "no finite p-values"}

    bonf = bonferroni(p, alpha)
    bh = benjamini_hochberg(p, alpha)
    table = pd.DataFrame({
        "label": (list(labels)[:m] if labels is not None else [f"test_{i + 1}" for i in range(m)]),
        "p_raw": p,
        "significant_uncorrected": p < alpha,
        "p_bonferroni": bonf["p_adjusted"].to_numpy(),
        "significant_bonferroni": bonf["reject"].to_numpy(),
        "p_bh": bh["p_adjusted"].to_numpy(),
        "significant_bh": bh["reject"].to_numpy(),
    })
    n_raw = int(table["significant_uncorrected"].sum())
    n_bonf = int(table["significant_bonferroni"].sum())
    n_bh = int(table["significant_bh"].sum())
    expected_fp = alpha * m

    summary = pd.DataFrame([
        {"correction": "none", "controls": "nothing (per-test error rate)",
         "n_significant": n_raw, "threshold_on_raw_p": alpha},
        {"correction": "Bonferroni", "controls": "FWER: probability of any false positive",
         "n_significant": n_bonf, "threshold_on_raw_p": alpha / m},
        {"correction": "Benjamini-Hochberg",
         "controls": f"FDR: expected fraction of false calls <= {alpha}",
         "n_significant": n_bh,
         "threshold_on_raw_p": float(np.max(p[bh["reject"].to_numpy()])) if n_bh else 0.0},
    ])
    return {
        "table": table, "summary": summary, "m_tests": m, "alpha": alpha,
        "expected_false_positives_uncorrected": float(expected_fp),
        "note": (
            f"{m} tests at alpha={alpha} yield about {expected_fp:.1f} significant results by "
            f"chance alone. {n_raw} uncorrected hits become {n_bh} under BH and {n_bonf} under "
            "Bonferroni."
        ),
    }


def max_statistic_null(observed: Sequence[float], null_matrix: np.ndarray,
                       alternative: str = "greater") -> pd.DataFrame:
    """
    Permutation FWER control by the max-statistic method.

    `null_matrix` is (n_permutations, n_tests): the statistic for every test in
    every permutation. Each observed statistic is compared against the
    distribution of the *best* statistic across all tests within each
    permutation, which handles multiplicity without assuming independence
    between tests -- the correlation structure is carried by the permutation
    itself.

    This is the right correction for a genome-wide scan in which neighbouring
    features are correlated, where Bonferroni is both conservative and wrong
    about the dependence.
    """
    observed = np.asarray(observed, float)
    null_matrix = np.asarray(null_matrix, float)
    if null_matrix.ndim != 2:
        raise ValueError("null_matrix must be (n_permutations, n_tests)")
    if alternative == "greater":
        max_null = null_matrix.max(axis=1)
        k = np.array([(max_null >= o).sum() for o in observed])
    elif alternative == "two-sided":
        max_null = np.abs(null_matrix).max(axis=1)
        k = np.array([(max_null >= abs(o)).sum() for o in observed])
    else:
        raise ValueError("alternative must be 'greater' or 'two-sided'")
    n_perm = null_matrix.shape[0]
    return pd.DataFrame({
        "observed": observed,
        "p_fwer": (k + 1) / (n_perm + 1),
        "n_perm": n_perm,
        "p_resolution": 1.0 / (n_perm + 1),
    })


def p_value_histogram(pvals: Sequence[float], bins: int = 20) -> pd.DataFrame:
    """
    Binned p-value distribution with the uniform expectation alongside.

    Read the shape: a spike at the left means real signal; a flat profile means
    mostly nulls; a spike anywhere else -- especially just under 0.05 -- is
    about the analysis, not the biology.
    """
    p = np.asarray(pvals, float)
    p = p[np.isfinite(p)]
    counts, edges = np.histogram(p, bins=bins, range=(0, 1))
    return pd.DataFrame({
        "bin_low": edges[:-1], "bin_high": edges[1:],
        "count": counts,
        "expected_if_uniform": p.size / bins,
    })


def just_below_rate(pvals: Sequence[float], window: tuple = SUSPECT_WINDOW,
                    reference: tuple = REFERENCE_WINDOW) -> dict:
    """
    Measure clustering of p-values immediately below a threshold.

    Compares the count in `window` (default 0.045-0.05) with the count in the
    adjacent, equally wide `reference` window (0.04-0.045). Under any continuous
    distribution of p-values -- null or not -- the two windows should hold
    similar counts, and the left one slightly more where there is real signal.
    An excess just under the threshold is the measurable signature of analytic
    choices made until something crossed.

    Interpretation limits the report must carry: this is a property of a
    *collection* of p-values and says nothing about any individual result; it
    needs a reasonable number of p-values (20+) to mean anything; and rounding
    conventions alone can produce apparent pile-up at 0.05 exactly.
    """
    p = np.asarray(pvals, float)
    p = p[np.isfinite(p)]
    n_win = int(((p >= window[0]) & (p < window[1])).sum())
    n_ref = int(((p >= reference[0]) & (p < reference[1])).sum())
    total = n_win + n_ref
    binom_p = float(stats.binomtest(n_win, total, 0.5, alternative="greater").pvalue) \
        if total else np.nan
    return {
        "n_in_suspect_window": n_win, "suspect_window": window,
        "n_in_reference_window": n_ref, "reference_window": reference,
        "ratio": float(n_win / n_ref) if n_ref else (np.inf if n_win else np.nan),
        "binomial_p_excess": binom_p,
        "n_pvalues_total": int(p.size),
        "exactly_at_threshold": int((p == 0.05).sum()),
        "interpretable": bool(p.size >= 20),
        "note": (
            "fewer than 20 p-values: this diagnostic is not interpretable at this scale; report "
            "the counts without drawing a conclusion."
            if p.size < 20 else
            f"{n_win} p-values in {window} vs {n_ref} in the equally wide {reference} "
            f"(one-sided binomial p = {binom_p:.3f})."
        ),
    }


def reported_vs_corrected_claims(claim_table: pd.DataFrame, p_col: str = "p_reported",
                                 alpha: float = 0.05) -> pd.DataFrame:
    """
    Apply `discovery_recount` to a paper's own extracted p-value table.

    Adds adjusted p-values and a `survives` column naming which corrections each
    claim survives, so the claim ledger can show exactly which headline results
    depend on no correction having been applied.
    """
    rows = claim_table[np.isfinite(claim_table[p_col].astype(float))].copy()
    if not len(rows):
        return rows
    rc = discovery_recount(rows[p_col].astype(float).to_numpy(), alpha=alpha)
    t = rc["table"]
    rows["p_bonferroni"] = t["p_bonferroni"].to_numpy()
    rows["p_bh"] = t["p_bh"].to_numpy()
    rows["survives"] = [
        ", ".join([n for n, flag in (("uncorrected", u), ("BH", b), ("Bonferroni", f)) if flag])
        or "none"
        for u, b, f in zip(t["significant_uncorrected"], t["significant_bh"],
                           t["significant_bonferroni"])
    ]
    rows.attrs["summary"] = rc["summary"]
    rows.attrs["note"] = rc["note"]
    return rows
