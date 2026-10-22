"""
Permutation, sign and bootstrap procedures, plus the resolution limits that
decide whether a p-value could ever have been small.

Two ideas from the lecture are implemented here as first-class outputs rather
than footnotes:

* **Resolution.** A permutation p-value cannot be finer than (k+1)/(N+1), and a
  design with few possible relabellings has a hard floor -- 3 vs 3 gives a
  minimum two-sided p of 0.10 however large the difference. Every function here
  returns its floor next to its p-value.
* **How you permute defines your null.** Shuffling across the whole dataset
  tests a different hypothesis from shuffling within patient or within cell
  type. `strata` makes that choice explicit and records it in the result.
"""

import math
from typing import Callable, Optional, Sequence, Union

import numpy as np
import pandas as pd
from scipy import stats

# 10k resolves p down to ~1e-4; raise it when the claim under review is smaller.
DEFAULT_N_PERM = 10000
DEFAULT_N_BOOT = 10000


def p_value_resolution(n_perm: int) -> float:
    """Smallest p-value obtainable from n_perm draws: 1/(n_perm+1)."""
    return 1.0 / (n_perm + 1)


def min_achievable_p(n1: int, n2: Optional[int] = None, paired: bool = False,
                     two_sided: bool = True) -> float:
    """
    The combinatorial floor on a permutation p-value for a given design.

    Two-sample: there are C(n1+n2, n1) distinct relabellings, so the most
    extreme arrangement gives a one-sided p of 1/C.
    Paired (sign-flip): there are 2**n_pairs sign assignments.

    A floor above 0.05 means the design cannot demonstrate significance by
    permutation, whatever the effect size. Report the floor, not just the p.

    Computed in log space. The number of relabellings overflows a float long
    before the sample sizes get interesting -- C(1154, 577) has 345 digits, and
    a direct 1.0/total raises OverflowError -- so the exponent is accumulated
    with lgamma and only then exponentiated, underflowing to a true 0.0 rather
    than raising. A floor of 0.0 is the honest answer at those sizes: the design
    imposes no practical limit, and the binding constraint is the number of
    permutations actually run (`p_value_resolution`).
    """
    if paired:
        log_total = int(n1) * math.log(2.0)
    else:
        if n2 is None:
            raise ValueError("two-sample floor needs n2")
        log_total = _log_comb(int(n1) + int(n2), int(n1))
    return _exp_p(-log_total + (math.log(2.0) if two_sided else 0.0))


def _log_comb(n: int, k: int) -> float:
    """log C(n, k) via lgamma, so huge binomials never touch a float."""
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def _exp_p(log_p: float) -> float:
    """exp() a log-probability, clamping to [0, 1] and underflowing to 0.0."""
    if log_p < -745.0:          # below the smallest positive double
        return 0.0
    return min(1.0, math.exp(log_p))


def min_sign_test_p(n_pairs: int, two_sided: bool = True) -> float:
    """
    Smallest p the sign test can return for n_pairs non-zero differences.

    With 5 pairs the smallest two-sided p is 0.0625 -- above 0.05, so a sign
    test on 5 pairs cannot reach significance however consistent the direction.
    Log space for the same reason as `min_achievable_p`.
    """
    if n_pairs <= 0:
        return np.nan
    return _exp_p(-int(n_pairs) * math.log(2.0) + (math.log(2.0) if two_sided else 0.0))


def _statistic_fn(statistic: Union[str, Callable]) -> Callable:
    if callable(statistic):
        return statistic
    fns = {
        "mean_diff": lambda a, b: float(np.mean(a) - np.mean(b)),
        "median_diff": lambda a, b: float(np.median(a) - np.median(b)),
        "abs_mean_diff": lambda a, b: float(abs(np.mean(a) - np.mean(b))),
        "var_ratio": lambda a, b: float(np.var(a, ddof=1) / np.var(b, ddof=1)),
        "t_welch": lambda a, b: float(stats.ttest_ind(a, b, equal_var=False).statistic),
        "rank_sum": lambda a, b: float(stats.rankdata(np.concatenate([a, b]))[:len(a)].sum()),
        "n_above_pooled_median": lambda a, b: float(
            (np.asarray(a) > np.median(np.concatenate([a, b]))).sum()
        ),
    }
    if statistic not in fns:
        raise ValueError(f"unknown statistic {statistic!r}; pass a callable or one of {sorted(fns)}")
    return fns[statistic]


def permutation_test(x: Sequence[float], y: Optional[Sequence[float]] = None,
                     statistic: Union[str, Callable] = "mean_diff",
                     n_perm: int = DEFAULT_N_PERM, paired: bool = False,
                     alternative: str = "two-sided", strata: Optional[Sequence] = None,
                     random_state: int = 0) -> dict:
    """
    Permutation test, with its own resolution reported alongside the p-value.

    paired=True sign-flips the within-pair differences (the null is "the label
    within each pair is arbitrary"). Otherwise group labels are shuffled; pass
    `strata` (one label per observation, x first then y) to shuffle *within*
    strata, which is how you encode a null such as "treatment is exchangeable
    within cell type".

    `statistic` may be a name ('mean_diff', 'median_diff', 't_welch',
    'rank_sum', 'n_above_pooled_median', 'var_ratio') or any callable f(a, b).
    Different statistics answer different questions and can disagree on the same
    data -- that freedom is the method's cost, so the chosen statistic is
    returned in the result and belongs in the report.

    Returns observed, p_value, null (array), n_perm, p_resolution,
    min_achievable_p and floor_binds (True when the design's floor exceeds 0.05).
    """
    rng = np.random.default_rng(random_state)
    x = np.asarray(x, float)

    if paired:
        if y is None:
            raise ValueError("paired permutation needs y")
        y = np.asarray(y, float)
        if x.size != y.size:
            raise ValueError("paired permutation needs equal-length x and y")
        mask = np.isfinite(x) & np.isfinite(y)
        d = x[mask] - y[mask]
        fn = statistic if callable(statistic) else (lambda arr: float(np.mean(arr)))
        observed = float(fn(d))
        signs = rng.choice(np.array([-1.0, 1.0]), size=(n_perm, d.size))
        null = np.array([float(fn(d * s)) for s in signs])
        floor = min_achievable_p(d.size, paired=True, two_sided=(alternative == "two-sided"))
        stat_name = getattr(statistic, "__name__", "mean_of_differences") if callable(statistic) \
            else "mean_of_differences"
        n1, n2 = int(d.size), None
    else:
        if y is None:
            raise ValueError("two-sample permutation needs y")
        x = x[np.isfinite(x)]
        y = np.asarray(y, float); y = y[np.isfinite(y)]
        fn = _statistic_fn(statistic)
        observed = float(fn(x, y))
        pooled = np.concatenate([x, y])
        n1 = int(x.size)
        null = np.empty(n_perm)
        if strata is None:
            idx = np.arange(pooled.size)
            for i in range(n_perm):
                perm = rng.permutation(idx)
                null[i] = fn(pooled[perm[:n1]], pooled[perm[n1:]])
        else:
            strata = np.asarray(strata)
            if strata.size != pooled.size:
                raise ValueError("strata must have one label per observation (x then y)")
            labels = np.concatenate([np.zeros(n1, bool), np.ones(y.size, bool)])
            uniq = np.unique(strata)
            for i in range(n_perm):
                lab = labels.copy()
                for s in uniq:
                    m = strata == s
                    lab[m] = rng.permutation(lab[m])
                null[i] = fn(pooled[~lab], pooled[lab])
        floor = min_achievable_p(x.size, y.size, two_sided=(alternative == "two-sided"))
        stat_name = getattr(statistic, "__name__", "custom") if callable(statistic) else statistic
        n2 = int(y.size)

    centre = float(np.mean(null))
    if alternative == "two-sided":
        k = int(np.sum(np.abs(null - centre) >= abs(observed - centre)))
    elif alternative == "greater":
        k = int(np.sum(null >= observed))
    elif alternative == "less":
        k = int(np.sum(null <= observed))
    else:
        raise ValueError("alternative must be two-sided, greater or less")

    p = (k + 1) / (n_perm + 1)
    return {
        "statistic_name": stat_name,
        "observed": observed,
        "p_value": float(p),
        "n_perm": int(n_perm),
        "n_extreme": k,
        "p_resolution": p_value_resolution(n_perm),
        "min_achievable_p": float(floor),
        "floor_binds": bool(floor > 0.05),
        "paired": bool(paired),
        "stratified": strata is not None,
        "alternative": alternative,
        "n1": n1, "n2": n2,
        "null": null,
        "note": (
            f"design floor is p >= {floor:.4f}: this comparison cannot reach p < 0.05 by "
            "permutation regardless of effect size."
            if floor > 0.05 else
            f"resolution limit {p_value_resolution(n_perm):.2e} at n_perm={n_perm}; raise n_perm "
            "to resolve smaller p-values."
        ),
    }


def sign_test(differences: Sequence[float], alternative: str = "two-sided") -> dict:
    """
    Exact sign test on paired differences.

    Keeps only the direction of each change, so it assumes nothing about the
    shape or symmetry of the differences -- only that each pair is an
    independent trial, Binomial(n, 1/2) under the null. Use it when the
    scientific question is "do most subjects improve", which is a different
    question from "does the mean change".
    """
    d = np.asarray(differences, float)
    d = d[np.isfinite(d)]
    n_zero = int((d == 0).sum())
    d = d[d != 0]
    n, n_pos = int(d.size), int((d > 0).sum())
    if n == 0:
        return {"n_pairs": 0, "n_zero_dropped": n_zero, "p_value": np.nan,
                "min_achievable_p": np.nan, "floor_binds": False,
                "note": "no non-zero differences"}
    res = stats.binomtest(n_pos, n, 0.5, alternative=alternative)
    floor = min_sign_test_p(n, two_sided=(alternative == "two-sided"))
    return {
        "n_pairs": n, "n_zero_dropped": n_zero,
        "n_positive": n_pos, "n_negative": n - n_pos,
        "proportion_positive": n_pos / n,
        "p_value": float(res.pvalue),
        "min_achievable_p": float(floor),
        "floor_binds": bool(floor > 0.05),
        "alternative": alternative,
        "note": (
            f"with {n} non-zero pairs the smallest two-sided p is {floor:.4f}; significance at "
            "0.05 is unreachable for this design."
            if floor > 0.05 else f"{n_pos} of {n} pairs moved positive."
        ),
    }


def bootstrap_ci(x: Sequence[float], statistic: Callable = np.mean,
                 n_boot: int = DEFAULT_N_BOOT, level: float = 0.95,
                 method: str = "BCa", random_state: int = 0) -> dict:
    """
    Bootstrap confidence interval for any statistic of one sample.

    Resamples with replacement to imitate repeating the experiment. The
    assumption is that the sample is representative of the population -- the
    bootstrap propagates a biased sample's bias, it does not repair it.

    method: 'BCa' (bias-corrected and accelerated, the default) or
    'percentile'. Tail quantiles need more resamples than the centre, so n_boot
    below ~2000 gives unstable interval endpoints.
    """
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if x.size < 2:
        return {"n": int(x.size), "ci_low": np.nan, "ci_high": np.nan,
                "note": "need at least 2 finite observations"}
    res = stats.bootstrap(
        (x,), statistic, n_resamples=n_boot, confidence_level=level,
        method=method.lower(), vectorized=False,
        rng=np.random.default_rng(random_state),
    )
    point = float(statistic(x))
    boot = np.asarray(res.bootstrap_distribution, float)
    return {
        "n": int(x.size),
        "statistic": getattr(statistic, "__name__", "statistic"),
        "point_estimate": point,
        "ci_low": float(res.confidence_interval.low),
        "ci_high": float(res.confidence_interval.high),
        "level": level, "method": method, "n_boot": int(n_boot),
        "bootstrap_se": float(np.std(boot, ddof=1)),
        "bootstrap_bias": float(np.mean(boot) - point),
        "distribution": boot,
        "note": (
            "tail quantiles are estimated from few resamples; raise n_boot above 2000 before "
            "quoting interval endpoints." if n_boot < 2000 else ""
        ),
    }


def three_routes(x: Sequence[float], y: Optional[Sequence[float]] = None,
                 paired: bool = False, n_perm: int = DEFAULT_N_PERM,
                 random_state: int = 0) -> pd.DataFrame:
    """
    Run the parametric, rank-based and permutation routes on the same data.

    Returns one row per route with its p-value, what it assumes, and its
    resolution floor. Agreement across routes makes a claim robust to the
    modelling choice; disagreement localises the claim's dependence on an
    assumption, which is exactly what a statistical review should report.

    This is a comparison of routes, not a menu to pick the smallest p from.
    Choosing among them after seeing the results is p-hacking; the report shows
    all of them.
    """
    x = np.asarray(x, float)
    rows = []

    if paired:
        y = np.asarray(y, float)
        mask = np.isfinite(x) & np.isfinite(y)
        xv, yv = x[mask], y[mask]
        d = xv - yv
        rows.append({"route": "parametric", "test": "paired t-test",
                     "p_value": float(stats.ttest_rel(xv, yv).pvalue),
                     "assumes": "differences are normal; pairs independent",
                     "min_achievable_p": 0.0})
        rows.append({"route": "non-parametric (ranks)", "test": "Wilcoxon signed-rank",
                     "p_value": float(stats.wilcoxon(xv, yv).pvalue),
                     "assumes": "differences are symmetric about the null median",
                     "min_achievable_p": 0.0})
        st = sign_test(d)
        rows.append({"route": "non-parametric (signs)", "test": "sign test",
                     "p_value": st["p_value"],
                     "assumes": "pairs independent; nothing about shape",
                     "min_achievable_p": st["min_achievable_p"]})
        pt = permutation_test(xv, yv, paired=True, n_perm=n_perm, random_state=random_state)
        rows.append({"route": "permutation", "test": "sign-flip permutation (mean difference)",
                     "p_value": pt["p_value"],
                     "assumes": "within-pair exchangeability under the null",
                     "min_achievable_p": pt["min_achievable_p"]})
    else:
        xv = x[np.isfinite(x)]
        yv = np.asarray(y, float); yv = yv[np.isfinite(yv)]
        rows.append({"route": "parametric", "test": "Welch's t-test",
                     "p_value": float(stats.ttest_ind(xv, yv, equal_var=False).pvalue),
                     "assumes": "group means normal; unequal variances allowed",
                     "min_achievable_p": 0.0})
        rows.append({"route": "parametric", "test": "Student's t-test",
                     "p_value": float(stats.ttest_ind(xv, yv, equal_var=True).pvalue),
                     "assumes": "normal; equal variances",
                     "min_achievable_p": 0.0})
        rows.append({"route": "non-parametric (ranks)", "test": "Mann-Whitney U",
                     "p_value": float(stats.mannwhitneyu(xv, yv, alternative="two-sided").pvalue),
                     "assumes": "independence; tests which group tends to be larger, "
                                "not equality of medians",
                     "min_achievable_p": 0.0})
        pt = permutation_test(xv, yv, statistic="mean_diff", n_perm=n_perm,
                              random_state=random_state)
        rows.append({"route": "permutation", "test": "label permutation (mean difference)",
                     "p_value": pt["p_value"],
                     "assumes": "exchangeability: identical distributions under the null",
                     "min_achievable_p": pt["min_achievable_p"]})

    out = pd.DataFrame(rows)
    out["reachable_at_0.05"] = out["min_achievable_p"] <= 0.05
    out["agrees_at_0.05"] = out["p_value"] < 0.05
    out.attrs["unanimous"] = bool(out["agrees_at_0.05"].nunique() == 1)
    return out
