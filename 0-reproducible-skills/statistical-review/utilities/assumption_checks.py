"""
Check what a test assumed, against what the data look like.

Design note, which the lecture is emphatic about and which this module enforces
by construction: a normality check is *not* a pass/fail gate. These functions
return descriptions -- shape statistics, Q-Q coordinates, variance ratios -- and
never a boolean "assumptions met". A Shapiro-Wilk p-value is reported where it
is computable, labelled as what it is: a test whose power depends on n, so it
rejects trivial departures at large n and misses real ones at small n.

The assumptions no diagnostic can check are flagged separately and loudly:
independence, representativeness of the sample, and absence of confounding.
Those are design questions, answerable only from the methods section.
"""

from typing import Dict, Sequence

import numpy as np
import pandas as pd
from scipy import stats

# Below this n, the t-distribution's heavier tails matter and a normality
# diagnostic has almost no power. From the lecture: close to normal at N >= 30.
SMALL_N_THRESHOLD = 30

# Shapiro-Wilk is only defined over this range of n.
SHAPIRO_MIN_N, SHAPIRO_MAX_N = 3, 5000

# Variance ratio beyond this makes Student's equal-variance t inadvisable;
# Welch costs almost nothing and is the safer default either way.
VARIANCE_RATIO_FLAG = 2.0

# Assumptions that no diagnostic in this module can evaluate.
UNCHECKABLE_ASSUMPTIONS = (
    "independence of observations",
    "the sample is representative of the stated population",
    "absence of unmeasured confounding",
    "measurements were not selected after seeing the outcome",
)


def qq_points(x: Sequence[float], dist: str = "norm") -> pd.DataFrame:
    """
    Theoretical vs observed quantiles, for a Q-Q plot.

    Returns a DataFrame with columns theoretical and sample; the fitted
    reference line is attached as .attrs['slope'] and .attrs['intercept'].
    Systematic curvature is the signal, not individual point deviation.
    """
    x = np.asarray(x, dtype=float)
    x = np.sort(x[np.isfinite(x)])
    n = x.size
    if n < 2:
        out = pd.DataFrame(columns=["theoretical", "sample"])
        out.attrs.update({"slope": np.nan, "intercept": np.nan, "n": int(n)})
        return out
    # Blom plotting positions
    a = 3.0 / 8 if n <= 10 else 0.5
    probs = (np.arange(1, n + 1) - a) / (n + 1 - 2 * a)
    theo = getattr(stats, dist).ppf(probs)
    out = pd.DataFrame({"theoretical": theo, "sample": x})
    if n > 2:
        slope, intercept = np.polyfit(theo, x, 1)
    else:
        slope, intercept = np.nan, np.nan
    out.attrs.update({"slope": float(slope), "intercept": float(intercept), "n": int(n)})
    return out


def distribution_shape(x: Sequence[float], name: str = "group") -> dict:
    """
    Describe the shape of one sample. No verdict.

    Returns n, mean, median, sd, IQR, skewness, excess kurtosis, the
    Shapiro-Wilk p-value where defined, and `interpretation_note` spelling out
    what that p-value can and cannot tell you at this n.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    out = {
        "name": name, "n": int(n),
        "mean": float(np.mean(x)) if n else np.nan,
        "median": float(np.median(x)) if n else np.nan,
        "sd": float(np.std(x, ddof=1)) if n > 1 else np.nan,
        "iqr": float(np.subtract(*np.percentile(x, [75, 25]))) if n else np.nan,
        "skewness": float(stats.skew(x, bias=False)) if n > 2 else np.nan,
        "excess_kurtosis": float(stats.kurtosis(x, bias=False)) if n > 3 else np.nan,
        "shapiro_p": np.nan,
        "n_below_small_threshold": bool(n < SMALL_N_THRESHOLD),
    }
    if SHAPIRO_MIN_N <= n <= SHAPIRO_MAX_N:
        out["shapiro_p"] = float(stats.shapiro(x).pvalue)
    if n < SMALL_N_THRESHOLD:
        out["interpretation_note"] = (
            f"n={n} < {SMALL_N_THRESHOLD}: a normality test here has little power, so a "
            "non-significant result is not evidence of normality. Read the Q-Q plot for shape, "
            "and prefer a method that does not lean on the shape assumption."
        )
    else:
        out["interpretation_note"] = (
            f"n={n}: averages are close to normal by the central limit theorem regardless of the "
            "shape of individual observations, so a significant normality test need not "
            "invalidate a test about means. Judge the Q-Q plot for tails and outliers instead."
        )
    return out


def shape_table(groups: Dict[str, Sequence[float]]) -> pd.DataFrame:
    """distribution_shape() for several groups, as one table."""
    return pd.DataFrame([distribution_shape(v, name=k) for k, v in groups.items()])


def variance_comparison(x: Sequence[float], y: Sequence[float]) -> dict:
    """
    Compare spread between two groups and say whether equal-variance t is safe.

    Levene (median-centred) is the default because it is robust to
    non-normality; Bartlett is reported alongside but is sensitive to it.
    The recommendation is Welch whenever the ratio exceeds VARIANCE_RATIO_FLAG
    or either group is small.
    """
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    y = np.asarray(y, float); y = y[np.isfinite(y)]
    v1, v2 = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    ratio = max(v1, v2) / min(v1, v2) if min(v1, v2) > 0 else np.inf
    out = {
        "n1": int(x.size), "n2": int(y.size),
        "sd1": float(np.sqrt(v1)), "sd2": float(np.sqrt(v2)),
        "variance_ratio": float(ratio),
        "levene_p": float(stats.levene(x, y, center="median").pvalue) if min(x.size, y.size) > 1 else np.nan,
        "bartlett_p": float(stats.bartlett(x, y).pvalue) if min(x.size, y.size) > 1 else np.nan,
    }
    flag = ratio > VARIANCE_RATIO_FLAG or min(x.size, y.size) < SMALL_N_THRESHOLD
    out["recommended_t_variant"] = "welch_t" if flag else "either (Welch still safe)"
    out["note"] = (
        f"variance ratio {ratio:.2f}: Student's equal-variance t-test is not appropriate here; "
        "Welch's t-test costs negligible power when the variances are in fact equal."
        if ratio > VARIANCE_RATIO_FLAG else
        f"variance ratio {ratio:.2f}: the equal-variance assumption is not contradicted."
    )
    return out


def paired_difference_check(x: Sequence[float], y: Sequence[float]) -> dict:
    """
    Diagnostics for paired data, run on the differences -- where the assumption
    actually lives.

    For a paired t-test the differences carry the normality assumption. For
    Wilcoxon signed-rank the differences must additionally be *symmetric*: a
    zero median is not sufficient. Both are reported, with the symmetry check
    that distinguishes signed-rank from the assumption-light sign test, and a
    flag for the case where direction and magnitude disagree (most subjects
    improve while the mean worsens).
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    d = x[mask] - y[mask]
    shape = distribution_shape(d, name="paired_differences")
    n_pos, n_neg, n_zero = int((d > 0).sum()), int((d < 0).sum()), int((d == 0).sum())
    skew = shape["skewness"]
    asym = bool(np.isfinite(skew) and abs(skew) > 1.0)
    med, mean = (float(np.median(d)), float(np.mean(d))) if d.size else (np.nan, np.nan)
    disagree = bool(d.size and np.isfinite(med) and np.sign(med) != np.sign(mean) and mean != 0)
    return {
        **shape,
        "n_pairs": int(d.size),
        "n_positive": n_pos, "n_negative": n_neg, "n_zero": n_zero,
        "median_difference": med,
        "mean_difference": mean,
        "symmetry_concern": asym,
        "signed_rank_note": (
            "differences are visibly asymmetric (|skew| > 1): Wilcoxon signed-rank assumes "
            "symmetry of the differences, so prefer the sign test or a permutation test on a "
            "statistic you can state."
            if asym else
            "differences show no gross asymmetry; the signed-rank symmetry assumption is not "
            "contradicted."
        ),
        "direction_vs_magnitude_note": (
            "direction and magnitude disagree: most pairs move one way while the mean moves the "
            "other. State which summary answers the scientific question before choosing a test."
            if disagree else ""
        ),
    }


def design_independence(data: pd.DataFrame, unit_col: str, group_col: str,
                        value_col: str) -> dict:
    """
    Detect pseudoreplication: is the n used by the test the n that was randomized?

    `unit_col` is the independent experimental unit (mouse, patient, dish,
    batch); rows are individual observations (cells, wells, technical repeats).

    Returns observations per unit, units per group, a one-way ICC, the design
    effect 1 + (k-1)*ICC and the implied effective n. A high ICC with many
    observations per unit means the effective sample size is far below the row
    count, and a test run on rows reports a p-value that is too small.
    """
    df = data[[unit_col, group_col, value_col]].dropna()
    per_unit = df.groupby(unit_col).size()
    units_per_group = df.groupby(group_col)[unit_col].nunique()
    k = float(per_unit.mean())

    grand = float(df[value_col].mean())
    unit_means = df.groupby(unit_col)[value_col].mean()
    unit_sizes = df.groupby(unit_col).size()
    ss_between = float((unit_sizes * (unit_means - grand) ** 2).sum())
    ss_within = float(df.groupby(unit_col)[value_col].transform(lambda s: s - s.mean()).pow(2).sum())
    n_units, n_obs = int(per_unit.size), int(len(df))
    ms_between = ss_between / (n_units - 1) if n_units > 1 else np.nan
    ms_within = ss_within / (n_obs - n_units) if n_obs > n_units else np.nan

    icc = np.nan
    if np.isfinite(ms_between) and np.isfinite(ms_within) and (ms_between + (k - 1) * ms_within) > 0:
        icc = float(np.clip((ms_between - ms_within) / (ms_between + (k - 1) * ms_within), 0.0, 1.0))

    design_effect = 1 + (k - 1) * icc if np.isfinite(icc) else np.nan
    effective_n = n_obs / design_effect if np.isfinite(design_effect) and design_effect > 0 else np.nan
    risk = bool(k > 1.5 and np.isfinite(icc) and icc > 0.1)
    return {
        "n_observations": n_obs,
        "n_units": n_units,
        "obs_per_unit_min": int(per_unit.min()), "obs_per_unit_max": int(per_unit.max()),
        "obs_per_unit_mean": k,
        "units_per_group": units_per_group.to_dict(),
        "icc": icc,
        "design_effect": float(design_effect) if np.isfinite(design_effect) else np.nan,
        "effective_n": float(effective_n) if np.isfinite(effective_n) else np.nan,
        "pseudoreplication_risk": risk,
        "note": (
            f"{n_obs} observations from {n_units} units ({k:.1f} per unit, ICC={icc:.2f}). "
            f"Effective n is about {effective_n:.1f}, not {n_obs}. A test treating rows as "
            "independent overstates precision."
            if risk else
            f"{n_obs} observations from {n_units} units ({k:.1f} per unit); clustering does not "
            "materially inflate precision at this ICC."
        ),
    }


def small_sample_notes(*ns: int) -> dict:
    """
    What the sample sizes alone imply, before any test is run.

    Includes the coarsest p-value each design can produce -- its resolution
    floor. A design that cannot produce a p-value below 0.05 cannot support a
    claim of significance at 0.05, whatever test is applied.
    """
    from .resampling import min_achievable_p, min_sign_test_p

    ns = tuple(int(n) for n in ns if n)
    notes = []
    if any(n < SMALL_N_THRESHOLD for n in ns):
        notes.append(
            f"n={ns}: below {SMALL_N_THRESHOLD}, the t-distribution's heavier tails matter and "
            "normality diagnostics have little power."
        )
    floors = {}
    if len(ns) == 2:
        floors["permutation_two_sample_min_p"] = min_achievable_p(ns[0], ns[1], two_sided=True)
    if len(ns) == 1:
        floors["sign_test_min_p"] = min_sign_test_p(ns[0])
        floors["paired_permutation_min_p"] = min_achievable_p(ns[0], paired=True, two_sided=True)
    for key, value in floors.items():
        if np.isfinite(value) and value > 0.05:
            notes.append(
                f"{key} = {value:.4f}: this design cannot produce p < 0.05 by that route, so a "
                "significant result from it rests on the parametric model rather than on "
                "something the data can show on their own."
            )
    return {"n": ns, "resolution_floors": floors, "notes": notes}


def uncheckable_assumptions_note() -> str:
    """The assumptions this module cannot test. Include verbatim in the report."""
    bullets = "\n".join(f"  - {a}" for a in UNCHECKABLE_ASSUMPTIONS)
    return (
        "Diagnostics in this report describe distribution shape and spread only. They cannot "
        "evaluate:\n" + bullets +
        "\nThese are design questions. Assess them from the methods section, and say so "
        "explicitly when the methods do not contain enough information to assess them."
    )
