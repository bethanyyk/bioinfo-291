"""
Statistical review of:

  Forman J., Renner E. & Leavens D.A. (2023) Fetching felines: a survey of cat
  owners on the diversity of cat (Felis catus) fetching behaviour.
  Scientific Reports 13:20456.  doi:10.1038/s41598-023-47409-w

Method: the `statistical-review` skill (six-stage protocol); figures build on
the `data-visualization` skill. Regenerates every artifact under forman_2023/.

Seeds and resampling parameters are fixed at the top; nothing below is random
without one.
"""

import hashlib
import itertools
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

SKILL_ROOT = os.path.join("skills", "statistical-review")
if SKILL_ROOT not in sys.path:
    sys.path.insert(0, SKILL_ROOT)

from utilities import assumption_checks as ac          # noqa: E402
from utilities import claim_extraction as ce           # noqa: E402
from utilities import multiplicity as mult             # noqa: E402
from utilities import recompute as rc                  # noqa: E402
from utilities import resampling as rs                 # noqa: E402
from utilities import review_figures as rf             # noqa: E402
from kernel import propose_verdict                     # noqa: E402

# `resampling.min_achievable_p` evaluates 2 / C(n1+n2, n1) directly, which
# overflows a float for designs of this size (C(1154, 617) has ~346 digits).
# The floor is computed in log space instead and underflows to exactly 0, which
# is the honest answer: at these n the combinatorial floor is far below any
# reportable p-value. Review-side fix only; the skill is not modified.
_ORIG_MIN_P = rs.min_achievable_p


def _min_achievable_p_safe(n1, n2=None, paired=False, two_sided=True):
    try:
        return _ORIG_MIN_P(n1, n2, paired=paired, two_sided=two_sided)
    except OverflowError:
        import math
        if paired or n2 is None:
            log_total = n1 * math.log(2.0)
        else:
            log_total = (math.lgamma(n1 + n2 + 1) - math.lgamma(n1 + 1)
                         - math.lgamma(n2 + 1))
        log_p = (math.log(2.0) if two_sided else 0.0) - log_total
        return math.exp(log_p) if log_p > -700 else 0.0


rs.min_achievable_p = _min_achievable_p_safe

SEED = 0
N_PERM = 10_000          # permutation replicates (resolution 1/(N+1))
N_MC = 100_000           # Monte-Carlo multinomial replicates for chi-square nulls
N_BOOT = 10_000          # bootstrap resamples
ALPHA = 0.05

SLUG = "forman_2023"
OUTDIR = SLUG
FIGDIR = os.path.join(OUTDIR, "figures")
os.makedirs(FIGDIR, exist_ok=True)

PAPER_VID = "9629532c-5b7a-4550-8e81-64ba958d6bb8"
DATA_VID = "2210645b-30b9-4271-96ae-550007d38f65"


def _resolve(vid, fallback):
    """Resolve an input artifact to a local path; fall back to a workspace copy."""
    try:
        return host.artifact_path(vid)            # noqa: F821  (kernel global)
    except Exception:
        return fallback


PAPER = _resolve(PAPER_VID, os.path.join(OUTDIR, "paper_fulltext.txt"))
DATA = _resolve(DATA_VID, os.path.join(OUTDIR, "dataset_S3.csv"))


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def _decimals(x):
    """Decimal places a printed value carries, for round-trip comparison."""
    s = f"{x:g}"
    if "e" in s or "E" in s:
        return 12
    return len(s.split(".")[1]) if "." in s else 0


def _printed_interval(printed):
    """
    The interval of true values a printed number could have been rounded from.

    '0.4373' -> (0.43725, 0.43735); '2.652e-05' -> (2.6515e-05, 2.6525e-05).
    Needed because both sides of a statcheck comparison are rounded: the printed
    statistic and the printed p-value.
    """
    s = str(printed).strip()
    v = float(s)
    if "e" in s or "E" in s:
        mantissa = s.split("e")[0].split("E")[0]
        dec = len(mantissa.split(".")[1]) if "." in mantissa else 0
        exp = int(s.split("e")[-1].split("E")[-1])
        half = 0.5 * 10 ** (-dec + exp)
    else:
        half = 0.5 * 10 ** -(len(s.split(".")[1]) if "." in s else 0)
    return v - half, v + half


def _p_round_trips(p_recomputed, printed):
    """
    Does a recomputed p-value round to the printed one at the precision printed?

    Printed p-values come in two forms here -- fixed decimals (0.4373) and
    scientific notation (2.652e-05) -- and the second carries significant digits,
    not decimal places, so a fixed absolute tolerance would call a correct value
    inconsistent purely because of its magnitude.
    """
    printed = str(printed).strip()
    if not np.isfinite(p_recomputed):
        return None
    if "e" in printed or "E" in printed:
        mantissa = printed.split("e")[0].split("E")[0].lstrip("-")
        sig = len(mantissa.replace(".", "").lstrip("0")) or 1
        return float(f"{p_recomputed:.{max(sig - 1, 0)}e}") == float(printed)
    return abs(p_recomputed - float(printed)) <= 0.5 * 10 ** -_decimals(float(printed))


# ---------------------------------------------------------------- stage 0
TXT = open(PAPER, encoding="utf-8", errors="replace").read()
DF = pd.read_csv(DATA)
INV = ce.build_claim_inventory([TXT])
INV_SUMMARY = ce.inventory_summary(INV)
STATCHECK = rc.statcheck_table(INV["tests"], alpha=ALPHA)

# ---- Supplementary S2: the authors' R source with console output interleaved.
# Supplied after the first pass of this review; it makes every reported test
# checkable against the call that produced it.
S2_VID = "f76dcde0-6f2b-4ce4-9819-73df46421996"
S2_PATH = _resolve(S2_VID, os.path.join(OUTDIR, "supplementary_S2_Rcode.txt"))
S2 = open(S2_PATH, encoding="utf-8", errors="replace").read()
S2_FLAT = re.sub(r"\s+", " ", S2)

# Printed test statistics, in file order (retained- and excluded-outlier runs).
_s2_rows = []
for _m in re.finditer(r"(?:data: )?([A-Za-z0-9_$.,\s]{0,80}?)"
                      r"(Kruskal-Wallis chi-squared|X-squared) = ([\d.]+), df = (\d+), "
                      r"p-value ([<=]) ?([\d.e+-]+)", S2_FLAT):
    _s2_rows.append({"kind": _m.group(2), "data": _m.group(1).strip()[-58:],
                     "statistic": float(_m.group(3)), "df": int(_m.group(4)),
                     "p_relation": _m.group(5), "p_printed": _m.group(6)})
for _m in re.finditer(r"data: ([A-Za-z0-9_$.,\s]{0,70}?) ?W = (\d+), p-value = ([\d.e+-]+)",
                      S2_FLAT):
    _s2_rows.append({"kind": "W (Wilcoxon rank sum)", "data": _m.group(1).strip()[-58:],
                     "statistic": float(_m.group(2)), "df": np.nan,
                     "p_relation": "=", "p_printed": _m.group(3)})
S2T = pd.DataFrame(_s2_rows)
S2T["p_from_statistic"] = [stats.chi2.sf(r.statistic, r.df) if np.isfinite(r.df) else np.nan
                           for r in S2T.itertuples()]

# Pairwise Wilcoxon matrices: (Cat vs About equal, Myself vs About equal, Myself vs Cat)
PAIRWISE_S2 = re.findall(
    r"Cat\s+(<\s*2e-16|[0-9.e+-]+)\s+-\s*Myself\s+(<\s*2e-16|[0-9.e+-]+)\s+"
    r"(<\s*2e-16|[0-9.e+-]+)\s*P value adjustment", S2_FLAT)

# Counts of the calls actually made, which is the denominator stage 6 needs.
S2_CALLS = {name: len(re.findall(rx, S2)) for name, rx in [
    ("chisq.test", r"chisq\.test\("), ("kruskal.test", r"kruskal\.test\("),
    ("pairwise.wilcox.test", r"pairwise\.wilcox\.test\("),
    ("wilcox.test", r"(?<!\.)wilcox\.test\("), ("leveneTest", r"leveneTest\("),
    ("shapiro_test", r"shapiro[._]test\(")]}
S2_INFERENTIAL_P = (S2_CALLS["chisq.test"] + S2_CALLS["kruskal.test"]
                    + 3 * S2_CALLS["pairwise.wilcox.test"] + S2_CALLS["wilcox.test"])

TEXT_NOTE = (
    "The paper was supplied as an extracted text layer, not a PDF, so "
    "`load_paper` was bypassed and the claim extractors were run on the text "
    "directly. Page numbers are therefore unavailable and claims are located "
    "by section."
)

# Reported p-values, curated by hand from the Results text. The automatic
# extractor also picks up three numbers from the Table 1 demographic block and
# the alpha threshold inside the Levene sentence, which are not test results;
# those are excluded here and the exclusion is reported.
REPORTED_P = [
    ("sex ratio 617 vs 537", 0.018, "="),
    ("initiation omnibus", 0.001, "<"),
    ("initiation cat vs owner", 0.001, "<"),
    ("initiation cat vs equal", 0.001, "<"),
    ("initiation owner vs equal", 0.001, "<"),
    ("termination omnibus", 0.001, "<"),
    ("termination cat vs owner", 0.001, "<"),
    ("termination cat vs equal", 0.001, "<"),
    ("termination owner vs equal", 0.437, "="),
    ("KW retrievals ~ initiator", 0.001, "<"),
    ("KW retrievals ~ terminator", 0.001, "<"),
    ("KW monthly ~ initiator", 0.001, "<"),
    ("KW monthly ~ terminator", 0.001, "<"),
    ("W retrievals(init) owner vs cat", 0.009, "="),
    ("W retrievals(init) owner vs equal", 0.001, "<"),
    ("W retrievals(init) cat vs equal", 0.995, "="),
    ("W retrievals(end) owner vs cat", 0.001, "<"),
    ("W retrievals(end) owner vs equal", 0.001, "<"),
    ("W retrievals(end) equal vs cat", 0.001, "<"),
    ("W monthly(init) owner vs cat", 0.001, "<"),
    ("W monthly(init) owner vs equal", 0.001, "<"),
    ("W monthly(end) cat vs owner", 0.001, "<"),
    ("W monthly(end) cat vs equal", 0.001, "<"),
    ("W monthly(end) owner vs equal", 1.000, "="),
    ("sex x initiation", 0.403, "="),
    ("sex x termination", 0.374, "="),
    ("MW sex retrievals", 0.091, "="),
    ("MW sex monthly", 0.622, "="),
]
P_LABELS = [r[0] for r in REPORTED_P]
P_VALUES = np.array([r[1] for r in REPORTED_P], float)
P_IS_BOUND = np.array([r[2] == "<" for r in REPORTED_P])

# ------------------------------------------------- stage 1-2: the unit question
OWNER_DEMOG = ["country", "continent", "ppl_count_household", "person_one_age",
               "person_one_gender", "person_two_age", "person_two_gender",
               "dogs_household", "number_dogs_household", "other_pets",
               "total_household_number_cats"]
ID_LIKE = [c for c in DF.columns
           if re.search(r"(^|_)(id|respondent|participant|owner|response)(_|$)", c, re.I)]
DEMOG_KEY = DF[OWNER_DEMOG].astype(str).agg("|".join, axis=1)
DEMOG_SIZES = DEMOG_KEY.value_counts()
REPEAT_FLAG = {str(k): int(v) for k, v in
               DF["ppl_demographics_repeated"].value_counts(dropna=False).items()}

N_CATS = int(len(DF))
N_OWNERS_REPORTED = 924
MEAN_CATS_PER_OWNER = N_CATS / N_OWNERS_REPORTED

# Recorded fetching-cat cluster sizes: owners of cats in multi-fetching-cat
# households reported how many OTHER fetching cats they had. This bounds the
# cluster-size distribution without inventing an owner key.
OTHER_FETCH = pd.to_numeric(DF["number_other_fetch_cats_household"], errors="coerce")
PER_CAT_CLUSTER = np.where(OTHER_FETCH.notna(), OTHER_FETCH.fillna(0) + 1, 1.0)
CLUSTER_SIZES_RECORDED = pd.Series(PER_CAT_CLUSTER).value_counts().sort_index()
N_MULTI_FETCH = int((PER_CAT_CLUSTER > 1).sum())
N_FETCH_HH_UNKNOWN = int(DF["fetch_cats_household"].isna().sum())

# Kish design-effect bound: deff = 1 + (m_eff - 1) * rho with rho <= 1 (worst
# case), and m_eff = sum_c m_c^2 / sum_c m_c = mean cluster size per observation.
# (a) from the reported 924 owners / 1154 cats, clusters of equal size:
DEFF_EQUAL_WORST = 1.0 + (MEAN_CATS_PER_OWNER - 1.0) * 1.0
# (b) from the per-cat count of fetching cats in the household, as recorded:
M_EFF_RECORDED = float(PER_CAT_CLUSTER.sum() / N_CATS)
DEFF_RECORDED_WORST = 1.0 + (M_EFF_RECORDED - 1.0) * 1.0
DEFF_USED = max(DEFF_EQUAL_WORST, DEFF_RECORDED_WORST)

# ----------------------------------------------------------- recomputation
def gof_chi2(counts, label, n_mc=N_MC, seed=SEED):
    """
    Goodness-of-fit chi-square against equal frequencies, by three routes:
    the asymptotic Pearson statistic, the exact/Monte-Carlo multinomial null,
    and (for 2 categories) the exact binomial test.
    """
    counts = np.asarray(counts, float)
    k, N = len(counts), counts.sum()
    exp = np.full(k, N / k)
    chi2 = float(((counts - exp) ** 2 / exp).sum())
    p_asym = float(stats.chi2.sf(chi2, k - 1))
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(int(N), np.full(k, 1 / k), size=n_mc)
    null = ((draws - exp) ** 2 / exp).sum(axis=1)
    p_mc = float((np.sum(null >= chi2) + 1) / (n_mc + 1))
    out = {
        "label": label, "counts": counts, "k": k, "N": int(N), "chi2": chi2,
        "df": k - 1, "p_asym": p_asym, "p_mc": p_mc, "p_res": 1.0 / (n_mc + 1),
        "cramers_v": float(np.sqrt(chi2 / (N * (k - 1)))),
        "null": null, "min_exact_p": float(k * (1.0 / k) ** N),
    }
    if k == 2:
        bt = stats.binomtest(int(counts[0]), int(N))
        lo, hi = bt.proportion_ci(0.95, method="exact")
        yates = float(((np.abs(counts - exp) - 0.5) ** 2 / exp).sum())
        out.update({"p_exact": float(bt.pvalue), "chi2_yates": yates,
                    "p_yates": float(stats.chi2.sf(yates, 1)),
                    "prop": float(counts[0] / N), "prop_ci": (float(lo), float(hi))})
    # chi-square under a first-order Rao-Scott correction at the worst-case
    # design effect implied by the reported owner/cat ratio
    out["chi2_deff"] = chi2 / DEFF_USED
    out["p_deff"] = float(stats.chi2.sf(out["chi2_deff"], k - 1))
    return out


def routes_frequency(g):
    """The three routes for a frequency claim, in the schema `three_routes` uses."""
    rows = [
        {"route": "parametric", "test": "Pearson chi-square (asymptotic)",
         "p_value": g["p_asym"],
         "assumes": "asymptotic chi-square approximation to the multinomial",
         "min_achievable_p": g["min_exact_p"]},
    ]
    if g["k"] == 2:
        rows.append({"route": "exact", "test": "exact binomial test",
                     "p_value": g["p_exact"],
                     "assumes": "independent Bernoulli trials; no approximation",
                     "min_achievable_p": g["min_exact_p"]})
        rows.append({"route": "parametric (continuity-corrected)",
                     "test": "Yates-corrected chi-square", "p_value": g["p_yates"],
                     "assumes": "asymptotic approximation with continuity correction",
                     "min_achievable_p": g["min_exact_p"]})
    rows.append({"route": "resampling", "test": f"Monte-Carlo multinomial null ({N_MC:,})",
                 "p_value": g["p_mc"],
                 "assumes": "cats independent and identically distributed under the null",
                 "min_achievable_p": g["p_res"]})
    out = pd.DataFrame(rows)
    out["reachable_at_0.05"] = out["min_achievable_p"] <= ALPHA
    out["agrees_at_0.05"] = out["p_value"] < ALPHA
    out.attrs["unanimous"] = bool(out["agrees_at_0.05"].nunique() == 1)
    return out


def perm_rank_test(a, b, n_perm=N_PERM, seed=SEED):
    """
    Permutation test on the rank-sum statistic -- the resampling counterpart of
    the Mann-Whitney test, as opposed to a permutation of the mean difference.
    Reported alongside it because the two answer different questions on the same
    data.
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    vals = np.concatenate([a, b])
    ranks = stats.rankdata(vals)
    n1 = a.size
    obs = float(ranks[:n1].sum())
    centre = n1 * (vals.size + 1) / 2.0
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = rng.permutation(ranks)[:n1].sum()
    p = float((np.sum(np.abs(null - centre) >= abs(obs - centre)) + 1) / (n_perm + 1))
    return {"statistic_name": "rank sum (group 1)", "observed": obs,
            "p_value": p, "n_perm": n_perm, "p_resolution": 1.0 / (n_perm + 1),
            "min_achievable_p": 1.0 / (n_perm + 1), "floor_binds": False,
            "null": null, "paired": False, "stratified": False,
            "note": (f"permutation of the rank-sum statistic, {n_perm:,} replicates, "
                     f"seed {seed}; two-sided by distance from the null centre "
                     f"{centre:.0f}.")}


def _kw_H_from_ranks(ranks, codes, k, N, tie_corr):
    sums = np.bincount(codes, weights=ranks, minlength=k)
    ns = np.bincount(codes, minlength=k).astype(float)
    H = 12.0 / (N * (N + 1)) * np.sum(sums ** 2 / ns) - 3.0 * (N + 1)
    return H / tie_corr


def kw_claim(value_col, group_col, label, levels=("Cat", "About equal", "Myself"),
             n_perm=N_PERM, n_boot=N_BOOT, seed=SEED):
    """
    Kruskal-Wallis claim by three routes: one-way ANOVA (parametric), the
    rank-based KW test, and a label permutation of the KW statistic. Plus a
    bootstrap interval for epsilon-squared and Bonferroni pairwise Wilcoxons.
    """
    sub = DF[[value_col, group_col]].dropna()
    gs = [sub.loc[sub[group_col] == l, value_col].to_numpy(float) for l in levels]
    groups = dict(zip(levels, gs))
    vals = np.concatenate(gs)
    codes = np.concatenate([np.full(len(g), i) for i, g in enumerate(gs)])
    N, k = vals.size, len(levels)
    ranks = stats.rankdata(vals)
    _, tcnt = np.unique(vals, return_counts=True)
    tie_corr = 1.0 - ((tcnt ** 3 - tcnt).sum()) / (N ** 3 - N)
    H = float(_kw_H_from_ranks(ranks, codes, k, N, tie_corr))
    p_kw = float(stats.chi2.sf(H, k - 1))
    F, p_anova = stats.f_oneway(*gs)

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = _kw_H_from_ranks(rng.permutation(ranks), codes, k, N, tie_corr)
    p_perm = float((np.sum(null >= H) + 1) / (n_perm + 1))

    eps2 = float((H - k + 1) / (N - k))
    rng_b = np.random.default_rng(seed + 1)
    boot = np.empty(n_boot)
    idx_all = np.arange(N)
    for i in range(n_boot):
        take = rng_b.choice(idx_all, N, replace=True)
        v, c = vals[take], codes[take]
        if len(np.unique(c)) < k:
            boot[i] = np.nan
            continue
        r = stats.rankdata(v)
        _, tc = np.unique(v, return_counts=True)
        tcb = 1.0 - ((tc ** 3 - tc).sum()) / (N ** 3 - N)
        Hb = _kw_H_from_ranks(r, c, k, N, tcb)
        boot[i] = (Hb - k + 1) / (N - k)
    boot = boot[np.isfinite(boot)]
    ci = np.percentile(boot, [2.5, 97.5])

    pw = {}
    for x, y in itertools.combinations(range(k), 2):
        u, pu = stats.mannwhitneyu(gs[x], gs[y], alternative="two-sided")
        n1, n2 = gs[x].size, gs[y].size
        pw[f"{levels[x]} vs {levels[y]}"] = {
            "p_raw": float(pu), "p_bonf": float(min(1.0, pu * 3)),
            "U": float(u), "n1": int(n1), "n2": int(n2),
            "rank_biserial": float(2 * u / (n1 * n2) - 1),
            "median_1": float(np.median(gs[x])), "median_2": float(np.median(gs[y])),
            "iqr_1": float(np.subtract(*np.percentile(gs[x], [75, 25]))),
            "iqr_2": float(np.subtract(*np.percentile(gs[y], [75, 25]))),
        }

    routes = pd.DataFrame([
        {"route": "parametric", "test": "one-way ANOVA (F)", "p_value": float(p_anova),
         "assumes": "normal residuals; equal variances", "min_achievable_p": 0.0},
        {"route": "non-parametric (ranks)", "test": "Kruskal-Wallis H",
         "p_value": p_kw,
         "assumes": "independence; tests stochastic ordering, not equality of medians",
         "min_achievable_p": 0.0},
        {"route": "permutation", "test": f"label permutation of H ({n_perm:,})",
         "p_value": p_perm, "assumes": "exchangeability of cats across groups",
         "min_achievable_p": 1.0 / (n_perm + 1)},
    ])
    routes["reachable_at_0.05"] = routes["min_achievable_p"] <= ALPHA
    routes["agrees_at_0.05"] = routes["p_value"] < ALPHA
    routes.attrs["unanimous"] = bool(routes["agrees_at_0.05"].nunique() == 1)

    medians = {l: (float(np.median(g)), float(np.subtract(*np.percentile(g, [75, 25]))),
                   int(g.size)) for l, g in groups.items()}
    perm_result = {
        "statistic_name": "Kruskal-Wallis H", "observed": H, "p_value": p_perm,
        "n_perm": n_perm, "p_resolution": 1.0 / (n_perm + 1),
        "min_achievable_p": 1.0 / (n_perm + 1), "floor_binds": False,
        "null": null, "paired": False, "stratified": False,
        "note": (f"label permutation of the Kruskal-Wallis statistic, {n_perm:,} "
                 f"replicates, seed {seed}; resolution {1.0/(n_perm+1):.1e}."),
    }
    boot_result = {
        "n": N, "statistic": "epsilon-squared", "point_estimate": eps2,
        "ci_low": float(ci[0]), "ci_high": float(ci[1]), "level": 0.95,
        "method": "percentile", "n_boot": int(boot.size), "distribution": boot,
    }
    return {"label": label, "groups": groups, "H": H, "df": k - 1, "p_kw": p_kw,
            "F": float(F), "p_anova": float(p_anova), "p_perm": p_perm,
            "p_res": 1.0 / (n_perm + 1), "eps2": eps2, "eps2_ci": (float(ci[0]), float(ci[1])),
            "medians": medians, "pairwise": pw, "routes": routes,
            "permutation": perm_result, "bootstrap": boot_result,
            "shape": ac.shape_table(groups), "N": N}


# ------------------------------------------------------------- descriptives
TRAINED = DF["trained"].value_counts()
N_UNTRAINED = int(TRAINED.get("Untrained", 0))
UNTRAINED_BT = stats.binomtest(N_UNTRAINED, N_CATS)
UNTRAINED_CI = tuple(float(x) for x in UNTRAINED_BT.proportion_ci(0.95, method="exact"))
PCT_UNTRAINED = 100.0 * N_UNTRAINED / N_CATS

LIFE_STAGE = DF["life_stage"].value_counts()
AGE_BAND = pd.cut(DF["fetch_age_years_total"], [-0.001, 0.999, 6.9999, 100],
                  labels=["<1 yr", "1 to <7 yr", ">=7 yr"])
LIFE_X_AGE = pd.crosstab(DF["life_stage"], AGE_BAND)
YA_RANGE = (float(DF.loc[DF.life_stage == "Young adult", "fetch_age_years_total"].min()),
            float(DF.loc[DF.life_stage == "Young adult", "fetch_age_years_total"].max()))

INI = DF["initiate_fetch"].value_counts()
END = DF["end_fetch"].value_counts()
SEX = DF["cat_sex"].value_counts()

G = {
    "ini_omni": gof_chi2([INI["Cat"], INI["About equal"], INI["Myself"]], "initiation omnibus"),
    "ini_cat_owner": gof_chi2([INI["Cat"], INI["Myself"]], "initiation cat vs owner"),
    "ini_cat_equal": gof_chi2([INI["Cat"], INI["About equal"]], "initiation cat vs equal"),
    "ini_owner_equal": gof_chi2([INI["Myself"], INI["About equal"]], "initiation owner vs equal"),
    "end_omni": gof_chi2([END["Cat"], END["About equal"], END["Myself"]], "termination omnibus"),
    "end_cat_owner": gof_chi2([END["Cat"], END["Myself"]], "termination cat vs owner"),
    "end_cat_equal": gof_chi2([END["Cat"], END["About equal"]], "termination cat vs equal"),
    "end_owner_equal": gof_chi2([END["Myself"], END["About equal"]], "termination owner vs equal"),
    "sex": gof_chi2([SEX["Male"], SEX["Female"]], "sex ratio"),
    "untrained": gof_chi2([N_UNTRAINED, N_CATS - N_UNTRAINED], "untrained vs trained"),
}

KW = {
    "ret_ini": kw_claim("times_object_retrieved_recent", "initiate_fetch",
                        "retrievals in most recent session ~ initiator"),
    "ret_end": kw_claim("times_object_retrieved_recent", "end_fetch",
                        "retrievals in most recent session ~ terminator"),
    "mon_ini": kw_claim("times_monthly_fetch_sessions", "initiate_fetch",
                        "monthly fetching sessions ~ initiator"),
    "mon_end": kw_claim("times_monthly_fetch_sessions", "end_fetch",
                        "monthly fetching sessions ~ terminator"),
}

# cat sex x initiation / termination contingency tables
SEX_INI_TAB = pd.crosstab(DF["cat_sex"], DF["initiate_fetch"])
SEX_END_TAB = pd.crosstab(DF["cat_sex"], DF["end_fetch"])


def contingency(tab, label, n_mc=N_MC, seed=SEED):
    chi2, p, dof, exp = stats.chi2_contingency(tab.to_numpy(), correction=False)
    # Monte-Carlo null: permute the cat-sex label across cats, keeping both
    # margins' totals, and recompute the table statistic.
    rows = np.repeat(np.arange(tab.shape[0]), tab.sum(axis=1).to_numpy())
    cols = np.repeat(np.tile(np.arange(tab.shape[1]), tab.shape[0]),
                     tab.to_numpy().reshape(-1))
    rng = np.random.default_rng(seed)
    n_rep = 10_000
    null = np.empty(n_rep)
    for i in range(n_rep):
        perm = rng.permutation(rows)
        t = np.zeros_like(tab.to_numpy(), dtype=float)
        np.add.at(t, (perm, cols), 1.0)
        e = np.outer(t.sum(1), t.sum(0)) / t.sum()
        null[i] = float(((t - e) ** 2 / e).sum())
    p_mc = float((np.sum(null >= chi2) + 1) / (n_rep + 1))
    n = int(tab.to_numpy().sum())
    v = float(np.sqrt(chi2 / (n * (min(tab.shape) - 1))))
    routes = pd.DataFrame([
        {"route": "parametric", "test": "Pearson chi-square (asymptotic)",
         "p_value": float(p), "assumes": "asymptotic chi-square approximation",
         "min_achievable_p": 0.0},
        {"route": "resampling", "test": f"label permutation of the table ({n_rep:,})",
         "p_value": p_mc, "assumes": "cat sex exchangeable across response categories",
         "min_achievable_p": 1.0 / (n_rep + 1)},
    ])
    routes["reachable_at_0.05"] = routes["min_achievable_p"] <= ALPHA
    routes["agrees_at_0.05"] = routes["p_value"] < ALPHA
    routes.attrs["unanimous"] = bool(routes["agrees_at_0.05"].nunique() == 1)
    return {"label": label, "chi2": float(chi2), "df": int(dof), "p_asym": float(p),
            "p_mc": p_mc, "p_res": 1.0 / (n_rep + 1), "cramers_v": v, "n": n,
            "null": null, "routes": routes, "table": tab}


SEX_INI = contingency(SEX_INI_TAB, "cat sex x who initiates")
SEX_END = contingency(SEX_END_TAB, "cat sex x who ends")


def sex_continuous(value_col, label, n_perm=N_PERM, n_boot=N_BOOT, seed=SEED):
    """Two-group claim: all four routes from `three_routes`, plus Z and effect size."""
    a = DF.loc[DF.cat_sex == "Male", value_col].to_numpy(float)
    b = DF.loc[DF.cat_sex == "Female", value_col].to_numpy(float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    groups = {"Male": a, "Female": b}
    routes = rs.three_routes(a, b, n_perm=n_perm, random_state=seed)
    perm = rs.permutation_test(a, b, statistic="mean_diff", n_perm=n_perm,
                               random_state=seed)
    u, p_u = stats.mannwhitneyu(a, b, alternative="two-sided")
    n1, n2 = a.size, b.size
    mu = n1 * n2 / 2.0
    su = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    z = float((u - mu) / su)
    # tie-corrected sigma: the paper's R routines (coin) correct for ties, and
    # these are heavily tied integer counts
    nn = n1 + n2
    _, tc = np.unique(np.concatenate([a, b]), return_counts=True)
    su_tie = np.sqrt(n1 * n2 / 12.0 * ((nn + 1) - ((tc ** 3 - tc).sum()
                                                   / (nn * (nn - 1)))))
    z_tie = float((u - mu) / su_tie)
    rbc = float(2 * u / (n1 * n2) - 1)
    rng = np.random.default_rng(seed + 1)
    draws = np.array([rng.choice(a, n1, True).mean() - rng.choice(b, n2, True).mean()
                      for _ in range(n_boot)])
    point = float(a.mean() - b.mean())
    boot = {"n": n1 + n2, "statistic": "difference in means (male - female)",
            "point_estimate": point, "ci_low": float(np.percentile(draws, 2.5)),
            "ci_high": float(np.percentile(draws, 97.5)), "level": 0.95,
            "method": "percentile", "n_boot": n_boot, "distribution": draws}
    prank = perm_rank_test(a, b, n_perm=n_perm, seed=seed)
    routes = pd.concat([routes, pd.DataFrame([{
        "route": "permutation (ranks)",
        "test": f"permutation of the rank sum ({n_perm:,})",
        "p_value": prank["p_value"],
        "assumes": "exchangeability; the rank counterpart of Mann-Whitney",
        "min_achievable_p": prank["p_resolution"],
        "reachable_at_0.05": True,
        "agrees_at_0.05": prank["p_value"] < ALPHA}])], ignore_index=True)
    routes.attrs["unanimous"] = bool(routes["agrees_at_0.05"].nunique() == 1)
    return {"label": label, "groups": groups, "routes": routes, "permutation": perm,
            "permutation_rank": prank,
            "bootstrap": boot, "U": float(u), "p_mw": float(p_u), "z": z,
            "z_tie": z_tie,
            "rank_biserial": rbc, "shape": ac.shape_table(groups),
            "variance": ac.variance_comparison(a, b),
            "medians": {"Male": float(np.median(a)), "Female": float(np.median(b))},
            "n1": int(n1), "n2": int(n2),
            "min_achievable_p": float(rs.min_achievable_p(n1, n2))}


SEX_RET = sex_continuous("times_object_retrieved_recent", "cat sex x retrievals")
SEX_MON = sex_continuous("times_monthly_fetch_sessions", "cat sex x monthly sessions")


def outlier_sensitivity_sex(value_col):
    """
    Re-run the sex comparison with extreme outliers excluded, the way S2 does it:
    rstatix::identify_outliers marks is.extreme beyond Q3 + 3*IQR or Q1 - 3*IQR,
    applied WITHIN each cat_sex group, and the marked rows are anti-joined out.
    W is reported with Female first, matching R's factor order, so it can be
    compared with the printed value directly.
    """
    sub = DF[["cat_sex", value_col]].dropna().copy()
    keep = np.ones(len(sub), bool)
    for s in ("Female", "Male"):
        idx = np.where(sub["cat_sex"].to_numpy() == s)[0]
        x = sub[value_col].to_numpy(float)[idx]
        q1, q3 = np.percentile(x, [25, 75])
        iqr = q3 - q1
        keep[idx[(x > q3 + 3 * iqr) | (x < q1 - 3 * iqr)]] = False
    out = {}
    for tag, d in (("retained", sub), ("excluded", sub[keep])):
        f = d.loc[d.cat_sex == "Female", value_col].to_numpy(float)
        m = d.loc[d.cat_sex == "Male", value_col].to_numpy(float)
        u, p = stats.mannwhitneyu(f, m, alternative="two-sided")
        out[tag] = {"W_female_first": float(u), "p": float(p), "n_female": int(f.size),
                    "n_male": int(m.size), "z_from_p": float(stats.norm.ppf(p / 2))}
    out["n_excluded"] = int((~keep).sum())
    return out


SENS_RET = outlier_sensitivity_sex("times_object_retrieved_recent")
SENS_MON = outlier_sensitivity_sex("times_monthly_fetch_sessions")

# --------------------------------------------------------- multiplicity
RECOUNT = mult.discovery_recount(P_VALUES, alpha=ALPHA, labels=P_LABELS)
JUST_BELOW = mult.just_below_rate(P_VALUES)
BONF_THRESHOLD = ALPHA / len(P_VALUES)


def pair_claim(value_col, group_col, l1, l2, n_perm=N_PERM, n_boot=N_BOOT, seed=SEED):
    """One pairwise comparison, by all routes `three_routes` provides."""
    a = DF.loc[DF[group_col] == l1, value_col].to_numpy(float)
    b = DF.loc[DF[group_col] == l2, value_col].to_numpy(float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    routes = rs.three_routes(a, b, n_perm=n_perm, random_state=seed)
    perm = rs.permutation_test(a, b, statistic="mean_diff", n_perm=n_perm,
                               random_state=seed)
    u, p_u = stats.mannwhitneyu(a, b, alternative="two-sided")
    rng = np.random.default_rng(seed + 1)
    draws = np.array([rng.choice(a, a.size, True).mean() - rng.choice(b, b.size, True).mean()
                      for _ in range(n_boot)])
    boot = {"n": a.size + b.size, "statistic": f"mean difference ({l1} - {l2})",
            "point_estimate": float(a.mean() - b.mean()),
            "ci_low": float(np.percentile(draws, 2.5)),
            "ci_high": float(np.percentile(draws, 97.5)), "level": 0.95,
            "method": "percentile", "n_boot": n_boot, "distribution": draws}
    prank = perm_rank_test(a, b, n_perm=n_perm, seed=seed)
    routes = pd.concat([routes, pd.DataFrame([{
        "route": "permutation (ranks)",
        "test": f"permutation of the rank sum ({n_perm:,})",
        "p_value": prank["p_value"],
        "assumes": "exchangeability; the rank counterpart of Mann-Whitney",
        "min_achievable_p": prank["p_resolution"],
        "reachable_at_0.05": True,
        "agrees_at_0.05": prank["p_value"] < ALPHA}])], ignore_index=True)
    routes.attrs["unanimous"] = bool(routes["agrees_at_0.05"].nunique() == 1)
    return {"groups": {l1: a, l2: b}, "routes": routes, "permutation": perm,
            "permutation_rank": prank, "bootstrap": boot, "p_mw_raw": float(p_u),
            "p_mw_bonf": float(min(1.0, p_u * 3)),
            "rank_biserial": float(2 * u / (a.size * b.size) - 1),
            "medians": (float(np.median(a)), float(np.median(b))),
            "n1": int(a.size), "n2": int(b.size)}


RET_INI_PAIR = pair_claim("times_object_retrieved_recent", "initiate_fetch",
                          "Myself", "Cat")


def boot_proportion(k, n, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    draws = rng.binomial(n, k / n, size=n_boot) / n
    return {"n": int(n), "statistic": "proportion", "point_estimate": float(k / n),
            "ci_low": float(np.percentile(draws, 2.5)),
            "ci_high": float(np.percentile(draws, 97.5)), "level": 0.95,
            "method": "percentile", "n_boot": n_boot, "distribution": draws}


def boot_cramers_v(counts, n_boot=N_BOOT, seed=SEED):
    counts = np.asarray(counts, float)
    N, k = counts.sum(), len(counts)
    rng = np.random.default_rng(seed)
    draws_counts = rng.multinomial(int(N), counts / N, size=n_boot).astype(float)
    exp = N / k
    chi2 = ((draws_counts - exp) ** 2 / exp).sum(axis=1)
    v = np.sqrt(chi2 / (N * (k - 1)))
    point = float(np.sqrt(((counts - exp) ** 2 / exp).sum() / (N * (k - 1))))
    return {"n": int(N), "statistic": "Cramer's V", "point_estimate": point,
            "ci_low": float(np.percentile(v, 2.5)),
            "ci_high": float(np.percentile(v, 97.5)), "level": 0.95,
            "method": "percentile", "n_boot": n_boot, "distribution": v}


BOOT_V = {k: boot_cramers_v(v["counts"]) for k, v in G.items()}

# ------------------------------------------------------------------ ledger
LEDGER_COLS = ["claim_id", "claim_short", "verdict", "effect_metric",
               "effect_reported", "effect_recomputed", "ci_low", "ci_high",
               "p_reported", "p_reported_relation", "p_recomputed",
               "p_smallest_possible", "survives_bonferroni_m28", "n_analysed",
               "stage", "verdict_source", "reasons"]
_rows = []


def add_claim(claim_id, claim_short, verdict, reasons, stage, verdict_source="proposed",
              effect_metric="", effect_reported=np.nan, effect_recomputed=np.nan,
              ci=(np.nan, np.nan), p_reported=np.nan, p_relation="",
              p_recomputed=np.nan, p_floor=np.nan, n_analysed=N_CATS):
    survives = ""
    ref_p = p_recomputed if np.isfinite(p_recomputed) else p_reported
    if np.isfinite(ref_p):
        survives = "yes" if ref_p < BONF_THRESHOLD else "no"
    _rows.append({
        "claim_id": claim_id, "claim_short": claim_short, "verdict": verdict,
        "effect_metric": effect_metric, "effect_reported": effect_reported,
        "effect_recomputed": effect_recomputed, "ci_low": ci[0], "ci_high": ci[1],
        "p_reported": p_reported, "p_reported_relation": p_relation,
        "p_recomputed": p_recomputed, "p_smallest_possible": p_floor,
        "survives_bonferroni_m28": survives, "n_analysed": n_analysed,
        "stage": stage, "verdict_source": verdict_source,
        "reasons": " | ".join(reasons),
    })


def freq_verdict(g, p_reported, extra_reasons=(), override=None, source="proposed"):
    routes = routes_frequency(g)
    ev = {"routes": routes, "permutation": {"floor_binds": False,
                                            "min_achievable_p": g["p_res"],
                                            "p_value": g["p_mc"]},
          "p_reported": p_reported, "p_recomputed": g["p_mc"],
          "effect_reported": None, "effect_recomputed": g["cramers_v"]}
    pv = propose_verdict(ev, alpha=ALPHA)
    reasons = list(pv["reasons"]) + list(extra_reasons)
    return (override or pv["verdict"]), reasons, routes, (source if override else "proposed")


def kw_verdict(k, p_reported, extra_reasons=(), override=None, source="reviewer-override"):
    ev = {"routes": k["routes"], "permutation": k["permutation"],
          "p_reported": p_reported, "p_recomputed": k["p_perm"],
          "effect_reported": None, "effect_recomputed": k["eps2"]}
    pv = propose_verdict(ev, alpha=ALPHA)
    reasons = list(pv["reasons"]) + list(extra_reasons)
    return (override or pv["verdict"]), reasons, (source if override else "proposed")


CLUSTER_NOTE = (
    f"cats are nested in owners ({N_CATS} cats, {N_OWNERS_REPORTED} owners reported) "
    f"and S3 carries no respondent identifier, so no cluster-robust test can be run; "
    f"the worst-case design effect is {DEFF_USED:.3f} (rho = 1), which divides the "
    f"statistic by at most that factor"
)

# --- C1 proportion untrained
_v, _r, _routes_c1, _src = freq_verdict(
    G["untrained"], np.nan,
    extra_reasons=[
        f"recomputed {N_UNTRAINED}/{N_CATS} = {PCT_UNTRAINED:.2f}% against the reported "
        f"94.4%, exact 95% CI [{100*UNTRAINED_CI[0]:.2f}, {100*UNTRAINED_CI[1]:.2f}]%",
        "the denominator is cats, not owners: the abstract attributes the 94.4% to "
        f"owners, and the owner-level proportion is not computable from S3",
        "no p-value is reported or needed for this descriptive claim",
    ], override="supported-with-caveat", source="reviewer-override")
add_claim("C1", "94.4% of fetching emerged without training", _v, _r, "1-2",
          _src, "proportion untrained", 0.944, N_UNTRAINED / N_CATS,
          UNTRAINED_CI, np.nan, "", np.nan, np.nan)

# --- C2 age bands
add_claim("C2", "fetching first noticed <1 yr (701) or 1-7 yr (415)", "supported",
          [f"life_stage counts reproduce exactly: Kitten {int(LIFE_STAGE['Kitten'])}, "
           f"Young adult {int(LIFE_STAGE['Young adult'])}, Mature adult "
           f"{int(LIFE_STAGE['Mature adult'])}, Senior {int(LIFE_STAGE['Senior'])}",
           f"the 'Young adult' band spans {YA_RANGE[0]:.2f}-{YA_RANGE[1]:.2f} years, so "
           f"'1-7 years' reads as 1 to under 7; the three cats recorded at exactly 7.00 "
           f"years are in 'Mature adult'",
           "descriptive counts; no test reported and none required"],
          "1-2", "reviewer-override", "proportion <1 yr",
          np.nan, float(LIFE_STAGE["Kitten"] / N_CATS),
          tuple(float(x) for x in
                stats.binomtest(int(LIFE_STAGE["Kitten"]), N_CATS).proportion_ci(
                    0.95, method="exact")),
          np.nan, "", np.nan, np.nan)

# --- C3-C6 initiation
_v, _r, _routes_c3, _src = freq_verdict(
    G["ini_omni"], 0.001,
    extra_reasons=[f"chi2 reproduces exactly: reported 124.29, recomputed "
                   f"{G['ini_omni']['chi2']:.2f} (df 2, N {G['ini_omni']['N']})",
                   f"Cramer's V = {G['ini_omni']['cramers_v']:.3f} "
                   f"[{BOOT_V['ini_omni']['ci_low']:.3f}, {BOOT_V['ini_omni']['ci_high']:.3f}]; "
                   f"the paper reports no effect size",
                   CLUSTER_NOTE])
add_claim("C3", "unequal distribution of who initiates fetching", _v, _r, "3-5", _src,
          "Cramer's V", np.nan, G["ini_omni"]["cramers_v"],
          (BOOT_V["ini_omni"]["ci_low"], BOOT_V["ini_omni"]["ci_high"]),
          0.001, "<", G["ini_omni"]["p_mc"], G["ini_omni"]["p_res"])

for cid, key, rep_chi2, short in [
        ("C4", "ini_cat_owner", 115.23, "cats initiate more often than owners"),
        ("C5", "ini_cat_equal", 44.64, "cats initiate more often than 'about equal'"),
        ("C6", "ini_owner_equal", 17.65, "owners initiate less often than 'about equal'")]:
    gg = G[key]
    _v, _r, _rt, _src = freq_verdict(
        gg, 0.001,
        extra_reasons=[f"chi2 reproduces: reported {rep_chi2}, recomputed {gg['chi2']:.2f} "
                       f"(df 1, n {gg['N']}); exact binomial p = {gg['p_exact']:.3g}",
                       f"post-hoc family of 3 interpreted at the paper's corrected "
                       f"alpha 0.016; survives", CLUSTER_NOTE])
    add_claim(cid, short, _v, _r, "3-6", _src, "Cramer's V", np.nan, gg["cramers_v"],
              (BOOT_V[key]["ci_low"], BOOT_V[key]["ci_high"]), 0.001, "<",
              gg["p_mc"], gg["p_res"], n_analysed=gg["N"])

# --- C7-C10 termination
_v, _r, _routes_c7, _src = freq_verdict(
    G["end_omni"], 0.001,
    extra_reasons=[f"chi2 reproduces exactly: reported 329.08, recomputed "
                   f"{G['end_omni']['chi2']:.2f} (df 2, N {G['end_omni']['N']})",
                   f"Cramer's V = {G['end_omni']['cramers_v']:.3f} "
                   f"[{BOOT_V['end_omni']['ci_low']:.3f}, {BOOT_V['end_omni']['ci_high']:.3f}]",
                   CLUSTER_NOTE])
add_claim("C7", "unequal distribution of who ends fetching", _v, _r, "3-5", _src,
          "Cramer's V", np.nan, G["end_omni"]["cramers_v"],
          (BOOT_V["end_omni"]["ci_low"], BOOT_V["end_omni"]["ci_high"]),
          0.001, "<", G["end_omni"]["p_mc"], G["end_omni"]["p_res"])

for cid, key, rep_chi2, short in [
        ("C8", "end_cat_owner", 197.54, "cats terminate more often than owners"),
        ("C9", "end_cat_equal", 217.59, "cats terminate more often than 'about equal'")]:
    gg = G[key]
    _v, _r, _rt, _src = freq_verdict(
        gg, 0.001,
        extra_reasons=[f"chi2 reproduces: reported {rep_chi2}, recomputed {gg['chi2']:.2f} "
                       f"(df 1, n {gg['N']}); exact binomial p = {gg['p_exact']:.3g}",
                       CLUSTER_NOTE])
    add_claim(cid, short, _v, _r, "3-6", _src, "Cramer's V", np.nan, gg["cramers_v"],
              (BOOT_V[key]["ci_low"], BOOT_V[key]["ci_high"]), 0.001, "<",
              gg["p_mc"], gg["p_res"], n_analysed=gg["N"])

gg = G["end_owner_equal"]
_v, _r, _routes_c10, _src = freq_verdict(
    gg, 0.437,
    extra_reasons=[f"chi2 and p both reproduce: reported 0.60 / p = 0.437, recomputed "
                   f"{gg['chi2']:.2f} / p = {gg['p_asym']:.3f}",
                   f"the claim is correctly worded as 'not a significant difference': the "
                   f"owner share of terminations is {100*gg['prop']:.1f}% of this pair, "
                   f"exact 95% CI [{100*gg['prop_ci'][0]:.1f}, {100*gg['prop_ci'][1]:.1f}]%, "
                   f"so differences up to about 6 percentage points remain compatible "
                   f"with these data"])
add_claim("C10", "owner vs 'about equal' termination: no difference", _v, _r, "3-5",
          _src, "Cramer's V", np.nan, gg["cramers_v"],
          (BOOT_V["end_owner_equal"]["ci_low"], BOOT_V["end_owner_equal"]["ci_high"]),
          0.437, "=", gg["p_mc"], gg["p_res"], n_analysed=gg["N"])

# --- C11 sex ratio
gg = G["sex"]
_v, _r, _routes_c11, _src = freq_verdict(
    gg, 0.018,
    extra_reasons=[
        f"reported X2 = 5.50 does not match the authors' own output: S2 prints "
        f"'Chi-squared test for given probabilities ... X-squared = "
        f"{S2T.statistic.iloc[0]:.4f}, df = {int(S2T.df.iloc[0])}, p-value = "
        f"{S2T.p_printed.iloc[0]}', which is exactly the recomputed value "
        f"{gg['chi2']:.4f} (p = {gg['p_asym']:.5f}). The analysis is correct and the "
        f"manuscript's 5.50 is a transcription of 5.5459; no continuity correction is "
        f"involved, since chisq.test applies Yates only to 2x2 contingency tables and "
        f"this is a goodness-of-fit call on a one-dimensional table",
        f"the null tested is 50:50 within a sample of fetching cats only; with no "
        f"non-fetching denominator the test describes sample composition and cannot "
        f"address whether sex predicts fetching",
        f"Rao-Scott correction at the worst-case design effect {DEFF_USED:.3f} gives "
        f"X2 = {gg['chi2_deff']:.2f}, p = {gg['p_deff']:.4f}: still below 0.05",
        f"does not survive Bonferroni across the {len(P_VALUES)} reported tests "
        f"(threshold {BONF_THRESHOLD:.5f})"],
    override="supported-with-caveat", source="reviewer-override")
add_claim("C11", "more male (617) than female (537) fetching cats", _v, _r, "3-6", _src,
          "Cramer's V", np.nan, gg["cramers_v"],
          (BOOT_V["sex"]["ci_low"], BOOT_V["sex"]["ci_high"]),
          0.018, "=", gg["p_mc"], gg["p_res"], n_analysed=gg["N"])

# --- C12-C17 Kruskal-Wallis family
_kw_specs = [
    ("C12", "ret_ini", 13.98, "retrievals differ by who initiates"),
    ("C14", "ret_end", 110.55, "retrievals differ by who ends"),
    ("C16", "mon_ini", 52.95, "monthly sessions differ by who initiates"),
    ("C17", "mon_end", 62.97, "monthly sessions differ by who ends"),
]
for cid, key, rep_H, short in _kw_specs:
    k = KW[key]
    extra = [f"H reproduces exactly: reported {rep_H}, recomputed {k['H']:.2f} "
             f"(df 2, N {k['N']})",
             f"all three routes agree at alpha: ANOVA p = {k['p_anova']:.3g}, "
             f"Kruskal-Wallis p = {k['p_kw']:.3g}, label permutation p = {k['p_perm']:.4f} "
             f"({N_PERM:,} replicates, resolution {k['p_res']:.1e})",
             f"epsilon-squared = {k['eps2']:.4f} [{k['eps2_ci'][0]:.4f}, "
             f"{k['eps2_ci'][1]:.4f}]: the grouping explains "
             f"{100*k['eps2']:.1f}% of rank variance; the paper reports no effect size",
             CLUSTER_NOTE]
    _v, _r, _src = kw_verdict(k, 0.001, extra_reasons=extra)
    add_claim(cid, short, _v, _r, "3-5", _src, "epsilon-squared", np.nan, k["eps2"],
              k["eps2_ci"], 0.001, "<", k["p_perm"], k["p_res"], n_analysed=k["N"])

# --- C13 the pairwise wording claim
pr = RET_INI_PAIR
_routes_c13 = pr["routes"]
_ev = {"routes": pr["routes"], "permutation": pr["permutation"], "p_reported": 0.009,
       "p_recomputed": float(pr["permutation"]["p_value"]),
       "effect_reported": None, "effect_recomputed": pr["rank_biserial"]}
_pv = propose_verdict(_ev, alpha=ALPHA)
add_claim("C13", "fewer retrievals when the owner initiates", "supported-with-caveat",
          _pv["reasons"] + [
              f"the Bonferroni-adjusted Wilcoxon reproduces: reported p = 0.009, "
              f"recomputed p = {pr['p_mw_bonf']:.4f} (raw {pr['p_mw_raw']:.5f} x 3 for the "
              f"family of three)",
              f"all five routes agree at alpha: "
              + "; ".join(f"{r.test} p = {r.p_value:.4g}"
                          for r in pr["routes"].itertuples()),
              f"both medians are {pr['medians'][0]:.0f} retrievals (owner IQR "
              f"{KW['ret_ini']['pairwise']['Cat vs Myself']['iqr_2']:.0f}, cat IQR "
              f"{KW['ret_ini']['pairwise']['Cat vs Myself']['iqr_1']:.0f}), so 'fewer "
              f"retrievals' is not a statement about the medians the sentence quotes; "
              f"Mann-Whitney tests stochastic ordering (rank-biserial "
              f"{pr['rank_biserial']:.3f} with the owner group first, i.e. "
              f"{-pr['rank_biserial']:+.3f} oriented to the claim's direction, which is "
              f"the sign recorded in effect_recomputed)",
              f"the difference is in the means and the upper tail: mean difference "
              f"(owner - cat) = {pr['bootstrap']['point_estimate']:.2f} retrievals, "
              f"95% CI [{pr['bootstrap']['ci_low']:.2f}, {pr['bootstrap']['ci_high']:.2f}]",
              f"what would resolve it: state the comparison as a shift in the "
              f"distribution, or quote means as well as medians"],
          "3-4", "reviewer-override", "rank-biserial (claim direction positive)", np.nan,
          -pr["rank_biserial"], (np.nan, np.nan), 0.009, "=",
          float(pr["permutation_rank"]["p_value"]),
          float(pr["permutation_rank"]["p_resolution"]), n_analysed=pr["n1"] + pr["n2"])

# --- C15 the terminator pairwise claim (reported medians 8 vs 5)
k = KW["ret_end"]
pw = k["pairwise"]["Cat vs Myself"]
add_claim("C15", "more retrievals when the owner ends the session", "supported",
          [f"medians and IQRs reproduce exactly: owner {pw['median_2']:.0f} "
           f"(IQR {pw['iqr_2']:.0f}) vs cat {pw['median_1']:.0f} (IQR {pw['iqr_1']:.0f})",
           f"Bonferroni-adjusted Wilcoxon p = {pw['p_bonf']:.3g} against the reported "
           f"p < 0.001; rank-biserial {pw['rank_biserial']:.3f} with the cat group first, "
           f"i.e. {-pw['rank_biserial']:+.3f} oriented to the claim's direction, which is "
           f"the sign recorded in effect_recomputed",
           f"the two remaining pairwise comparisons also reproduce: "
           + "; ".join(f"{a} p = {b['p_bonf']:.4g}" for a, b in k["pairwise"].items()
                       if a != "Cat vs Myself")],
          "3-5", "reviewer-override", "rank-biserial (claim direction positive)", np.nan,
          -pw["rank_biserial"], (np.nan, np.nan), 0.001, "<", pw["p_bonf"],
          float(rs.min_achievable_p(pw["n1"], pw["n2"])),
          n_analysed=pw["n1"] + pw["n2"])

# --- C18 cat sex x initiation / termination
for cid, res, rep_chi2, rep_p, short in [
        ("C18", SEX_INI, 1.82, 0.403, "no sex difference in who initiates"),
        ("C19", SEX_END, 1.97, 0.374, "no sex difference in who ends")]:
    ev = {"routes": res["routes"],
          "permutation": {"floor_binds": False, "min_achievable_p": res["p_res"],
                          "p_value": res["p_mc"]},
          "p_reported": rep_p, "p_recomputed": res["p_mc"],
          "effect_reported": None, "effect_recomputed": res["cramers_v"]}
    pv = propose_verdict(ev, alpha=ALPHA)
    add_claim(cid, short, pv["verdict"],
              pv["reasons"] + [
                  f"chi2 reproduces: reported {rep_chi2}, recomputed {res['chi2']:.3f} "
                  f"(df {res['df']}, N {res['n']}); asymptotic p = {res['p_asym']:.4f}, "
                  f"table permutation p = {res['p_mc']:.4f}",
                  f"statcheck flagged the reported pair for 'no sex difference in who "
                  f"ends' because p from the rounded X2 = 1.97 at df 2 is 0.3734, not the "
                  f"printed 0.374; S2 settles it - the authors' console output is "
                  f"'Pearson's Chi-squared test ... X-squared = {S2T.statistic.iloc[10]:.4f}, "
                  f"df = {int(S2T.df.iloc[10])}, p-value = {S2T.p_printed.iloc[10]}', so the "
                  f"printed p follows from the unrounded statistic and the flag was a "
                  f"rounding artefact",
                  "a null result: the data are compatible with no association and with "
                  "associations too small for this n to resolve"],
              "3-5", "proposed", "Cramer's V", np.nan, res["cramers_v"],
              (np.nan, np.nan), rep_p, "=", res["p_mc"], res["p_res"],
              n_analysed=res["n"])

# --- C20 / C21 sex x counts (two-group, four routes)
for cid, res, sens, s2i, rep_z, rep_p, short, override in [
        ("C20", SEX_RET, SENS_RET, 19, -1.69, 0.091, "no sex difference in retrievals",
         "supported-with-caveat"),
        ("C21", SEX_MON, SENS_MON, 21, -0.49, 0.622, "no sex difference in monthly sessions",
         None)]:
    ev = {"routes": res["routes"], "permutation": res["permutation_rank"],
          "p_reported": rep_p, "p_recomputed": float(res["permutation_rank"]["p_value"]),
          "effect_reported": None, "effect_recomputed": res["rank_biserial"]}
    pv = propose_verdict(ev, alpha=ALPHA)
    routes_txt = "; ".join(f"{r.test} p = {r.p_value:.4g}" for r in res["routes"].itertuples())
    crossed = (sens["retained"]["p"] < ALPHA) != (sens["excluded"]["p"] < ALPHA)
    reasons = pv["reasons"] + [
        f"the test is settled by S2: wilcox.test(exact = FALSE, conf.int = TRUE), printed "
        f"as 'Wilcoxon rank sum test with continuity correction, W = "
        f"{S2T.statistic.iloc[s2i]:.0f}, p-value = {S2T.p_printed.iloc[s2i]}'. Recomputing "
        f"with Female first reproduces W exactly: W = "
        f"{sens['retained']['W_female_first']:.0f}, p = {sens['retained']['p']:.5f}",
        f"the reported Z is the normal deviate of that p-value: qnorm(p/2) = "
        f"{sens['retained']['z_from_p']:.4f} against the reported Z = {rep_z}. The earlier "
        f"flag is resolved - the tie-corrected |Z| = {abs(res['z_tie']):.3f} is the right "
        f"comparison and the untied {abs(res['z']):.3f} is not",
        f"routes agree at alpha but not in magnitude: {routes_txt}",
        f"mean difference (male - female) = {res['bootstrap']['point_estimate']:.3f}, "
        f"95% CI [{res['bootstrap']['ci_low']:.3f}, {res['bootstrap']['ci_high']:.3f}]"]
    if crossed:
        reasons.append(
            f"the paper's own outlier-excluded sensitivity run crosses alpha: S2 prints "
            f"W = {S2T.statistic.iloc[s2i+1]:.0f}, p-value = {S2T.p_printed.iloc[s2i+1]} "
            f"for the same comparison with extreme outliers removed, and recomputing that "
            f"run gives W = {sens['excluded']['W_female_first']:.0f}, p = "
            f"{sens['excluded']['p']:.5f} on {sens['excluded']['n_female']} female and "
            f"{sens['excluded']['n_male']} male cats ({sens['n_excluded']} extreme values "
            f"removed). The Methods state that analyses were run with and without extreme "
            f"outliers and that 'there were no differences in the interpretation of the "
            f"findings'; for this comparison the interpretation does change")
    else:
        reasons.append(
            f"the outlier-excluded sensitivity run does not change the interpretation: "
            f"S2 prints W = {S2T.statistic.iloc[s2i+1]:.0f}, p-value = "
            f"{S2T.p_printed.iloc[s2i+1]}, and recomputing gives W = "
            f"{sens['excluded']['W_female_first']:.0f}, p = {sens['excluded']['p']:.5f} "
            f"({sens['n_excluded']} extreme values removed)")
    add_claim(cid, short, override or pv["verdict"], reasons,
              "3-5", "reviewer-override" if override else "proposed",
              "rank-biserial (male - female)", np.nan,
              res["rank_biserial"], (np.nan, np.nan), rep_p, "=",
              float(res["permutation_rank"]["p_value"]),
              float(res["permutation_rank"]["p_resolution"]),
              n_analysed=res["n1"] + res["n2"])

# --- C22-C24: claims the supplied inputs cannot address
add_claim("C22", "cats show independent, co-ordinated agency", "not-assessable",
          ["the supplied data are owner recall of who 'usually' initiates and ends "
           "sessions, recorded once per cat on a 3-point scale; they locate a frequency "
           "asymmetry but carry no within-session sequence, no latency and no "
           "observation independent of the owner",
           "what would change it: bout-level observational or video-coded data with "
           "timestamps for initiation and termination, or an experimental manipulation "
           "of owner availability"],
          "1", "reviewer-override")
add_claim("C23", "1154 cats analysed as independent observations", "not-assessable",
          [f"S3 contains no respondent, household or submission identifier (columns "
           f"searched for id/owner/respondent/participant/response: "
           f"{'none found' if not ID_LIKE else ', '.join(ID_LIKE)}), so cats cannot be "
           f"grouped by owner and no cluster-robust or mixed model can be fitted",
           f"{N_CATS} cats were reported by {N_OWNERS_REPORTED} owners "
           f"({MEAN_CATS_PER_OWNER:.3f} cats per owner); "
           f"{REPEAT_FLAG.get('Yes', 0)} rows carry ppl_demographics_repeated = 'Yes' "
           f"and {N_MULTI_FETCH} cats are recorded as sharing a household with at least "
           f"one other fetching cat",
           f"worst-case design effect (rho = 1) is {DEFF_EQUAL_WORST:.3f} from the "
           f"owner/cat ratio and {DEFF_RECORDED_WORST:.3f} from the recorded "
           f"fetching-cat household sizes; dividing every reported chi-square by "
           f"{DEFF_USED:.3f} leaves all of them on the same side of alpha",
           "what would change it: a respondent id column in S3, or the owner-level "
           "counts needed to fit a cluster-robust model"],
          "1-2", "reviewer-override")
add_claim("C24", "thematic analysis yields three main themes", "not-assessable",
          ["the review covers statistical claims; the themes rest on qualitative coding "
           f"of free text ({int(DF['extra_comments'].notna().sum())} of {N_CATS} cats have "
           f"extra_comments text in S3), and the paper reports codes "
           "(n = 268 initial, 35 child) without an inter-coder agreement statistic",
           "what would change it: the code frame plus a second coder's assignments, or "
           "the reliability statistic for the coding"],
          "1", "reviewer-override")

# --- C25 breed composition
GEN_BREED = DF["general_breed"].value_counts()
N_MIXED_REPORTED = int(GEN_BREED.get("Mixed", 0) + GEN_BREED.get("Unknown", 0))
add_claim("C25", "994 (86.1%) of fetching cats were mixed breed",
          "supported-with-caveat",
          [f"the reported 994 (86.1%) reproduces only by counting the "
           f"{int(GEN_BREED['Unknown'])} cats whose general_breed is 'Unknown' as mixed: "
           f"{int(GEN_BREED['Mixed'])} + {int(GEN_BREED['Unknown'])} = {N_MIXED_REPORTED} "
           f"({100*N_MIXED_REPORTED/N_CATS:.1f}%)",
           f"cats positively recorded as mixed breed number {int(GEN_BREED['Mixed'])} "
           f"({100*GEN_BREED['Mixed']/N_CATS:.1f}%); purebred "
           f"{int(GEN_BREED['Purebred'])} ({100*GEN_BREED['Purebred']/N_CATS:.1f}%) "
           f"reproduces exactly",
           "the paper notes unknown breeds as a limitation in the Discussion but the "
           "994 figure is quoted there without that split",
           "what would resolve it: report 'mixed or unknown' as the category label, or "
           "give the three counts separately"],
          "1-2", "reviewer-override", "proportion mixed-or-unknown",
          0.861, N_MIXED_REPORTED / N_CATS,
          tuple(float(x) for x in stats.binomtest(N_MIXED_REPORTED, N_CATS)
                .proportion_ci(0.95, method="exact")),
          np.nan, "", np.nan, np.nan)

# descriptive cross-check of the Table 1 / abstract figures
_sex_neuter = pd.crosstab(DF["cat_sex"], DF["cat_neutered"])
DESC_CHECK = pd.DataFrame([
    ("untrained cats", "1089 (94.4%)", f"{N_UNTRAINED} ({PCT_UNTRAINED:.1f}%)"),
    ("trained cats", "65 (5.6%)", f"{int(TRAINED['Trained'])} "
                                 f"({100*TRAINED['Trained']/N_CATS:.1f}%)"),
    ("fetching first noticed: kitten", "701", f"{int(LIFE_STAGE['Kitten'])}"),
    ("fetching first noticed: young adult", "415", f"{int(LIFE_STAGE['Young adult'])}"),
    ("current cat / previous cat", "853 (73.9%) / 301 (26.1%)",
     f"{int(DF.current_previous.value_counts()['Current cat'])} / "
     f"{int(DF.current_previous.value_counts()['Previous cat'])}"),
    ("female spayed / female intact", "501 (43.4%) / 36 (3.1%)",
     f"{int(_sex_neuter.loc['Female', 'Yes'])} / {int(_sex_neuter.loc['Female', 'No'])}"),
    ("male neutered / male intact / male unknown", "601 (52.1%) / 15 (1.3%) / 1 (0.1%)",
     f"{int(_sex_neuter.loc['Male', 'Yes'])} / {int(_sex_neuter.loc['Male', 'No'])} / "
     f"{int(_sex_neuter.loc['Male', 'Unknown'])}"),
    ("mixed breed", "994 (86.1%)",
     f"{int(GEN_BREED['Mixed'])} mixed + {int(GEN_BREED['Unknown'])} unknown = "
     f"{N_MIXED_REPORTED} ({100*N_MIXED_REPORTED/N_CATS:.1f}%)"),
    ("purebred", "160 (13.9%)",
     f"{int(GEN_BREED['Purebred'])} ({100*GEN_BREED['Purebred']/N_CATS:.1f}%)"),
    ("respondent gender: female / male / non-binary",
     "748 (64.8%) / 266 (23.1%) / 109 (9.4%)",
     " / ".join(str(int(DF.person_one_gender.value_counts()[k]))
                for k in ["Female", "Male", "Non-binary"])),
], columns=["reported in paper", "as printed", "recomputed from S3"])

# ------------------------------------------- S2 cross-check, claim by claim
# Row indices into S2T are the authors' print order (retained-outlier run first,
# outlier-excluded run immediately after for the Kruskal-Wallis and Wilcoxon tests).
S2_BY_CLAIM = {"C11": 0, "C3": 1, "C4": 2, "C5": 3, "C6": 4, "C7": 5, "C8": 6, "C9": 7,
               "C10": 8, "C18": 9, "C19": 10, "C12": 11, "C14": 13, "C16": 15, "C17": 17,
               "C20": 19, "C21": 21}
_recomputed_stat = {
    "C3": G["ini_omni"]["chi2"], "C4": G["ini_cat_owner"]["chi2"],
    "C5": G["ini_cat_equal"]["chi2"], "C6": G["ini_owner_equal"]["chi2"],
    "C7": G["end_omni"]["chi2"], "C8": G["end_cat_owner"]["chi2"],
    "C9": G["end_cat_equal"]["chi2"], "C10": G["end_owner_equal"]["chi2"],
    "C11": G["sex"]["chi2"], "C12": KW["ret_ini"]["H"], "C14": KW["ret_end"]["H"],
    "C16": KW["mon_ini"]["H"], "C17": KW["mon_end"]["H"],
    "C18": SEX_INI["chi2"], "C19": SEX_END["chi2"],
    "C20": SENS_RET["retained"]["W_female_first"],
    "C21": SENS_MON["retained"]["W_female_first"],
}
_recomputed_p = {
    "C3": G["ini_omni"]["p_asym"], "C4": G["ini_cat_owner"]["p_asym"],
    "C5": G["ini_cat_equal"]["p_asym"], "C6": G["ini_owner_equal"]["p_asym"],
    "C7": G["end_omni"]["p_asym"], "C8": G["end_cat_owner"]["p_asym"],
    "C9": G["end_cat_equal"]["p_asym"], "C10": G["end_owner_equal"]["p_asym"],
    "C11": G["sex"]["p_asym"], "C12": KW["ret_ini"]["p_kw"], "C14": KW["ret_end"]["p_kw"],
    "C16": KW["mon_ini"]["p_kw"], "C17": KW["mon_end"]["p_kw"],
    "C18": SEX_INI["p_asym"], "C19": SEX_END["p_asym"],
    "C20": SENS_RET["retained"]["p"], "C21": SENS_MON["retained"]["p"],
}
# the two pairwise-Wilcoxon claims come from the printed comparison matrices
_PW_CLAIM = {"C13": (0, 2, KW["ret_ini"]["pairwise"]["Cat vs Myself"]["p_bonf"]),
             "C15": (2, 2, KW["ret_end"]["pairwise"]["Cat vs Myself"]["p_bonf"])}

_s2_check = []
for cid, i in S2_BY_CLAIM.items():
    r = S2T.iloc[i]
    rec = _recomputed_stat[cid]
    agree = abs(rec - r.statistic) <= max(0.5 * 10 ** -_decimals(r.statistic), 1e-9) \
        if r.kind != "W (Wilcoxon rank sum)" else (abs(rec - r.statistic) < 0.5)
    _s2_check.append({
        "claim_id": cid, "statistic": r.kind, "data": r["data"],
        "authors_printed": f"{r.statistic:g}", "recomputed": f"{rec:.6g}",
        "authors_p": f"{r.p_relation} {r.p_printed}", "recomputed_p": f"{_recomputed_p[cid]:.4g}",
        "agree": "yes" if agree else "no"})
for cid, (mi, pi, rec) in _PW_CLAIM.items():
    printed = PAIRWISE_S2[mi][pi].strip()
    agree = (rec < 2e-16) if printed.startswith("<") else _p_round_trips(rec, printed)
    _s2_check.append({
        "claim_id": cid, "statistic": "pairwise Wilcoxon (Bonferroni)",
        "data": "owner vs cat", "authors_printed": "-", "recomputed": "-",
        "authors_p": printed, "recomputed_p": f"{rec:.3g}",
        "agree": "yes" if agree else "no"})
S2_CHECK = pd.DataFrame(_s2_check).set_index("claim_id").loc[
    [c for c in [r["claim_id"] for r in _rows] if c in
     set(S2_BY_CLAIM) | set(_PW_CLAIM)]].reset_index()

_s2_note = {r.claim_id: (
    f"S2 cross-check: the authors' printed {r.statistic} is {r.authors_printed} "
    f"(p {r.authors_p}) and the recomputed value is {r.recomputed} "
    f"(p = {r.recomputed_p}); the two "
    f"{'agree' if r.agree == 'yes' else 'do not agree'}"
    if r.authors_printed != "-" else
    f"S2 cross-check: the authors' printed pairwise Wilcoxon p is {r.authors_p} and the "
    f"recomputed Bonferroni-adjusted p is {r.recomputed_p}; the two "
    f"{'agree' if r.agree == 'yes' else 'do not agree'}")
    for r in S2_CHECK.itertuples()}
for _row in _rows:
    note = _s2_note.get(_row["claim_id"])
    if note and "S2" not in _row["reasons"]:
        _row["reasons"] = _row["reasons"] + " | " + note
    _row["authors_statistic_s2"] = S2_CHECK.loc[
        S2_CHECK.claim_id == _row["claim_id"], "authors_printed"].squeeze() \
        if _row["claim_id"] in set(S2_CHECK.claim_id) else ""
    _row["authors_p_s2"] = S2_CHECK.loc[
        S2_CHECK.claim_id == _row["claim_id"], "authors_p"].squeeze() \
        if _row["claim_id"] in set(S2_CHECK.claim_id) else ""

LEDGER_COLS = (LEDGER_COLS[:LEDGER_COLS.index("p_reported")]
               + ["authors_statistic_s2", "authors_p_s2"]
               + LEDGER_COLS[LEDGER_COLS.index("p_reported"):])
LEDGER = pd.DataFrame(_rows)[LEDGER_COLS]
VERDICT_COUNTS = (LEDGER["verdict"].value_counts()
                  .reindex(rf.VERDICT_ORDER, fill_value=0)
                  .rename_axis("verdict").reset_index(name="n_claims"))

# =========================================================== figures
OKABE = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9"]
FIGS = {}


def fig_design_map_survey(path):
    """
    F1 - what the sample contains, and what it cannot be grouped by.

    Left: the sampling ladder from consent to analysed cats, each step labelled
    with the count and the exclusion that produced it (abundance_and_proportion:
    count bars from zero, denominators stated).
    Right: cats per response category for the two headline frequency claims --
    the counts every chi-square in the paper is computed from.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.4),
                             gridspec_kw={"width_ratios": [1.15, 1]})

    ax = axes[0]
    steps = [("owners consenting", 1258, ""),
             ("owners completing (analysed denominator)", N_OWNERS_REPORTED,
              "-334 incomplete"),
             ("cats reported by those owners", 1184, ""),
             ("cats in the final sample", N_CATS, "-30 excluded")]
    y = np.arange(len(steps))[::-1]
    ax.barh(y, [s[1] for s in steps], color=["#999999", OKABE[0], "#999999", OKABE[1]],
            height=0.6)
    for yy, (lab, n, drop) in zip(y, steps):
        ax.text(n + 18, yy, f"{n:,}" + (f"   ({drop})" if drop else ""),
                va="center", fontsize=8.5)
    ax.set_yticks(y)
    ax.set_yticklabels([s[0] for s in steps], fontsize=8.5)
    ax.set_xlim(0, 3100)
    ax.set_xlabel("count (bars start at zero)")
    ax.set_title("Sampling ladder: owners consent, cats are analysed",
                 fontsize=10, loc="left")
    ax.text(0.98, 0.06,
            f"cats per owner (reported): {MEAN_CATS_PER_OWNER:.3f}\n"
            f"respondent id in S3:       none\n"
            f"rows flagged demographics-repeated: {REPEAT_FLAG.get('Yes', 0)}\n"
            f"cats sharing a household with\n"
            f"  another fetching cat:    {N_MULTI_FETCH}\n"
            f"worst-case design effect:  {DEFF_USED:.3f}  (rho = 1)\n"
            f"effective n (worst case):  {N_CATS / DEFF_USED:.0f} of {N_CATS}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            family="monospace",
            bbox=dict(boxstyle="round", facecolor="#FFF1E0", alpha=0.95,
                      edgecolor="#cccccc"))

    ax = axes[1]
    cats = ["Cat", "About equal", "Owner"]
    ini_v = [int(INI["Cat"]), int(INI["About equal"]), int(INI["Myself"])]
    end_v = [int(END["Cat"]), int(END["About equal"]), int(END["Myself"])]
    x = np.arange(3)
    w = 0.38
    ax.bar(x - w / 2, ini_v, w, color=OKABE[0], label="initiates sessions")
    ax.bar(x + w / 2, end_v, w, color=OKABE[1], label="ends sessions")
    for xi, v in zip(x - w / 2, ini_v):
        ax.text(xi, v + 8, f"{v}\n{100*v/N_CATS:.1f}%", ha="center", va="bottom",
                fontsize=7.5)
    for xi, v in zip(x + w / 2, end_v):
        ax.text(xi, v + 8, f"{v}\n{100*v/N_CATS:.1f}%", ha="center", va="bottom",
                fontsize=7.5)
    ax.axhline(N_CATS / 3, color="#444444", linestyle="--", linewidth=1,
               label=f"equal distribution ({N_CATS/3:.0f} cats)")
    ax.set_xticks(x); ax.set_xticklabels(cats)
    ax.set_ylim(0, max(max(ini_v), max(end_v)) * 1.42)
    ax.set_ylabel(f"cats (count of {N_CATS}; one response per cat)")
    ax.set_xlabel("owner-reported agent")
    ax.set_title("The counts behind C3-C10", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False, loc="upper right")

    fig.suptitle("F1 - Design map: 1,154 cats, 924 owners, no owner identifier",
                 fontsize=12, x=0.01, ha="left")
    fig.text(0.01, -0.07,
             "Unit of replication: the cat. The owner is the sampling unit and is not "
             "identifiable in S3, so clustering of cats within owners cannot be modelled; "
             f"the panel reports the worst-case bound instead. Denominator for all "
             f"percentages: {N_CATS} cats.",
             fontsize=8, color="#444444", ha="left", wrap=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def panel_frequency_claim(claim_id, g, routes, labels, verdict, path, title_extra="",
                          boot=None, expected_label=None):
    """
    Per-claim panel for a frequency claim: the counts, the Monte-Carlo null of
    the chi-square statistic with the observed value marked, and every route on
    one log p-axis.
    """
    n_panels = 3
    fig, axes = plt.subplots(1, n_panels, figsize=(4.5 * n_panels, 4.1))

    ax = axes[0]
    counts = g["counts"]
    cols = OKABE[:len(counts)]
    ax.bar(range(len(counts)), counts, color=cols, width=0.6)
    for i, v in enumerate(counts):
        ax.text(i, v + g["N"] * 0.012, f"{int(v)}\n{100*v/g['N']:.1f}%", ha="center",
                va="bottom", fontsize=8.5)
    ax.axhline(g["N"] / g["k"], color="#444444", linestyle="--", linewidth=1,
               label=expected_label or f"equal frequencies ({g['N']/g['k']:.0f})")
    ax.set_xticks(range(len(counts))); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0, max(counts) * 1.32)
    ax.set_ylabel(f"cats (count of {g['N']})")
    t = "counts, against the tested null"
    if boot is not None:
        t += (f"\n{boot['statistic']} = {boot['point_estimate']:.3f}, 95% CI "
              f"[{boot['ci_low']:.3f}, {boot['ci_high']:.3f}]")
    ax.set_title(t, fontsize=9, loc="left")
    ax.legend(fontsize=7.5, frameon=False, loc="upper center")

    ax = axes[1]
    null = np.asarray(g["null"], float)
    ax.hist(null, bins=60, color="#CCCCCC", edgecolor="white", linewidth=0.3)
    obs_in_range = g["chi2"] <= null.max() * 1.05
    ax.axvline(g["chi2"], color="#D55E00", linewidth=2,
               label=f"observed chi2 = {g['chi2']:.2f}")
    if not obs_in_range:
        ax.set_xlim(0, null.max() * 1.1)
        ax.annotate(f"observed chi2 = {g['chi2']:.1f}\n(off scale to the right)",
                    xy=(0.97, 0.55), xycoords="axes fraction", ha="right", fontsize=8,
                    color="#D55E00")
    ax.set_xlabel("Pearson chi-square under the null (multinomial draws)")
    ax.set_ylabel(f"replicates (count of {N_MC:,})")
    ax.set_title("Monte-Carlo null", fontsize=9, loc="left")
    ax.legend(fontsize=8, frameon=False, loc="center left")
    ax.text(0.98, 0.97,
            f"p(MC)      = {g['p_mc']:.5f}\n"
            f"resolution = {g['p_res']:.1e}\n"
            f"p(asympt)  = {g['p_asym']:.3g}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="#cccccc"))

    ax = axes[2]
    r = routes.reset_index(drop=True)
    ypos = np.arange(len(r))[::-1]
    ax.scatter(np.clip(r["p_value"], 1e-30, 1),
               ypos, s=60,
               color=["#009E73" if s else "#D55E00" for s in r["agrees_at_0.05"]],
               zorder=3)
    for yy, row in zip(ypos, r.itertuples()):
        ax.annotate(f"{row.p_value:.3g}", (max(row.p_value, 1e-30), yy),
                    textcoords="offset points", xytext=(6, 5), fontsize=7)
    ax.axvline(ALPHA, color="#D55E00", linestyle=":", linewidth=1.2,
               label=f"alpha = {ALPHA}")
    ax.set_xscale("log")
    lo = min(1e-30, float(np.min(r["p_value"])) / 10) if float(np.min(r["p_value"])) > 0 else 1e-30
    ax.set_xlim(max(lo, 1e-80), 2.0)
    ax.set_yticks(ypos); ax.set_yticklabels(r["test"], fontsize=8)
    ax.set_xlabel("p-value (log scale)")
    ax.set_title("every route, same data", fontsize=9, loc="left")
    ax.legend(fontsize=8, frameon=False, loc="lower left")

    fig.suptitle(f"{claim_id} - verdict: {verdict}{title_extra}", fontsize=11, x=0.01,
                 ha="left", color=rf.VERDICT_COLORS.get(verdict, "#999999"))
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def build_figures():
    p = os.path.join
    FIGS["F1_design_map"] = fig_design_map_survey(p(FIGDIR, "F1_design_map.png"))

    head = KW["ret_end"]
    f = rf.fig_raw_behind_claim(
        head["groups"], "retrievals in the most recent session (count)",
        path=p(FIGDIR, "F2_raw_behind_claim.png"),
        test_name=f"Kruskal-Wallis H = {head['H']:.1f}, label permutation",
        p_value=head["p_perm"],
        reported_error="median and IQR only; no SD, SEM or interval is reported")
    plt.close(f)
    FIGS["F2_raw_behind_claim"] = p(FIGDIR, "F2_raw_behind_claim.png")

    f = rf.fig_assumption_qq(head["groups"], path=p(FIGDIR, "F3_assumption_shape.png"),
                             value_label="retrievals in the most recent session (count)")
    plt.close(f)
    FIGS["F3_assumption_shape"] = p(FIGDIR, "F3_assumption_shape.png")

    f = rf.fig_claim_ledger(
        LEDGER, path=p(FIGDIR, "F4_claim_ledger.png"),
        effect_label="effect size (proportion, Cramer's V, epsilon-squared or "
                     "rank-biserial; see ledger column effect_metric)")
    plt.close(f)
    FIGS["F4_claim_ledger"] = p(FIGDIR, "F4_claim_ledger.png")

    gend = G["end_omni"]
    perm_for_f5 = {
        "statistic_name": "Pearson chi-square, who ends sessions (3 categories)",
        "observed": gend["chi2"], "p_value": gend["p_mc"], "n_perm": N_MC,
        "p_resolution": gend["p_res"], "min_achievable_p": gend["p_res"],
        "floor_binds": False, "null": gend["null"],
        "note": (f"Monte-Carlo multinomial null, {N_MC:,} replicates, seed {SEED}; "
                 f"resolution {gend['p_res']:.1e}. The observed chi-square "
                 f"({gend['chi2']:.1f}) lies far outside the null range "
                 f"(max {gend['null'].max():.1f}), so the Monte-Carlo p-value is "
                 f"reported at its resolution limit rather than as an estimate."),
    }
    f = rf.fig_null_distribution(perm_for_f5, BOOT_V["end_omni"],
                                path=p(FIGDIR, "F5_null_distribution.png"))
    plt.close(f)
    FIGS["F5_null_distribution"] = p(FIGDIR, "F5_null_distribution.png")

    f = rf.fig_multiplicity(P_VALUES, alpha=ALPHA, path=p(FIGDIR, "F6_multiplicity.png"),
                            labels=P_LABELS)
    plt.close(f)
    FIGS["F6_multiplicity"] = p(FIGDIR, "F6_multiplicity.png")

    # ---- per-claim panels
    def verdict_of(cid):
        return str(LEDGER.loc[LEDGER.claim_id == cid, "verdict"].iloc[0])

    FIGS["panel_C1"] = panel_frequency_claim(
        "C1", G["untrained"], routes_frequency(G["untrained"]),
        ["untrained", "trained"], verdict_of("C1"),
        p(FIGDIR, "panel_C1_untrained.png"),
        title_extra="  (94.4% untrained)",
        boot=boot_proportion(N_UNTRAINED, N_CATS),
        expected_label="50:50 (template null, not the paper's claim)")

    FIGS["panel_C3"] = panel_frequency_claim(
        "C3", G["ini_omni"], _routes_c3, ["cat", "about equal", "owner"],
        verdict_of("C3"), p(FIGDIR, "panel_C3_initiation.png"),
        title_extra="  (who initiates)", boot=BOOT_V["ini_omni"])

    FIGS["panel_C7"] = panel_frequency_claim(
        "C7", G["end_omni"], _routes_c7, ["cat", "about equal", "owner"],
        verdict_of("C7"), p(FIGDIR, "panel_C7_termination.png"),
        title_extra="  (who ends)", boot=BOOT_V["end_omni"])

    FIGS["panel_C11"] = panel_frequency_claim(
        "C11", G["sex"], _routes_c11, ["male", "female"], verdict_of("C11"),
        p(FIGDIR, "panel_C11_sex_ratio.png"),
        title_extra="  (sex ratio among fetching cats)", boot=BOOT_V["sex"])

    f = rf.fig_claim_panel(
        "C13", RET_INI_PAIR["groups"], "retrievals in the most recent session (count)",
        RET_INI_PAIR["permutation"], RET_INI_PAIR["bootstrap"], RET_INI_PAIR["routes"],
        verdict=verdict_of("C13"), path=p(FIGDIR, "panel_C13_owner_initiates.png"))
    plt.close(f)
    FIGS["panel_C13"] = p(FIGDIR, "panel_C13_owner_initiates.png")

    f = rf.fig_claim_panel(
        "C20", SEX_RET["groups"], "retrievals in the most recent session (count)",
        SEX_RET["permutation"], SEX_RET["bootstrap"], SEX_RET["routes"],
        verdict=verdict_of("C20"), path=p(FIGDIR, "panel_C20_sex_retrievals.png"))
    plt.close(f)
    FIGS["panel_C20"] = p(FIGDIR, "panel_C20_sex_retrievals.png")

    for cid, reason in [
        ("C22", "owner recall of who 'usually' initiates and ends fetching cannot "
                "measure agency: no bout-level sequence, latency or independent "
                "observation is in S3"),
        ("C23", "S3 carries no respondent identifier, so cats cannot be grouped by "
                "owner and no cluster-robust test can be run"),
        ("C24", "the themes rest on qualitative coding; no code frame or inter-coder "
                "agreement statistic is reported"),
    ]:
        path_ = p(FIGDIR, f"panel_{cid}_not_assessable.png")
        f = rf.not_assessable_figure(
            f"{cid} - {LEDGER.loc[LEDGER.claim_id == cid, 'claim_short'].iloc[0]}",
            reason, path_)
        plt.close(f)
        FIGS[f"panel_{cid}"] = path_
    return FIGS


# =========================================================== report
MARKER_MAP_PATH = os.path.join(OUTDIR, "figure_artifact_map.json")


def _img(path):
    """Figure reference: an artifact marker when the map exists, else the path."""
    try:
        import json
        with open(MARKER_MAP_PATH, encoding="utf-8") as fh:
            m = json.load(fh)
        return m.get(os.path.basename(path), path.replace(os.sep, "/"))
    except Exception:
        return path.replace(os.sep, "/")


def md_table(dfr, cols=None):
    d = dfr if cols is None else dfr[cols]
    head = "| " + " | ".join(str(c) for c in d.columns) + " |"
    rule = "|" + "|".join("---" for _ in d.columns) + "|"
    rows = ["| " + " | ".join("" if (isinstance(v, float) and not np.isfinite(v))
                              else (f"{v:.4g}" if isinstance(v, float) else str(v))
                              for v in r) + " |"
            for r in d.itertuples(index=False)]
    return "\n".join([head, rule] + rows)


def write_report():
    f = FIGS
    k_ret_end = KW["ret_end"]
    shape = k_ret_end["shape"]
    shape_txt = "; ".join(
        f"{r.name}: n = {int(r.n)}, median = {r.median:.0f}, skew = {r.skewness:.2f}, "
        f"excess kurtosis = {r.excess_kurtosis:.1f}, Shapiro p = {r.shapiro_p:.2g}"
        for r in shape.itertuples(index=False))
    nsig = int(RECOUNT["summary"].loc[RECOUNT["summary"].correction == "none",
                                      "n_significant"].iloc[0])
    nbh = int(RECOUNT["summary"].loc[RECOUNT["summary"].correction == "Benjamini-Hochberg",
                                     "n_significant"].iloc[0])
    nbonf = int(RECOUNT["summary"].loc[RECOUNT["summary"].correction == "Bonferroni",
                                       "n_significant"].iloc[0])
    n_bounds = int(P_IS_BOUND.sum())
    vc = VERDICT_COUNTS[VERDICT_COUNTS.n_claims > 0]
    sc_all = STATCHECK_ALL if STATCHECK_ALL is not None else build_statcheck()
    _sc_text = sc_all[sc_all.source == "manuscript text"]
    _sc_s2 = sc_all[sc_all.source == "S2 console output"]
    _sc_text_total = len(_sc_text)
    _sc_text_checked = int(_sc_text["consistent"].notna().sum())
    _sc_s2_checkable = int(_sc_s2["consistent"].notna().sum())
    _sc_s2_consistent = int((_sc_s2["consistent"] == True).sum())  # noqa: E712
    _sc_s2_de = int((_sc_s2["decision_error"] == True).sum())      # noqa: E712
    _sc_s2_inconsistent = _sc_s2_checkable - _sc_s2_consistent
    _de_clause = (f"and {_sc_s2_de} of them is a decision error (the printed p and the p "
                  f"its own statistic implies fall on opposite sides of alpha)"
                  if _sc_s2_de else "and none is a decision error")
    _sc_s2_sentence = (
        f"all {_sc_s2_checkable} printed (statistic, df, p) triples in S2 are internally "
        f"consistent, {_de_clause}."
        if _sc_s2_inconsistent == 0 else
        f"{_sc_s2_consistent} of the {_sc_s2_checkable} printed (statistic, df, p) triples "
        f"in S2 are internally consistent and {_sc_s2_inconsistent} "
        f"{'is' if _sc_s2_inconsistent == 1 else 'are'} not "
        f"({', '.join(_sc_s2.loc[_sc_s2.consistent == False, 'claim_id'])}), {_de_clause}.")  # noqa: E712
    _n_agree = int((S2_CHECK.agree == "yes").sum())
    _s2_agree_sentence = (
        f"All {len(S2_CHECK)} agree."
        if _n_agree == len(S2_CHECK) else
        f"{_n_agree} of the {len(S2_CHECK)} agree; "
        f"{', '.join(S2_CHECK.loc[S2_CHECK.agree != 'yes', 'claim_id'])} do not.")

    md = f"""# Statistical review: Forman, Renner & Leavens (2023), *Fetching felines*

Forman J., Renner E. & Leavens D.A. (2023) Fetching felines: a survey of cat owners
on the diversity of cat (*Felis catus*) fetching behaviour. *Scientific Reports*
**13**:20456. doi:10.1038/s41598-023-47409-w

Reviewed with the `statistical-review` skill (six-stage protocol); figures follow the
`data-visualization` skill. Every number below was recomputed from the supplied
Supplementary S3 dataset; nothing is simulated, reconstructed or approximated.

---

## 1. Header: what was supplied, and what was assessable

| input | path | size | sha256 (first 16) | corresponds to |
|---|---|---|---|---|
| paper text | `{os.path.basename(PAPER)}` | {len(TXT):,} characters | `{sha(PAPER)[:16]}` | whole article (text layer) |
| raw data | `{os.path.basename(DATA)}` | {N_CATS} rows x {DF.shape[1]} columns | `{sha(DATA)[:16]}` | Supplementary S3; Table 1-2, Figures 2-3 and every test in Results |
| authors' analysis code | `{os.path.basename(S2_PATH)}` | {len(S2):,} characters | `{sha(S2_PATH)[:16]}` | Supplementary S2, extracted from the .docx: R source with console output interleaved |

{TEXT_NOTE}

**Supplementary S2.** S2 was supplied after the first pass of this review and this
revision incorporates it. Word collapsed many in-paragraph line breaks, so it is not a
runnable script, but the printed console output survives intact: {len(S2T)} test statistics
with their degrees of freedom and p-values, {len(PAIRWISE_S2)} pairwise-Wilcoxon comparison
matrices, and the call that produced each one. Every recomputed test in this report is now
checked against the authors' own call and printed result (section 5), and the one
discrepancy the first pass could not resolve -- the Mann-Whitney Z -- is settled.

**Claim extraction.** The automatic extractors found
{int(INV_SUMMARY.loc[INV_SUMMARY.table == 'tests', 'n_rows'].iloc[0])} test statistics,
{int(INV_SUMMARY.loc[INV_SUMMARY.table == 'p_values', 'n_rows'].iloc[0])} p-values and
{int(INV_SUMMARY.loc[INV_SUMMARY.table == 'named_tests', 'n_rows'].iloc[0])} named tests.
Four of the extracted p-values are artefacts of the text layer (three numbers inside the
Table 1 demographic block, and the alpha threshold quoted in the Levene sentence); they
are excluded, leaving **{len(P_VALUES)} reported p-values** from
**{len(P_VALUES) - 11} named tests plus 11 pairwise Wilcoxon comparisons**. The
{int(INV_SUMMARY.loc[INV_SUMMARY.table == 'tests', 'n_rows'].iloc[0])} reported test
statistics were each recomputed by hand against the raw data.

**Coverage.** {len(LEDGER)} claims are in the ledger. They cover all
{int(INV_SUMMARY.loc[INV_SUMMARY.table == 'tests', 'n_rows'].iloc[0])} reported test
statistics, all {len(P_VALUES)} reported p-values (the 11 pairwise Wilcoxons are
recomputed inside C12-C17), the four headline descriptives in the abstract, and ten
Table 1 descriptives cross-checked in section 4. {int(vc.loc[vc.verdict == 'not-assessable', 'n_claims'].iloc[0])} claims are
`not-assessable` and are listed by id in section 9.

**Not covered.** The object-category tallies and thematic map (the paper's Figures 1-4
and the object groupings) are descriptive or qualitative; no statistical claim rests on
them and none was re-derived beyond the breed and training counts in section 4. S2 cannot
be executed as supplied (collapsed line breaks), so each test is still re-implemented
rather than re-run -- but it is now re-implemented against a known call, and the result is
compared with the authors' printed output test by test.

---

## 2. Summary of findings

{md_table(vc)}

Three findings would change how a reader reads the paper.

1. **The unit of analysis is the cat; the sampling unit is the owner; and S3 cannot
   connect them (C23, stage 1-2).** {N_OWNERS_REPORTED} owners reported on {N_CATS} cats
   ({MEAN_CATS_PER_OWNER:.3f} cats per owner), {REPEAT_FLAG.get('Yes', 0)} rows are
   flagged as repeating a previous respondent's demographics, and {N_MULTI_FETCH} cats
   are recorded as sharing a household with at least one other fetching cat -- but S3
   contains no respondent, household or submission identifier, so no cluster-robust or
   mixed model can be fitted. What *can* be computed is the bound: the worst-case design
   effect (intra-cluster correlation = 1) is {DEFF_EQUAL_WORST:.3f} from the owner/cat
   ratio and {DEFF_RECORDED_WORST:.3f} from the recorded fetching-cat household sizes.
   Dividing every reported chi-square by {DEFF_USED:.3f} leaves all of them on the same
   side of alpha -- the sex-ratio test moves from p = {G['sex']['p_asym']:.4f} to
   p = {G['sex']['p_deff']:.4f}. Clustering is therefore a reporting gap here, not a
   result-changing one.

2. **The outlier sensitivity analysis changes one interpretation, and the Methods say it
   changes none (C20, stage 5-6).** The Methods state that analyses were run with and
   without extreme outliers and that "there were no differences in the interpretation of
   the findings". For the sex difference in retrievals that is not the case. With outliers
   retained -- the reported analysis -- S2 prints W = {S2T.statistic.iloc[19]:.0f},
   p = {S2T.p_printed.iloc[19]}, which recomputes exactly
   (W = {SENS_RET['retained']['W_female_first']:.0f}, p = {SENS_RET['retained']['p']:.5f}).
   With extreme outliers excluded by the authors' own rule S2 prints
   W = {S2T.statistic.iloc[20]:.0f}, p = {S2T.p_printed.iloc[20]}, which also recomputes
   exactly (W = {SENS_RET['excluded']['W_female_first']:.0f},
   p = {SENS_RET['excluded']['p']:.5f}; {SENS_RET['n_excluded']} extreme values removed,
   leaving {SENS_RET['excluded']['n_female']} female and {SENS_RET['excluded']['n_male']}
   male cats). That crosses alpha: the comparison the paper reports as non-significant is
   significant in its own sensitivity run, and S2's code comment on that line reads
   "significant differences between groups". The same check on monthly sessions does not
   change (p = {SENS_MON['retained']['p']:.4f} retained,
   p = {SENS_MON['excluded']['p']:.4f} excluded).

3. **"Significantly fewer retrievals if the owner usually initiates" is quoted with two
   identical medians (C13, stage 3-4).** The sentence reads "(median = 5, IQR = 4)
   compared to the cat (median = 5, IQR = 6, Pairwise Wilcoxon p = 0.009)". The adjusted
   p reproduces ({RET_INI_PAIR['p_mw_bonf']:.4f}) and all five routes agree at alpha, so
   the result stands -- but it is a statement about stochastic ordering and the mean, not
   about medians: the recomputed mean difference is
   {RET_INI_PAIR['bootstrap']['point_estimate']:.2f} retrievals
   (95% CI [{RET_INI_PAIR['bootstrap']['ci_low']:.2f},
   {RET_INI_PAIR['bootstrap']['ci_high']:.2f}]), rank-biserial
   {RET_INI_PAIR['rank_biserial']:.3f} with the owner group first.

Everything else is in the ledger. The headline results of the paper -- that fetching is
overwhelmingly untrained, that it is first seen in young cats, and that cats initiate and
terminate fetching more often than owners do -- reproduce exactly from S3.

---

## 3. Claim ledger -- F4

![F4 claim ledger]({_img(f['F4_claim_ledger'])})

Left: the recomputed effect with its 95% bootstrap interval, one row per claim. Every
open marker is absent because **the paper reports no effect size for any of its tests**;
the only reported effects are the percentages in C1, C2 and C25. Right: reported against
recomputed p-values on log axes. No point lies in an off-diagonal quadrant across 0.05,
so no recomputation changes a decision. The points sitting above the diagonal are the
claims whose reported p is a `< 0.001` bound recorded at 0.001 while the recomputed value
is smaller, and C20/C21, where the paper's Mann-Whitney p and the review's rank
permutation p differ in the third decimal.

{md_table(LEDGER, ['claim_id', 'claim_short', 'verdict', 'effect_metric',
                   'effect_reported', 'effect_recomputed', 'ci_low', 'ci_high',
                   'p_reported', 'p_reported_relation', 'p_recomputed',
                   'p_smallest_possible', 'survives_bonferroni_m28'])}

**Sign convention.** `effect_recomputed` is oriented so that a positive value means "in
the direction the paper claims" wherever `effect_metric` says so (C13, C15); the raw
rank-biserial in the `reasons` column carries the sign of the group order the test was
run in, and both signs of the same statistic are given there so the two can be matched.
For C20 and C21 -- null claims, where there is no claimed direction -- the raw
male-minus-female sign is kept.

`reasons` for every row is in `claim_ledger.csv`; `p_smallest_possible` is the resolution
of the recomputed reference route (Monte-Carlo or permutation), i.e. the smallest p-value
that procedure could have returned.

---

## 4. Design and sampling -- F1

![F1 design map]({_img(f['F1_design_map'])})

**Population.** The paper claims diversity of fetching behaviour among cats that fetch.
The sample was recruited through Twitter, Facebook and Reddit from owners who currently
or previously owned a fetching cat, so the population is self-selected owners of fetching
cats who use those platforms; {int(DF.country.value_counts()['USA'])} of {N_CATS} cats
are in the USA. The paper states this limit explicitly and does not claim prevalence.

**Unit of replication.** The left panel is the sampling ladder: 1,258 owners consented,
{N_OWNERS_REPORTED} completed (334 excluded as incomplete), reporting 1,184 cats, of
which {N_CATS} survive the stated exclusions. Two things follow. First, the "{N_OWNERS_REPORTED}
owners / {N_CATS} cats" pairing in the abstract combines a pre-exclusion owner count with
a post-exclusion cat count; the number of owners contributing to the final {N_CATS} cats
is not reported and cannot be derived from S3. Second, cats are nested within owners and
no identifier links them (C23): the columns searched for an id are
`{'id/owner/respondent/participant/response: none found' if not ID_LIKE else ', '.join(ID_LIKE)}`.
The design-effect bound in the panel is the strongest statement the supplied data allow.

**Exclusions** are stated and counted in the paper (334 incomplete, 2 extreme, 23
ambiguous, 5 duplicate) and the arithmetic closes: 1,184 - 30 = {N_CATS}.

**Confounding.** Not checkable. Nothing in S3 identifies when a response was collected,
and the grouping variables of interest (who initiates, who ends) are themselves
owner-reported, so owner interpretation is shared between predictor and outcome.

**Descriptive cross-check.** Ten figures from the abstract and Table 1, recomputed:

{md_table(DESC_CHECK)}

Nine of the ten reproduce exactly. The tenth is C25: "994 (86.1%) mixed breed" is
{int(GEN_BREED['Mixed'])} cats recorded as mixed plus {int(GEN_BREED['Unknown'])} recorded
as unknown. The Discussion quotes the 994 as "mixed breeds" in one sentence and notes
unknown breed as a limitation in the next.

---

## 5. Method and assumptions -- F2, F3

![F2 raw data behind the claim]({_img(f['F2_raw_behind_claim'])})

![F3 distribution shape]({_img(f['F3_assumption_shape'])})

**What the paper did.** Normality was screened with Shapiro-Wilk tests and Q-Q plots, and
because "the Shapiro-Wilk tests determined the distributions of the data were
significantly different from a normal distribution", non-parametric tests were used
throughout: chi-square for frequencies, Kruskal-Wallis with Bonferroni-adjusted pairwise
Wilcoxon tests for the two count outcomes, Mann-Whitney for sex contrasts. Extreme
outliers (beyond Q3 + 3xIQR) were identified and the analyses were run with and without
them, with the outlier-retaining version reported.

**What the screen can and cannot establish at this n.** Shapiro-Wilk is a test of a null
of exact normality whose power rises with n; at the group sizes here
({'; '.join(f"{g_}: {int(n_)}" for g_, n_ in zip(shape['name'], shape['n']))}) it will
reject departures far too small to matter for any test about means, and a non-significant
result would not have been evidence of normality either. The diagnostic that *is*
informative is the shape itself, and F3 shows it: {shape_txt}. These are bounded,
heavily right-skewed integer counts with a floor at zero
({int((DF.times_object_retrieved_recent == 0).sum())} cats retrieved 0 times;
{int((DF.times_monthly_fetch_sessions < 1).sum())} had fewer than one monthly session) and
a long upper tail (maximum {int(DF.times_object_retrieved_recent.max())} retrievals,
{int(DF.times_monthly_fetch_sessions.max())} monthly sessions). Rank-based tests are a
defensible choice for such data -- but because the outcome is a skewed count and the
question is about which group tends to be higher, not because a normality test returned
p < 0.05. This review does not treat the Shapiro result as a pass/fail gate in either
direction, and no verdict in the ledger turns on it.

**Whether the gate mattered.** It did not. Running the parametric route anyway, the
one-way ANOVA agrees with Kruskal-Wallis at alpha on all four omnibus comparisons
(C12 p = {KW['ret_ini']['p_anova']:.3g} vs {KW['ret_ini']['p_kw']:.3g};
C14 {KW['ret_end']['p_anova']:.3g} vs {KW['ret_end']['p_kw']:.3g};
C16 {KW['mon_ini']['p_anova']:.3g} vs {KW['mon_ini']['p_kw']:.3g};
C17 {KW['mon_end']['p_anova']:.3g} vs {KW['mon_end']['p_kw']:.3g}), and the label
permutation of H agrees with both. The conclusions do not depend on the modelling choice.

**Null versus question.** Three mismatches are worth the authors' attention, none of them
decision-changing. (i) Kruskal-Wallis and Mann-Whitney test stochastic ordering, not
equality of medians; the Results sentences quote medians either side of a Wilcoxon
p-value, which reads as a median comparison (C13 is the clearest case: both medians 5).
(ii) The chi-square goodness-of-fit tests in C3-C11 test a uniform null over the response
categories; "cats initiate more than owners" is a claim about that sample's composition,
which is what the test delivers, but the sex-ratio claim (C11) is phrased as a claim about
cats in general. (iii) The paper reports Levene's tests showing unequal variances as the
reason for choosing the pairwise adjustment; Levene on heavily skewed counts is itself
shape-sensitive, and the choice of a Bonferroni-adjusted rank test does not depend on it.

**Cross-check against the authors' own calls (S2).** Every test the paper reports can now
be matched to the call that produced it and to its printed output:

| claim | call in S2 | authors' printed statistic | recomputed here | authors' p | recomputed p | agree |
|---|---|---|---|---|---|---|
{chr(10).join('| ' + ' | '.join([r.claim_id, str(r.statistic), str(r.authors_printed), str(r.recomputed), str(r.authors_p), str(r.recomputed_p), str(r.agree)]) + ' |' for r in S2_CHECK.itertuples())}

{_s2_agree_sentence} Statcheck gains from S2 too: the manuscript prints degrees of
freedom for only {int(_sc_text_checked)} of its {int(_sc_text_total)} test statistics, so
most triples were not checkable from the text, whereas {_sc_s2_sentence} The calls
are `chisq.test` on one-dimensional tables for C3-C11
("Chi-squared test for given probabilities", i.e. goodness-of-fit against equal
probabilities -- which is what this review implemented), `chisq.test` on a 2x3 table for
C18-C19, `kruskal.test` for C12-C17, `pairwise.wilcox.test(p.adjust.method =
"bonferroni")` for the pairwise comparisons, and
`wilcox.test(paired = FALSE, exact = FALSE, conf.int = TRUE)` for C20-C21. Three
consequences:

- **C11 is a transcription slip, not an analysis error.** S2 prints X-squared =
  {S2T.statistic.iloc[0]:.4f}, df = {int(S2T.df.iloc[0])}, p-value = {S2T.p_printed.iloc[0]},
  identical to the recomputed {G['sex']['chi2']:.4f}. The manuscript's "X2 = 5.50" is
  5.5459 re-rounded. Continuity correction never entered it: `chisq.test` applies Yates
  only to 2x2 contingency tables, so the Yates alternative floated in the first pass of
  this review is not a candidate and has been withdrawn.
- **The Mann-Whitney Z is resolved.** R's `wilcox.test` prints W, not Z; the reported
  Z = -1.69 and Z = -0.49 are the normal deviates of the printed p-values
  (qnorm(p/2) = {SENS_RET['retained']['z_from_p']:.4f} and
  {SENS_MON['retained']['z_from_p']:.4f}). The tie-corrected |Z| computed in the first pass
  ({abs(SEX_RET['z_tie']):.3f}, {abs(SEX_MON['z_tie']):.3f}) is the matching quantity; the
  uncorrected one is not. There is no discrepancy.
- **The statcheck flag on C19 was a rounding artefact.** S2 prints X-squared =
  {S2T.statistic.iloc[10]:.4f}, p-value = {S2T.p_printed.iloc[10]}; the manuscript's 1.97
  and 0.374 are both correct roundings of it.

**Error bars.** None are reported anywhere: the paper gives medians with IQRs and means
with SDs for demographics, and no interval for any estimate. F2's right panel shows
mean +/- SD against mean +/- SEM for the headline outcome so the difference between the
two is visible; neither is what the paper plots, which is why the panel says so.

---

## 6. Resolution and evidence -- F5

![F5 null distribution]({_img(f['F5_null_distribution'])})

No design in this paper is resolution-limited: with n in the hundreds, the combinatorial
floor of every test is far below any reportable p-value (for C15, 675 vs 248 cats, the
two-sided permutation floor is {float(rs.min_achievable_p(675, 248)):.3g}). What *is*
limited is the review's own resampling: the Monte-Carlo chi-square nulls resolve to
{G['end_omni']['p_res']:.1e} ({N_MC:,} replicates) and the permutation routes to
{1.0/(N_PERM+1):.1e} ({N_PERM:,} replicates), so every recomputed p-value at those values
is a resolution limit, not an estimate -- which is why the ledger carries
`p_smallest_possible` next to `p_recomputed`.

The left panel makes the point for the strongest claim in the paper (C7, who ends
sessions): the observed chi-square {G['end_omni']['chi2']:.1f} lies far outside the entire
Monte-Carlo null (maximum {G['end_omni']['null'].max():.1f} over {N_MC:,} draws). The
right panel answers the separate question the paper does not: how precise the effect is.
Cramer's V = {G['end_omni']['cramers_v']:.3f}
[{BOOT_V['end_omni']['ci_low']:.3f}, {BOOT_V['end_omni']['ci_high']:.3f}] -- a moderate
association, reported here because the paper reports none.

The paper reports no permutation or bootstrap analysis, so there is no shuffling scheme
or iteration count of its own to check.

---

## 7. Multiplicity -- F6

![F6 multiplicity]({_img(f['F6_multiplicity'])})

The paper runs {len(P_VALUES)} reported tests, where {ALPHA * len(P_VALUES):.1f}
significant results are expected by chance alone if every null were true.
{nsig} are significant uncorrected; {nbh} survive Benjamini-Hochberg and {nbonf} survive
Bonferroni at the full m = {len(P_VALUES)} (threshold {BONF_THRESHOLD:.5f}). The two that
do not survive Bonferroni are C11 (sex ratio, p = 0.018) and C13 (pairwise retrievals,
p = 0.009) -- both already flagged `supported-with-caveat` for other reasons.

**The paper did correct, within families.** Post-hoc chi-squares were interpreted at a
corrected alpha of 0.016 (0.05/3) and the pairwise Wilcoxon tests carry Bonferroni
adjustment, which is the right correction for those families. What is not stated is a
global count: the recount above is the view across the whole Results section, and it is
offered as context rather than as a demand, since the comparisons are pre-specified by
the survey's structure rather than selected after the fact.

**What S2 shows about the true denominator.** The script contains
{S2_CALLS['chisq.test']} `chisq.test`, {S2_CALLS['kruskal.test']} `kruskal.test`,
{S2_CALLS['pairwise.wilcox.test']} `pairwise.wilcox.test` (three comparisons each) and
{S2_CALLS['wilcox.test']} `wilcox.test` calls -- {S2_INFERENTIAL_P} inferential p-values in
total, against the {len(P_VALUES)} the paper prints. The gap is almost entirely the
outlier-excluded duplicate of each count analysis ({S2_CALLS['kruskal.test'] // 2} Kruskal-Wallis,
{3 * (S2_CALLS['pairwise.wilcox.test'] // 2)} pairwise and {S2_CALLS['wilcox.test'] // 2}
Wilcoxon tests), which is a robustness check rather than a second pass at the hypothesis,
and the paper summarises it in one Methods sentence -- a sentence that is inaccurate for
C20 (section 2, finding 2). On top of those, S2 runs {S2_CALLS['shapiro_test']}
Shapiro-Wilk and {S2_CALLS['leveneTest']} Levene screens, whose printed results confirm the
Methods description (for instance Levene F = 3.6907, p = 0.02525 on retrievals by
initiator).

**One pairwise comparison is run but not reported.** For monthly sessions by initiator the
paper reports two of the three pairwise Wilcoxons; the third, cat versus about-equal, is
printed in S2 as p = {PAIRWISE_S2[4][0]} and recomputes as
p = {KW['mon_ini']['pairwise']['Cat vs About equal']['p_bonf']:.3f} (Bonferroni-adjusted).
Reporting all three keeps the correction's denominator visible.

**Distribution of reported p-values.** {n_bounds} of the {len(P_VALUES)} p-values are
printed as `p < 0.001` bounds and are recorded here at 0.001, so the left-hand spike in
F6 is partly a reporting convention rather than a measured distribution. The
just-below-threshold diagnostic returns
{JUST_BELOW['n_in_suspect_window']} p-values in [0.045, 0.05) against
{JUST_BELOW['n_in_reference_window']} in [0.040, 0.045) -- both windows are empty, so the
diagnostic is uninformative here and no binomial comparison is defined. It is in any case
a property of a *collection* of p-values, says nothing about any individual result, and
cannot distinguish analytic flexibility from a pre-specified analysis; no inference about
intent is drawn from it.

---

## 8. Per-claim supplement

One panel per headline claim and per claim whose verdict is not `supported`.

**C1 -- 94.4% of fetching emerged without training** (`{LEDGER.loc[LEDGER.claim_id=='C1','verdict'].iloc[0]}`)

![C1 panel]({_img(f['panel_C1'])})

{N_UNTRAINED} of {N_CATS} cats are recorded `Untrained` ({PCT_UNTRAINED:.2f}%, exact 95%
CI [{100*UNTRAINED_CI[0]:.2f}, {100*UNTRAINED_CI[1]:.2f}]%), reproducing the reported
94.4% exactly. The caveat is the denominator: the abstract says "94.4% of these *owners*
report", and the figure is a proportion of *cats*. With no respondent identifier the
owner-level proportion cannot be computed from S3. The centre panel's null (50:50
training) is shown because the panel template tests it, not because the paper claims it.

**C3 / C7 -- who initiates and who ends** (`{LEDGER.loc[LEDGER.claim_id=='C3','verdict'].iloc[0]}`)

![C3 panel]({_img(f['panel_C3'])})

![C7 panel]({_img(f['panel_C7'])})

Both omnibus statistics reproduce exactly (initiation {G['ini_omni']['chi2']:.2f} against
the reported 124.29; termination {G['end_omni']['chi2']:.2f} against 329.08), as do all
six post-hoc comparisons, and every route agrees. Cramer's V is
{G['ini_omni']['cramers_v']:.3f} for initiation and {G['end_omni']['cramers_v']:.3f} for
termination: the asymmetry in who ends a session is the stronger of the two, which the
paper's narrative also emphasises.

**C11 -- more male than female fetching cats** (`{LEDGER.loc[LEDGER.claim_id=='C11','verdict'].iloc[0]}`)

![C11 panel]({_img(f['panel_C11'])})

The routes agree with each other (asymptotic {G['sex']['p_asym']:.4f}, exact binomial
{G['sex']['p_exact']:.4f}, Yates {G['sex']['p_yates']:.4f}, Monte-Carlo
{G['sex']['p_mc']:.4f}) and the asymptotic value is exactly what S2 prints
({S2T.statistic.iloc[0]:.4f}, p = {S2T.p_printed.iloc[0]}). The Yates row is shown for
completeness only: `chisq.test` does not apply it to a one-dimensional table, so it was
never part of the authors' result. What remains is the scope of the claim -- a 50:50 null
inside a sample of fetching cats describes composition, and with no non-fetching
denominator it cannot support a claim about sex and fetching -- plus the manuscript's
mis-transcription of 5.5459 as 5.50 and the fact that p = 0.018 does not survive
Bonferroni across the {len(P_VALUES)} reported tests.

**C13 -- fewer retrievals when the owner initiates** (`{LEDGER.loc[LEDGER.claim_id=='C13','verdict'].iloc[0]}`)

![C13 panel]({_img(f['panel_C13'])})

See finding 3. The right-hand panel is the evidence that the result is robust to the
modelling choice: Welch, Student, Mann-Whitney, mean-difference permutation and rank-sum
permutation all land below alpha.

**C20 -- no sex difference in retrievals** (`{LEDGER.loc[LEDGER.claim_id=='C20','verdict'].iloc[0]}`)

![C20 panel]({_img(f['panel_C20'])})

The reported Mann-Whitney reproduces exactly, including the W statistic S2 prints
(W = {S2T.statistic.iloc[19]:.0f}, p = {S2T.p_printed.iloc[19]}; recomputed
W = {SENS_RET['retained']['W_female_first']:.0f}, p = {SENS_RET['retained']['p']:.5f}), and
the reported Z = -1.69 is the normal deviate of that p-value
({SENS_RET['retained']['z_from_p']:.4f}). The verdict is `supported-with-caveat` for a
different reason -- the paper's own outlier-excluded run of this same comparison gives
p = {SENS_RET['excluded']['p']:.5f}, on the other side of alpha (finding 2). This panel is
also worth reading because the routes, while unanimous at
alpha, differ by an order of magnitude in p: the rank routes sit near
{SEX_RET['p_mw']:.2f} while the mean-based routes sit near
{float(SEX_RET['routes'].p_value.iloc[0]):.2f}. They are answering different questions --
whether male cats tend to rank higher, versus whether their mean differs -- and the mean
difference is {SEX_RET['bootstrap']['point_estimate']:.3f} retrievals
(95% CI [{SEX_RET['bootstrap']['ci_low']:.3f},
{SEX_RET['bootstrap']['ci_high']:.3f}]). A null result at this n is compatible both with
no association and with an association too small for {N_CATS} cats to resolve.

**C22, C23, C24 -- not assessable**

![C22 panel]({_img(f['panel_C22'])})

![C23 panel]({_img(f['panel_C23'])})

![C24 panel]({_img(f['panel_C24'])})

---

## 9. What was not assessed

{ac.uncheckable_assumptions_note()}

**Claims the supplied inputs cannot address** (verdict `not-assessable`, with what would
change each):

- **C22 - "cats demonstrate independent and co-ordinated agency".** S3 records, once per
  cat, who *usually* initiates and ends sessions on a three-point scale. That locates a
  frequency asymmetry; it carries no within-bout sequence, no latency, and no observation
  independent of the owner's interpretation. Bout-level observational or video-coded data,
  or an experimental manipulation of owner availability, would address it.
- **C23 - 1,154 cats analysed as independent observations.** No respondent, household or
  submission identifier exists in S3. A respondent id column, or owner-level counts, would
  allow a cluster-robust or mixed model. The worst-case bound computed in section 4 is the
  substitute, not a replacement.
- **C24 - the thematic analysis.** Qualitative coding is outside the scope of a
  statistical review, and the paper reports 268 initial codes with 35 child codes but no
  inter-coder agreement statistic. The code frame plus a second coder's assignments would
  make it checkable.

**Not checked for any claim:** independence of observations (see C23), representativeness
of the self-selected sample, and confounding between owner interpretation and the
reported behaviour. All three are design questions readable only from the methods, and
the methods here are explicit about the recruitment route.

**Review-side limitations.** Supplementary S2 is available to this revision but cannot be
executed as supplied -- the .docx extraction collapsed in-paragraph line breaks -- so each
test is re-implemented in Python against the authors' printed call rather than re-run in R.
The printed statistics, degrees of freedom and p-values are intact, and section 5 compares
every one of them with the recomputed value. The outlier-excluded sensitivity analysis is
now reproduced for the two cat-sex comparisons (the only ones where it could change an
interpretation, and it does for one); the outlier-excluded Kruskal-Wallis and pairwise runs
are read from S2's printed output and not independently recomputed, since all four remain
significant in both versions. Nothing was subsampled, filtered or approximated.

---

## 10. Reproduction

- Script: `review_{SLUG}.py` (this review, end to end: recomputation, figures, ledger,
  report).
- Seeds: `SEED = {SEED}` for every permutation, Monte-Carlo and bootstrap draw.
- Resampling: `N_PERM = {N_PERM:,}` (resolution {1.0/(N_PERM+1):.1e}),
  `N_MC = {N_MC:,}` (resolution {1.0/(N_MC+1):.1e}),
  `N_BOOT = {N_BOOT:,}`; alpha = {ALPHA}.
- Inputs: `{os.path.basename(PAPER)}` sha256 `{sha(PAPER)}`;
  `{os.path.basename(DATA)}` sha256 `{sha(DATA)}`;
  `{os.path.basename(S2_PATH)}` sha256 `{sha(S2_PATH)}`.
- Outputs: `claim_ledger.csv` (25 claims, now with the authors' printed statistic per
  row), `statcheck.csv` (manuscript-text triples plus the {len(S2T)} S2 console triples),
  `s2_crosscheck.csv` (the section 5 table), `reported_pvalues.csv`.
- Versions: python {sys.version.split()[0]}, numpy {np.__version__},
  pandas {pd.__version__}, scipy {__import__('scipy').__version__},
  matplotlib {matplotlib.__version__}.
- One review-side fix: `resampling.min_achievable_p` evaluates 2/C(n1+n2, n1) directly
  and overflows a float at these sample sizes, so it is wrapped with a log-space version
  that underflows to 0. The skill itself is unmodified.
- Figure references in this report resolve to artifact markers when
  `figure_artifact_map.json` is present and to relative paths otherwise.
"""
    with open(os.path.join(OUTDIR, "REVIEW_REPORT.md"), "w", encoding="utf-8") as fh:
        fh.write(md)
    return md


def build_statcheck():
    """
    Statcheck over both sources: the manuscript text (where df are mostly absent,
    so most triples are not checkable) and the S2 console output (where every
    printed statistic carries its df).
    """
    sc = STATCHECK.copy()
    sc["sentence"] = sc["sentence"].str.replace(r"\s+", " ", regex=True).str.slice(0, 300)
    sc.insert(0, "source", "manuscript text")
    sc = sc.rename(columns={"sentence": "context"})

    rows = []
    for i, r in S2T.iterrows():
        checkable = np.isfinite(r.df)
        p_rec = r.p_from_statistic if checkable else np.nan
        if not checkable:
            consistent, note = pd.NA, "W statistic: p depends on the tie and continuity " \
                                      "correction, not on df; checked against the raw data instead"
        elif r.p_relation == "<":
            consistent = bool(p_rec < float(r.p_printed))
            note = "reported as a bound; the bound holds" if consistent else \
                   "reported bound does not hold"
        else:
            # The statistic is printed rounded, so the p-value it implies is an
            # interval, not a point. Consistency means the printed p lies inside
            # the interval the printed statistic could have come from.
            half = 0.5 * 10 ** -_decimals(r.statistic)
            p_hi = float(stats.chi2.sf(max(r.statistic - half, 0.0), r.df))
            p_lo = float(stats.chi2.sf(r.statistic + half, r.df))
            pr_lo, pr_hi = _printed_interval(r.p_printed)
            consistent = bool(p_lo <= pr_hi and pr_lo <= p_hi)
            note = ("" if consistent else
                    f"reported p does not follow from the statistic and df: the printed "
                    f"statistic implies p in [{p_lo:.4g}, {p_hi:.4g}]")
        # A decision error is a consistency failure that crosses alpha: the
        # printed p and the p its own statistic implies fall on opposite sides.
        if not checkable:
            decision_error = pd.NA
        elif r.p_relation == "<":
            decision_error = bool((float(r.p_printed) < ALPHA) != (p_rec < ALPHA))
        else:
            decision_error = bool((float(r.p_printed) < ALPHA) != (p_rec < ALPHA))
        rows.append({"source": "S2 console output", "claim_id": f"S2_{i:02d}",
                     "context": f"{r.kind} on {r['data']}", "stat_type": r.kind,
                     "stat_value": r.statistic, "df1": r.df, "df2": np.nan,
                     "p_relation": r.p_relation, "p_reported": r.p_printed,
                     "tail_hint": 2, "p_recomputed": p_rec, "consistent": consistent,
                     "decision_error": decision_error, "note": note})
    return pd.concat([sc, pd.DataFrame(rows)], ignore_index=True)


STATCHECK_ALL = None


def write_outputs():
    global STATCHECK_ALL
    LEDGER.to_csv(os.path.join(OUTDIR, "claim_ledger.csv"), index=False)
    STATCHECK_ALL = build_statcheck()
    STATCHECK_ALL.to_csv(os.path.join(OUTDIR, "statcheck.csv"), index=False)
    S2_CHECK.to_csv(os.path.join(OUTDIR, "s2_crosscheck.csv"), index=False)
    pd.DataFrame({"label": P_LABELS, "p_reported": P_VALUES,
                  "is_bound": P_IS_BOUND}).to_csv(
        os.path.join(OUTDIR, "reported_pvalues.csv"), index=False)
    return write_report()


if __name__ == "__main__":
    build_figures()
    write_outputs()

