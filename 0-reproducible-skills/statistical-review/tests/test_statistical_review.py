"""
Tests for the statistical-review utilities.

Several of these pin the specific numerical facts from the lecture, so the
module cannot silently drift away from them:

  * 3 vs 3 permutation has a two-sided floor of exactly 0.10
  * a sign test on 5 pairs has a smallest two-sided p of exactly 0.0625
  * a dataset can have identical means (t-test p = 1.0) and still differ by a
    count-based statistic
  * direction and magnitude can disagree: 9 of 10 subjects improve while the
    mean change is negative

Run from the skill root:   python -m pytest tests/ -q
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utilities import (
    assumption_checks as ac,
    claim_extraction as ce,
    multiplicity as mult,
    recompute as rc,
    resampling as rs,
    review_figures as rf,
)


# --------------------------------------------------------------------------
# claim extraction
# --------------------------------------------------------------------------

SAMPLE_TEXT = (
    "Tumour volume was reduced in treated animals (t(18) = 2.31, p = 0.033). "
    "Expression differed across the three groups (F(2, 27) = 4.10, p = 0.028). "
    "Marker abundance correlated with survival (r(48) = 0.31, p = 0.029). "
    "Counts differed between conditions (chi2(1) = 6.20, p = 0.013). "
    "Data are shown as mean +/- SEM from n = 3 independent experiments. "
    "The hazard ratio was 0.62 (95% CI 0.41 to 0.94). "
    "P-values were corrected for multiple comparisons using the Benjamini-Hochberg procedure. "
    "Groups were compared with Welch's t-test; paired samples used a Wilcoxon signed-rank test."
)


def test_find_reported_tests_parses_all_four_statistics():
    tests = ce.find_reported_tests(SAMPLE_TEXT, page=1)
    assert set(tests["stat_type"]) >= {"t", "F", "r", "chi2"}
    t_row = tests[tests["stat_type"] == "t"].iloc[0]
    assert t_row["stat_value"] == pytest.approx(2.31)
    assert t_row["df1"] == 18
    assert t_row["p_reported"] == pytest.approx(0.033)
    f_row = tests[tests["stat_type"] == "F"].iloc[0]
    assert (f_row["df1"], f_row["df2"]) == (2, 27)


def test_find_sample_sizes_captures_the_unit():
    sizes = ce.find_sample_sizes(SAMPLE_TEXT, page=1)
    assert 3 in list(sizes["n"])
    row = sizes[sizes["n"] == 3].iloc[0]
    assert "independent experiment" in (row["replication_language"] or [])


def test_find_intervals_and_named_tests():
    intervals = ce.find_intervals(SAMPLE_TEXT, page=1)
    ci = intervals[intervals["kind"] == "interval"].iloc[0]
    assert (ci["lo"], ci["hi"], ci["level"]) == (0.41, 0.94, 95)
    named = set(ce.find_named_tests(SAMPLE_TEXT, page=1)["test_named"])
    assert {"welch_t", "wilcoxon_signed_rank"} <= named


def test_build_claim_inventory_and_text_quality():
    inv = ce.build_claim_inventory([SAMPLE_TEXT, SAMPLE_TEXT])
    assert set(inv) >= {"tests", "p_values", "sample_sizes", "intervals", "named_tests",
                        "error_bars", "multiplicity", "replication", "text_quality"}
    assert len(inv["tests"]) == 8  # four statistics on each of two pages
    assert inv["tests"]["claim_id"].is_unique
    assert bool(inv["text_quality"]["usable_text_layer"].iloc[0])


def test_empty_pdf_is_reported_as_unusable_not_as_no_claims():
    quality = ce.pdf_text_quality(["", "   ", ""])
    assert quality["usable_text_layer"] is False
    assert quality["empty_pages"] == 3


# --------------------------------------------------------------------------
# recomputation
# --------------------------------------------------------------------------

def test_recompute_p_matches_scipy_for_each_statistic():
    from scipy import stats
    assert rc.recompute_p("t", 2.31, 18) == pytest.approx(2 * stats.t.sf(2.31, 18))
    assert rc.recompute_p("F", 4.10, 2, 27) == pytest.approx(stats.f.sf(4.10, 2, 27))
    assert rc.recompute_p("chi2", 6.20, 1) == pytest.approx(stats.chi2.sf(6.20, 1))
    assert rc.recompute_p("z", 1.96) == pytest.approx(2 * stats.norm.sf(1.96), rel=1e-6)
    assert np.isnan(rc.recompute_p("t", 2.31, np.nan))


def test_statcheck_flags_a_decision_error_but_not_a_rounding_difference():
    tests = pd.DataFrame([
        # consistent: t(18)=2.31 really does give p ~ 0.033
        {"stat_type": "t", "stat_value": 2.31, "df1": 18, "df2": np.nan,
         "p_relation": "=", "p_reported": 0.033, "tail_hint": 2},
        # decision error: the true p is ~0.33, reported as 0.033
        {"stat_type": "t", "stat_value": 1.00, "df1": 18, "df2": np.nan,
         "p_relation": "=", "p_reported": 0.033, "tail_hint": 2},
        # not checkable: no df
        {"stat_type": "t", "stat_value": 2.31, "df1": np.nan, "df2": np.nan,
         "p_relation": "=", "p_reported": 0.033, "tail_hint": 2},
    ])
    out = rc.statcheck_table(tests)
    assert out.loc[0, "consistent"] is True or bool(out.loc[0, "consistent"])
    assert bool(out.loc[0, "decision_error"]) is False
    assert bool(out.loc[1, "consistent"]) is False
    assert bool(out.loc[1, "decision_error"]) is True
    assert pd.isna(out.loc[2, "consistent"])
    assert "not checkable" in out.loc[2, "note"]


def test_statcheck_handles_bounded_p_values():
    tests = pd.DataFrame([
        {"stat_type": "t", "stat_value": 8.0, "df1": 18, "df2": np.nan,
         "p_relation": "<", "p_reported": 0.001, "tail_hint": 2},
        {"stat_type": "t", "stat_value": 1.2, "df1": 18, "df2": np.nan,
         "p_relation": "<", "p_reported": 0.001, "tail_hint": 2},
    ])
    out = rc.statcheck_table(tests)
    assert bool(out.loc[0, "consistent"]) is True
    assert bool(out.loc[1, "consistent"]) is False


def test_grim_detects_an_impossible_mean():
    # with n = 7 integer observations, means are multiples of 1/7
    assert rc.grim_check(3.00, 7, decimals=2)["consistent"] is True
    assert rc.grim_check(3.15, 7, decimals=2)["consistent"] is False
    assert rc.grim_check(3.15, 500, decimals=2)["applicable"] is False


def test_error_bar_identity_separates_sd_from_sem():
    sd, n = 10.0, 25
    assert rc.error_bar_identity(10.0, sd, n)["best_match"] == "SD"
    assert rc.error_bar_identity(2.0, sd, n)["best_match"] == "SEM"
    assert rc.error_bar_identity(2.0, sd, n)["sd_to_sem_factor"] == pytest.approx(5.0)


def test_effect_from_summary_matches_a_direct_welch_test():
    from scipy import stats
    rng = np.random.default_rng(1)
    a = rng.normal(10, 2, 30); b = rng.normal(8, 3, 28)
    out = rc.effect_from_summary(a.mean(), a.std(ddof=1), a.size,
                                 b.mean(), b.std(ddof=1), b.size)
    direct = stats.ttest_ind(a, b, equal_var=False)
    assert out["welch_t"] == pytest.approx(direct.statistic, rel=1e-9)
    assert out["p_two_sided"] == pytest.approx(direct.pvalue, rel=1e-9)
    assert out["ci_low"] < out["mean_difference"] < out["ci_high"]


def test_reproduction_verdict_categories():
    assert rc.reproduction_verdict(0.033, 0.0331) == "reproduced"
    assert rc.reproduction_verdict(0.033, 0.33) == "decision-changing"
    assert rc.reproduction_verdict(0.001, 0.01) == "minor-discrepancy"
    assert rc.reproduction_verdict(np.nan, 0.01) == "not-checkable"


# --------------------------------------------------------------------------
# assumptions and design
# --------------------------------------------------------------------------

def test_distribution_shape_never_returns_a_verdict():
    rng = np.random.default_rng(2)
    shape = ac.distribution_shape(rng.normal(size=12), name="g")
    assert "interpretation_note" in shape
    assert not any(k in shape for k in ("passes", "normal", "assumptions_met"))
    assert shape["n_below_small_threshold"] is True
    assert "little power" in shape["interpretation_note"]


def test_variance_comparison_recommends_welch_on_unequal_spread():
    rng = np.random.default_rng(3)
    out = ac.variance_comparison(rng.normal(0, 1, 40), rng.normal(0, 4, 40))
    assert out["variance_ratio"] > ac.VARIANCE_RATIO_FLAG
    assert out["recommended_t_variant"] == "welch_t"
    assert "not appropriate" in out["note"]


def test_paired_check_flags_direction_versus_magnitude():
    # the lecture's example: 9 of 10 improve slightly, one worsens a lot
    before = np.zeros(10)
    after = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, -30], float)
    out = ac.paired_difference_check(after, before)
    assert out["n_positive"] == 9
    assert out["median_difference"] > 0 > out["mean_difference"]
    assert out["direction_vs_magnitude_note"]
    assert out["symmetry_concern"] is True


def test_design_independence_detects_pseudoreplication():
    rng = np.random.default_rng(4)
    rows = []
    for unit in range(6):
        offset = rng.normal(0, 5)           # large between-unit variance
        grp = "treated" if unit < 3 else "control"
        for _ in range(50):                 # 50 cells per animal
            rows.append({"mouse": f"m{unit}", "group": grp,
                         "value": offset + rng.normal(0, 0.5)})
    df = pd.DataFrame(rows)
    out = ac.design_independence(df, "mouse", "group", "value")
    assert out["n_observations"] == 300 and out["n_units"] == 6
    assert out["icc"] > 0.5
    assert out["effective_n"] < 20
    assert out["pseudoreplication_risk"] is True
    assert out["units_per_group"] == {"control": 3, "treated": 3}


def test_small_sample_notes_report_the_resolution_floor():
    out = ac.small_sample_notes(3, 3)
    assert out["resolution_floors"]["permutation_two_sample_min_p"] == pytest.approx(0.10)
    assert any("cannot produce p < 0.05" in n for n in out["notes"])


def test_uncheckable_assumptions_note_lists_design_questions():
    note = ac.uncheckable_assumptions_note()
    assert "independence of observations" in note
    assert "confounding" in note


# --------------------------------------------------------------------------
# resampling and resolution
# --------------------------------------------------------------------------

def test_permutation_floor_for_three_versus_three_is_exactly_one_tenth():
    assert rs.min_achievable_p(3, 3, two_sided=True) == pytest.approx(0.10)
    assert rs.min_achievable_p(3, 3, two_sided=False) == pytest.approx(0.05)


def test_floors_survive_survey_sized_samples():
    # Regression: C(1154, 577) has 345 digits, so a direct 1.0/total raised
    # OverflowError and took down any review of a survey-sized dataset.
    # Found by the Forman et al. (2023) review, n = 1154.
    assert rs.min_achievable_p(577, 577) == 0.0
    assert rs.min_achievable_p(617, 537) == 0.0
    assert rs.min_sign_test_p(1154) == 0.0
    assert rs.min_achievable_p(1154, paired=True) == 0.0
    # and the small-n answers are unchanged by the log-space rewrite
    assert rs.min_achievable_p(3, 3) == pytest.approx(0.10)
    assert rs.min_achievable_p(4, 4) == pytest.approx(2 / 70)
    assert rs.min_achievable_p(10, paired=True) == pytest.approx(2 / 1024)
    assert rs.min_sign_test_p(5) == pytest.approx(0.0625)


def test_sign_test_floor_for_five_pairs_is_exactly_0_0625():
    assert rs.min_sign_test_p(5) == pytest.approx(0.0625)
    assert rs.min_sign_test_p(5) > 0.05


def test_permutation_reports_a_binding_floor_even_for_a_huge_effect():
    out = rs.permutation_test([100.0, 101.0, 102.0], [0.0, 1.0, 2.0], n_perm=2000)
    assert out["min_achievable_p"] == pytest.approx(0.10)
    assert out["floor_binds"] is True
    assert out["p_value"] >= 0.05
    assert "cannot reach p < 0.05" in out["note"]


def test_p_value_resolution_formula():
    assert rs.p_value_resolution(999) == pytest.approx(0.001)
    out = rs.permutation_test(np.arange(20.0), np.arange(20.0) + 3, n_perm=500)
    assert out["p_value"] >= 1 / 501
    assert out["p_resolution"] == pytest.approx(1 / 501)


def test_identical_means_can_still_differ_by_a_count_statistic():
    # The lecture's example: means match exactly, so the t-test sees nothing,
    # but the groups differ in how the observations are arranged.
    a = np.array([-10.0, 0.0, 10.0, -10.0, 0.0, 10.0, -10.0, 0.0, 10.0, 0.0])
    b = np.array([-1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    b = b - b.mean() + a.mean()
    from scipy import stats
    assert stats.ttest_ind(a, b, equal_var=False).pvalue > 0.9
    var_out = rs.permutation_test(a, b, statistic="var_ratio", n_perm=5000, random_state=7)
    assert var_out["p_value"] < 0.05
    assert var_out["statistic_name"] == "var_ratio"


def test_paired_permutation_shuffles_within_pairs():
    rng = np.random.default_rng(5)
    before = rng.normal(10, 2, 20)
    after = before + rng.normal(1.5, 0.5, 20)
    out = rs.permutation_test(after, before, paired=True, n_perm=4000, random_state=0)
    assert out["paired"] is True
    assert out["n1"] == 20 and out["n2"] is None
    assert out["p_value"] < 0.05
    assert out["min_achievable_p"] == pytest.approx(2 / 2 ** 20)


def test_stratified_permutation_respects_strata():
    rng = np.random.default_rng(6)
    # Strong stratum effect, no treatment effect within stratum: unstratified
    # permutation can be fooled by imbalance, stratified permutation is not.
    x = np.concatenate([rng.normal(0, 1, 15), rng.normal(10, 1, 5)])
    y = np.concatenate([rng.normal(0, 1, 5), rng.normal(10, 1, 15)])
    strata = np.array(["A"] * 15 + ["B"] * 5 + ["A"] * 5 + ["B"] * 15)
    plain = rs.permutation_test(x, y, n_perm=2000, random_state=0)
    strat = rs.permutation_test(x, y, n_perm=2000, strata=strata, random_state=0)
    assert strat["stratified"] is True
    assert strat["p_value"] > plain["p_value"]


def test_sign_test_counts_directions_and_drops_zeros():
    out = rs.sign_test([1, 2, 3, 4, -1, 0, 0])
    assert out["n_pairs"] == 5 and out["n_zero_dropped"] == 2
    assert out["n_positive"] == 4
    assert out["min_achievable_p"] == pytest.approx(0.0625)
    assert out["floor_binds"] is True


def test_bootstrap_ci_brackets_the_estimate_and_reports_bias():
    rng = np.random.default_rng(8)
    x = rng.normal(5, 1, 60)
    out = rs.bootstrap_ci(x, np.mean, n_boot=3000, random_state=0)
    assert out["ci_low"] < out["point_estimate"] < out["ci_high"]
    assert abs(out["bootstrap_bias"]) < 0.1
    assert out["method"] == "BCa"
    small = rs.bootstrap_ci(x, np.mean, n_boot=500, random_state=0)
    assert "raise n_boot" in small["note"]


def test_three_routes_returns_one_row_per_route_with_assumptions():
    rng = np.random.default_rng(9)
    out = rs.three_routes(rng.normal(1, 1, 25), rng.normal(0, 1, 25), n_perm=2000)
    assert len(out) == 4
    assert {"route", "test", "p_value", "assumes", "min_achievable_p"} <= set(out.columns)
    assert out["assumes"].str.len().min() > 10
    paired = rs.three_routes(rng.normal(1, 1, 12), rng.normal(0, 1, 12), paired=True,
                             n_perm=2000)
    assert set(paired["test"]) == {"paired t-test", "Wilcoxon signed-rank", "sign test",
                                   "sign-flip permutation (mean difference)"}


# --------------------------------------------------------------------------
# multiplicity
# --------------------------------------------------------------------------

def test_bh_and_bonferroni_against_statsmodels():
    from statsmodels.stats.multitest import multipletests
    rng = np.random.default_rng(10)
    p = np.concatenate([rng.uniform(0, 0.01, 10), rng.uniform(0, 1, 90)])
    bh_ref = multipletests(p, alpha=0.05, method="fdr_bh")[1]
    bonf_ref = multipletests(p, alpha=0.05, method="bonferroni")[1]
    assert np.allclose(mult.benjamini_hochberg(p)["p_adjusted"].to_numpy(), bh_ref)
    assert np.allclose(mult.bonferroni(p)["p_adjusted"].to_numpy(), bonf_ref)


def test_discovery_recount_orders_corrections_and_states_expected_false_positives():
    rng = np.random.default_rng(11)
    p = np.concatenate([rng.uniform(0, 0.001, 8), rng.uniform(0, 1, 92)])
    out = mult.discovery_recount(p, alpha=0.05)
    n = dict(zip(out["summary"]["correction"], out["summary"]["n_significant"]))
    assert n["none"] >= n["Benjamini-Hochberg"] >= n["Bonferroni"]
    assert out["expected_false_positives_uncorrected"] == pytest.approx(5.0)
    assert out["m_tests"] == 100


def test_just_below_rate_detects_a_planted_spike():
    rng = np.random.default_rng(12)
    clean = rng.uniform(0, 1, 200)
    spiked = np.concatenate([clean, rng.uniform(0.045, 0.05, 25)])
    assert mult.just_below_rate(spiked)["binomial_p_excess"] < 0.01
    assert mult.just_below_rate(clean)["binomial_p_excess"] > 0.05
    assert mult.just_below_rate([0.01, 0.2])["interpretable"] is False


def test_max_statistic_null_controls_across_correlated_tests():
    rng = np.random.default_rng(13)
    null_matrix = rng.normal(size=(2000, 50))
    observed = np.array([4.5, 2.0, 0.1])
    out = mult.max_statistic_null(observed, null_matrix, alternative="two-sided")
    assert out["p_fwer"].is_monotonic_increasing
    assert out.loc[0, "p_fwer"] < 0.05 < out.loc[2, "p_fwer"]


def test_reported_vs_corrected_claims_names_what_survives():
    claims = pd.DataFrame({
        "claim_id": ["C1", "C2", "C3"],
        "p_reported": [0.0001, 0.03, 0.6],
    })
    out = mult.reported_vs_corrected_claims(claims)
    assert list(out["survives"]) == ["uncorrected, BH, Bonferroni", "uncorrected, BH", "none"]


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

@pytest.fixture
def demo_data():
    rng = np.random.default_rng(14)
    rows = []
    for unit in range(8):
        grp = "treated" if unit < 4 else "control"
        offset = rng.normal(2 if grp == "treated" else 0, 1.0)
        for _ in range(12):
            rows.append({"mouse": f"m{unit}", "group": grp,
                         "value": offset + rng.normal(0, 1.0)})
    return pd.DataFrame(rows)


def test_six_figures_render(tmp_path, demo_data):
    groups = {g: d["value"].to_numpy() for g, d in demo_data.groupby("group")}
    a, b = groups["treated"], groups["control"]

    f1 = rf.fig_design_map(demo_data, "mouse", "group", "value",
                           path=str(tmp_path / "f1.png"))
    f2 = rf.fig_raw_behind_claim(groups, "tumour volume (mm3)", test_name="Welch's t-test",
                                 p_value=0.03, reported_error="SEM",
                                 path=str(tmp_path / "f2.png"))
    f3 = rf.fig_assumption_qq(groups, path=str(tmp_path / "f3.png"))
    perm = rs.permutation_test(a, b, n_perm=1000, random_state=0)
    boot = rs.bootstrap_ci(a, np.mean, n_boot=1000, random_state=0)
    f5 = rf.fig_null_distribution(perm, boot, path=str(tmp_path / "f5.png"))
    rng = np.random.default_rng(15)
    f6 = rf.fig_multiplicity(np.concatenate([rng.uniform(0, 0.01, 5), rng.uniform(0, 1, 45)]),
                             path=str(tmp_path / "f6.png"))
    ledger = pd.DataFrame({
        "claim_id": ["C1", "C2"],
        "claim_short": ["treated tumours smaller", "marker correlates with survival"],
        "verdict": ["supported", "not-reproducible"],
        "effect_reported": [2.0, 0.31],
        "effect_recomputed": [1.8, 0.11],
        "ci_low": [0.9, -0.1], "ci_high": [2.7, 0.33],
        "p_reported": [0.03, 0.029], "p_recomputed": [0.04, 0.44],
    })
    f4 = rf.fig_claim_ledger(ledger, path=str(tmp_path / "f4.png"))

    for name in ["f1", "f2", "f3", "f4", "f5", "f6"]:
        assert (tmp_path / f"{name}.png").stat().st_size > 5000
    assert len(f1.axes) == 2 and len(f2.axes) == 2 and len(f3.axes) == 2
    assert len(f4.axes) == 2 and len(f5.axes) == 2 and len(f6.axes) == 2


def test_claim_panel_renders_with_three_routes(tmp_path, demo_data):
    groups = {g: d["value"].to_numpy() for g, d in demo_data.groupby("group")}
    a, b = groups["treated"], groups["control"]
    perm = rs.permutation_test(a, b, n_perm=1000, random_state=0)
    boot = rs.bootstrap_ci(a, np.mean, n_boot=1000, random_state=0)
    routes = rs.three_routes(a, b, n_perm=1000)
    fig = rf.fig_claim_panel("C1", groups, "tumour volume (mm3)", perm, boot, routes,
                             verdict="supported-with-caveat", path=str(tmp_path / "c1.png"))
    assert (tmp_path / "c1.png").stat().st_size > 5000
    assert len(fig.axes) == 3


def test_not_assessable_figure_states_the_reason(tmp_path):
    rf.not_assessable_figure("F6 - Multiplicity", "paper reports no p-values",
                             path=str(tmp_path / "na.png"))
    assert (tmp_path / "na.png").stat().st_size > 2000


def test_every_verdict_has_a_colour():
    assert set(rf.VERDICT_ORDER) == set(rf.VERDICT_COLORS)


def test_figures_find_the_data_visualization_skill():
    # The two skills live side by side in 0-reproducible-skills/, so the
    # palette helpers should resolve. If this fails the figures still render,
    # using the local fallback palette.
    dv = rf._dataviz()
    if dv is not None:
        assert len(rf._palette(3)) == 3
