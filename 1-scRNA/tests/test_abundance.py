"""Differential-abundance tests, including a regression test for the design matrix."""

from __future__ import annotations

import numpy as np
import pandas as pd

from icb_scrna import abundance as ab


def _toy(n_per_group=8, seed=0):
    """Two groups where cell type A really is enriched in the test group."""
    rng = np.random.default_rng(seed)
    rows, meta = [], []
    for g, pA in [("test", 0.40), ("ref", 0.15)]:
        for i in range(n_per_group):
            total = 200
            nA = rng.binomial(total, pA)
            rows.append({"A": nA, "B": total - nA})
            meta.append({"sample": f"{g}{i}", "grp": g, "patient": f"pat_{g}{i}"})
    counts = pd.DataFrame(rows, index=[m["sample"] for m in meta])
    smeta = pd.DataFrame(meta).set_index("sample")
    return counts, smeta


def test_design_matrix_keeps_the_intercept():
    _, smeta = _toy()
    X, _ = ab._design(smeta, "grp", "test", "ref", covariates=None)
    assert "Intercept" in X.columns
    assert (X["Intercept"] == 1.0).all()


def test_design_matrix_keeps_intercept_with_constant_covariate():
    _, smeta = _toy()
    smeta = smeta.assign(therapy="anti-PD1")  # no variation at all
    X, _ = ab._design(smeta, "grp", "test", "ref", covariates=["therapy"])
    assert "Intercept" in X.columns
    assert list(X.columns) == ["Intercept", "grp[test]"]


def test_known_enrichment_is_recovered_with_correct_sign():
    counts, smeta = _toy()
    res = ab.test_abundance(counts, smeta, "grp", "test", "ref").set_index("cell_type")
    assert res.loc["A", "log_odds_ratio"] > 0
    assert res.loc["B", "log_odds_ratio"] < 0
    assert res.loc["A", "p_adj_BH"] < 0.05
    # the fitted odds ratio must match the empirical one
    emp = (0.40 / 0.60) / (0.15 / 0.85)
    assert abs(res.loc["A", "odds_ratio"] - emp) < 0.5 * emp


def test_null_data_is_not_significant():
    rng = np.random.default_rng(1)
    rows, meta = [], []
    for g in ("test", "ref"):
        for i in range(8):
            nA = rng.binomial(200, 0.3)
            rows.append({"A": nA, "B": 200 - nA})
            meta.append({"sample": f"{g}{i}", "grp": g, "patient": f"pat_{g}{i}"})
    counts = pd.DataFrame(rows, index=[m["sample"] for m in meta])
    smeta = pd.DataFrame(meta).set_index("sample")
    res = ab.test_abundance(counts, smeta, "grp", "test", "ref")
    assert (res["p_adj_BH"] > 0.05).all()


def test_sign_of_fit_agrees_with_difference_in_mean_proportions():
    counts, smeta = _toy(seed=3)
    res = ab.test_abundance(counts, smeta, "grp", "test", "ref")
    diff = res["mean_prop_test"] - res["mean_prop_ref"]
    assert (np.sign(diff) == np.sign(res["log_odds_ratio"])).all()


def test_proportion_table_excludes_enriched_sorts_and_small_biopsies():
    obs = pd.DataFrame({
        "sample": ["s1"] * 30 + ["s2"] * 30 + ["s3"] * 5,
        "cell_type": ["A"] * 15 + ["B"] * 15 + ["A"] * 30 + ["A"] * 5,
        "enrichment": ["CD45"] * 30 + ["T_enriched"] * 30 + ["CD45"] * 5,
    })
    counts, props = ab.proportion_table(obs, min_cells=20)
    assert list(counts.index) == ["s1"]           # s2 enriched, s3 too small
    np.testing.assert_allclose(props.loc["s1"].sum(), 1.0)


def test_paired_change_collapses_repeat_biopsies_to_one_pair_per_patient():
    idx = ["a_pre", "a_post1", "a_post2", "b_pre", "b_post",
           "c_pre", "c_post", "d_pre", "d_post", "e_pre", "e_post"]
    props = pd.DataFrame({"A": [0.1, 0.5, 0.5, 0.1, 0.5, 0.2, 0.6, 0.1, 0.4, 0.2, 0.5],
                          "B": [0.9, 0.5, 0.5, 0.9, 0.5, 0.8, 0.4, 0.9, 0.6, 0.8, 0.5]},
                         index=idx)
    smeta = pd.DataFrame({
        "patient": ["a", "a", "a", "b", "b", "c", "c", "d", "d", "e", "e"],
        "timepoint": ["Pre", "Post", "Post", "Pre", "Post", "Pre", "Post",
                      "Pre", "Post", "Pre", "Post"]}, index=idx)
    res = ab.paired_change(props, smeta).set_index("cell_type")
    assert (res["n_patients_paired"] == 5).all()
    assert res.loc["A", "median_change"] > 0
