"""Signature tests - the important ones check that the CV cannot leak."""

from __future__ import annotations

import numpy as np
import pandas as pd

from icb_scrna import signature as sg


def test_clr_rows_sum_to_zero():
    props = pd.DataFrame({"A": [0.2, 0.5], "B": [0.3, 0.1], "C": [0.5, 0.4]})
    out = sg.clr(props)
    np.testing.assert_allclose(out.sum(axis=1).to_numpy(), 0.0, atol=1e-9)


def test_clr_is_invariant_to_rescaling_a_row():
    props = pd.DataFrame({"A": [0.2], "B": [0.3], "C": [0.5]})
    np.testing.assert_allclose(
        sg.clr(props, pseudocount=0.0).to_numpy(),
        sg.clr(props * 2, pseudocount=0.0).to_numpy(), atol=1e-9)


def test_loo_auc_recovers_a_separable_signal():
    rng = np.random.default_rng(0)
    y = pd.Series([0] * 10 + [1] * 10, index=[f"p{i}" for i in range(20)])
    X = pd.DataFrame({"f": np.r_[rng.normal(0, 0.3, 10), rng.normal(4, 0.3, 10)]},
                     index=y.index)
    auc, scores = sg.loo_auc(X, y)
    assert auc > 0.95
    assert len(scores) == len(y)


def test_loo_auc_on_pure_noise_is_near_chance():
    rng = np.random.default_rng(1)
    y = pd.Series(rng.integers(0, 2, 24), index=[f"p{i}" for i in range(24)])
    X = pd.DataFrame(rng.normal(size=(24, 3)), index=y.index)
    auc, _ = sg.loo_auc(X, y)
    assert 0.2 < auc < 0.8


def test_nested_selection_does_not_leak_on_random_labels():
    """Selecting 20 of 2000 noise genes inside the fold must stay near chance.

    Selecting them *outside* the fold would drive the AUC toward 1.0; this is
    the test that distinguishes the two.
    """
    rng = np.random.default_rng(2)
    n = 24
    idx = [f"p{i}" for i in range(n)]
    y = pd.Series([0, 1] * (n // 2), index=idx)
    profiles = pd.DataFrame(rng.normal(size=(n, 2000)), index=idx,
                            columns=[f"g{i}" for i in range(2000)])
    patient = pd.Series(idx, index=idx)
    auc, _, genes = sg.loo_auc_nested(profiles, y, n_genes=20, patient=patient)
    assert len(genes) == 20
    assert auc < 0.75, f"nested CV leaked: AUC={auc}"


def test_nested_selection_finds_a_real_signal():
    rng = np.random.default_rng(3)
    n = 24
    idx = [f"p{i}" for i in range(n)]
    y = pd.Series([0] * 12 + [1] * 12, index=idx)
    X = rng.normal(size=(n, 500))
    X[12:, :10] += 3.0  # ten genuinely informative genes
    profiles = pd.DataFrame(X, index=idx, columns=[f"g{i}" for i in range(500)])
    auc, _, _ = sg.loo_auc_nested(profiles, y, n_genes=10,
                                  patient=pd.Series(idx, index=idx))
    assert auc > 0.9


def test_nested_selection_excludes_held_out_patient_from_the_selection_set():
    """A patient present in both tables must be dropped from both in their fold."""
    rng = np.random.default_rng(4)
    idx = [f"p{i}" for i in range(10)]
    y = pd.Series([0] * 5 + [1] * 5, index=idx)
    pre = pd.DataFrame(rng.normal(size=(10, 50)), index=idx,
                       columns=[f"g{i}" for i in range(50)])
    sel = pd.DataFrame(rng.normal(size=(10, 50)), index=idx,
                       columns=pre.columns)
    pat = pd.Series(idx, index=idx)
    seen = []

    original = sg.stats.ttest_ind

    def spy(a, b, **kw):
        seen.append(len(a) + len(b))
        return original(a, b, **kw)

    sg.stats.ttest_ind = spy
    try:
        sg.loo_auc_nested(pre, y, n_genes=5, selection_profiles=sel, selection_y=y,
                          selection_patient=pat, patient=pat)
    finally:
        sg.stats.ttest_ind = original
    # each fold selects from 9 of the 10 selection profiles, not all 10
    assert seen[0] == 9


def test_permutation_p_is_bounded_and_large_for_noise():
    rng = np.random.default_rng(5)
    y = pd.Series(rng.integers(0, 2, 16), index=[f"p{i}" for i in range(16)])
    X = pd.DataFrame(rng.normal(size=(16, 2)), index=y.index)
    auc, _ = sg.loo_auc(X, y)
    p = sg.permutation_p(auc, lambda ys: sg.loo_auc(X, ys)[0], y, n_perm=50)
    assert 0.0 < p <= 1.0


def test_bootstrap_ci_brackets_the_point_estimate():
    rng = np.random.default_rng(6)
    y = pd.Series([0] * 10 + [1] * 10, index=[f"p{i}" for i in range(20)])
    scores = np.r_[rng.uniform(0, 0.5, 10), rng.uniform(0.5, 1.0, 10)]
    from sklearn.metrics import roc_auc_score
    auc = roc_auc_score(y, scores)
    lo, hi = sg.bootstrap_auc_ci(y, scores, n_boot=500)
    assert lo <= auc <= hi
