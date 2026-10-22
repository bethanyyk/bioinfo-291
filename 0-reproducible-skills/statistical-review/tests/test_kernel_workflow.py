"""
Tests for the kernel-level workflow: verdict assignment, per-claim review, and
the fixed figure spine.
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import kernel as K
from utilities import resampling as rs


# --------------------------------------------------------------------------
# verdict cascade
# --------------------------------------------------------------------------

def _routes(pvals):
    """Minimal three-routes frame with the columns propose_verdict reads."""
    return pd.DataFrame({
        "route": [f"r{i}" for i in range(len(pvals))],
        "test": [f"t{i}" for i in range(len(pvals))],
        "p_value": pvals,
        "assumes": ["something"] * len(pvals),
        "min_achievable_p": [0.0] * len(pvals),
        "reachable_at_0.05": [True] * len(pvals),
        "agrees_at_0.05": [p < 0.05 for p in pvals],
    })


def test_no_recomputation_is_not_assessable():
    out = K.propose_verdict({})
    assert out["verdict"] == "not-assessable"


def test_opposite_sign_effect_is_contradicted():
    out = K.propose_verdict({
        "routes": _routes([0.01, 0.01, 0.01]),
        "permutation": {"floor_binds": False},
        "effect_reported": 2.0, "effect_recomputed": -1.8,
    })
    assert out["verdict"] == "contradicted"
    assert "opposite sign" in out["reasons"][0]


def test_crossing_alpha_with_unanimous_routes_is_not_reproducible():
    out = K.propose_verdict({
        "routes": _routes([0.4, 0.4, 0.4]),
        "permutation": {"floor_binds": False},
        "p_reported": 0.001, "p_recomputed": 0.4,
    })
    assert out["verdict"] == "not-reproducible"
    assert "opposite sides" in out["reasons"][0]
    assert "every route agrees" in out["reasons"][0]


def test_binding_floor_outranks_the_discrepancy_it_explains():
    # A design that cannot reach alpha explains the mismatch rather than merely
    # co-occurring with it, so the floor is the verdict and the mismatch is a reason.
    out = K.propose_verdict({
        "routes": _routes([0.09, 0.0625, 0.09]),
        "permutation": {"floor_binds": True, "min_achievable_p": 0.0625},
        "p_reported": 0.04, "p_recomputed": 0.09,
    })
    assert out["verdict"] == "underpowered-resolution"
    assert "floor" in out["reasons"][0]
    assert any("0.04" in r for r in out["reasons"])


def test_disagreeing_routes_are_assumption_dependent_not_irreproducible():
    # The paper's own route (p = 0.03) reproduces; the rank and permutation
    # routes do not. The claim depends on the modelling choice, which is a
    # different finding from the paper's number failing to check out.
    out = K.propose_verdict({
        "routes": _routes([0.03, 0.08, 0.07]),
        "permutation": {"floor_binds": False},
        "p_reported": 0.03, "p_recomputed": 0.07,
    })
    assert out["verdict"] == "assumption-dependent"
    assert "disagree at alpha" in out["reasons"][0]
    assert any("agrees with some routes" in r for r in out["reasons"])


def test_clustering_downgrades_to_supported_with_caveat():
    out = K.propose_verdict({
        "routes": _routes([0.01, 0.01, 0.01]),
        "permutation": {"floor_binds": False},
        "p_reported": 0.01, "p_recomputed": 0.012,
        "independence": {"pseudoreplication_risk": True, "note": "ICC is high"},
    })
    assert out["verdict"] == "supported-with-caveat"
    assert out["reasons"] == ["ICC is high"]


def test_sem_error_bars_downgrade_to_supported_with_caveat():
    out = K.propose_verdict({
        "routes": _routes([0.01, 0.01, 0.01]),
        "permutation": {"floor_binds": False},
        "error_bar": {"resolved": True, "best_match": "SEM", "sd_to_sem_factor": 5.0},
    })
    assert out["verdict"] == "supported-with-caveat"
    assert "SEM" in out["reasons"][0]


def test_clean_evidence_is_supported():
    out = K.propose_verdict({
        "routes": _routes([0.001, 0.002, 0.0015]),
        "permutation": {"floor_binds": False},
        "p_reported": 0.001, "p_recomputed": 0.0015,
        "independence": {"pseudoreplication_risk": False},
    })
    assert out["verdict"] == "supported"


def test_every_proposed_verdict_is_in_the_vocabulary():
    cases = [
        {},
        {"routes": _routes([0.4]), "permutation": {"floor_binds": False},
         "p_reported": 0.001, "p_recomputed": 0.4},
        {"routes": _routes([0.09]), "permutation": {"floor_binds": True,
                                                    "min_achievable_p": 0.0625}},
        {"routes": _routes([0.001]), "permutation": {"floor_binds": False}},
    ]
    for case in cases:
        assert K.propose_verdict(case)["verdict"] in K.VERDICT_ORDER


# --------------------------------------------------------------------------
# review_claim
# --------------------------------------------------------------------------

@pytest.fixture
def clustered():
    """96 fields from 8 animals, with large between-animal variance."""
    rng = np.random.default_rng(0)
    rows = []
    for unit in range(8):
        arm = "treated" if unit < 4 else "control"
        effect = rng.normal(-0.35 if arm == "treated" else 0.0, 0.55)
        for field in range(12):
            rows.append({"mouse": f"M{unit}", "arm": arm,
                         "value": float(np.exp(5.2 + effect + rng.normal(0, 0.12)))})
    return pd.DataFrame(rows)


def test_review_claim_returns_every_diagnostic(clustered):
    out = K.review_claim(clustered, "value", "arm", unit_col="mouse", claim_id="C1",
                         claim_short="treated is lower", p_reported=0.001,
                         n_perm=1000, n_boot=1000, random_state=0)
    for key in ("routes", "permutation", "bootstrap", "independence", "shape",
                "variance", "verdict", "reasons", "ledger_row"):
        assert key in out, key
    assert out["verdict"] in K.VERDICT_ORDER
    assert len(out["routes"]) == 4
    assert out["independence"]["n_units"] == 8
    assert out["independence"]["pseudoreplication_risk"] is True
    row = out["ledger_row"].iloc[0]
    assert row["claim_id"] == "C1"
    assert np.isfinite(row["p_recomputed"])
    assert row["ci_low"] < row["effect_recomputed"] < row["ci_high"]


def test_unit_level_retest_loses_the_significance_clustering_manufactured(clustered):
    per_unit = clustered.groupby(["mouse", "arm"], as_index=False)["value"].mean()
    field_level = K.review_claim(clustered, "value", "arm", claim_id="F",
                                 n_perm=2000, n_boot=500, random_state=0)
    unit_level = K.review_claim(per_unit, "value", "arm", claim_id="U",
                                n_perm=2000, n_boot=500, random_state=0)
    assert field_level["p_recomputed"] < unit_level["p_recomputed"]
    assert unit_level["permutation"]["min_achievable_p"] == pytest.approx(
        rs.min_achievable_p(4, 4))


def test_paired_review_uses_the_pair_column():
    rng = np.random.default_rng(1)
    base = rng.normal(180, 25, 5)
    change = np.array([-18.0, -12.0, -9.0, -15.0, 4.0])
    df = pd.DataFrame({
        "subject": np.repeat([f"S{i}" for i in range(5)], 2),
        "arm": ["baseline", "treated"] * 5,
        "value": np.ravel([[b, b + c] for b, c in zip(base, change)]),
    })
    out = K.review_claim(df, "value", "arm", paired=True, pair_col="subject",
                         claim_id="C2", p_reported=0.040, n_perm=2000, n_boot=1000,
                         random_state=0)
    assert out["paired"] is True
    assert out["differences"].size == 5
    # 5 pairs: the sign-flip floor is 2/2**5 = 0.0625, above alpha
    assert out["permutation"]["min_achievable_p"] == pytest.approx(0.0625)
    assert out["verdict"] == "underpowered-resolution"
    assert set(out["routes"]["test"]) == {"paired t-test", "Wilcoxon signed-rank",
                                          "sign test",
                                          "sign-flip permutation (mean difference)"}


def test_review_claim_rejects_a_three_level_grouping(clustered):
    three = clustered.copy()
    three.loc[three.index[:10], "arm"] = "third"
    with pytest.raises(ValueError, match="exactly 2 levels"):
        K.review_claim(three, "value", "arm", claim_id="X", n_perm=100)


# --------------------------------------------------------------------------
# the figure spine
# --------------------------------------------------------------------------

def test_figure_spine_is_fixed_and_complete(tmp_path, clustered):
    claim = K.review_claim(clustered, "value", "arm", unit_col="mouse", claim_id="C1",
                           n_perm=500, n_boot=500, random_state=0,
                           figure_path=str(tmp_path / "panel_C1.png"))
    rng = np.random.default_rng(2)
    audit = {"pvalues": np.concatenate([rng.uniform(0, 0.01, 5), rng.uniform(0, 1, 40)])}
    figs = K.build_review_figures(clustered, "mouse", "arm", "value", claims=[claim],
                                  audit=audit, outdir=str(tmp_path / "figs"),
                                  value_label="volume (mm3)")
    for key, _ in K.FIGURE_SPINE:
        assert key in figs, key
        assert os.path.getsize(figs[key]) > 5000
    assert "panel_C1" in figs


def test_missing_inputs_still_render_every_slot(tmp_path):
    # No data, no claims, no p-values: every slot must still exist, each one
    # naming why it could not be built.
    figs = K.build_review_figures(None, None, None, None, claims=[], audit=None,
                                  outdir=str(tmp_path / "figs"))
    for key, _ in K.FIGURE_SPINE:
        assert os.path.getsize(figs[key]) > 2000


def test_verdict_counts_covers_the_whole_vocabulary(clustered):
    claim = K.review_claim(clustered, "value", "arm", claim_id="C1", n_perm=300,
                           n_boot=300, random_state=0)
    counts = K.verdict_counts([claim])
    assert list(counts["verdict"]) == list(K.VERDICT_ORDER)
    assert counts["n_claims"].sum() == 1


# --------------------------------------------------------------------------
# PDF round trip
# --------------------------------------------------------------------------

def test_load_paper_reads_a_generated_pdf(tmp_path):
    sys.path.insert(0, os.path.join(_ROOT, "examples"))
    from run_demo_review import build_demo_pdf

    pdf = str(tmp_path / "m.pdf")
    build_demo_pdf(pdf)
    paper = K.load_paper(pdf)
    assert paper["text_quality"]["usable_text_layer"] is True
    assert len(paper["inventory"]["tests"]) >= 3
    assert len(paper["inventory"]["p_values"]) >= 20
    assert "welch_t" not in set(paper["inventory"]["named_tests"]["test_named"])
    assert "student_t" in set(paper["inventory"]["named_tests"]["test_named"])

    audit = K.audit_reported_statistics(paper["inventory"])
    assert audit["multiplicity"]["m_tests"] >= 20
    assert audit["just_below"]["interpretable"] is True
    # the planted cluster of p-values just under 0.05
    assert audit["just_below"]["n_in_suspect_window"] >= 4
