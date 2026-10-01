"""Config loading, the preprocessing scalar conversion, and the DE helper."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

from icb_scrna import de as DE
from icb_scrna import preprocess as pp
from icb_scrna.config import load_config


def test_config_loads_and_resolves_paths():
    cfg = load_config()
    assert cfg["dataset"]["accession"] == "GSE120575"
    assert cfg.path("results").is_dir()
    assert cfg.path("results", "x.csv").name == "x.csv"


def test_config_rejects_unknown_path_kind():
    cfg = load_config()
    with pytest.raises(KeyError, match="unknown path kind"):
        cfg.path("nowhere")


def test_log2_to_natural_log_conversion_is_exact():
    tpm = np.array([[0.0, 9.0, 999.0]])
    adata = AnnData(X=sp.csr_matrix(np.log2(tpm + 1).astype(np.float32)))
    pp.to_natural_log(adata)
    np.testing.assert_allclose(adata.X.toarray(), np.log(tpm + 1), rtol=1e-5)
    assert adata.uns["log_base"] == "e"


def test_to_natural_log_is_idempotent():
    adata = AnnData(X=sp.csr_matrix(np.array([[1.0, 2.0]], dtype=np.float32)))
    pp.to_natural_log(adata)
    first = adata.X.toarray().copy()
    pp.to_natural_log(adata)
    np.testing.assert_allclose(adata.X.toarray(), first)


def test_qc_metrics_computes_mito_fraction_on_the_linear_scale():
    # one cell: MT gene at TPM 300, other at TPM 700 -> 30% mitochondrial
    tpm = np.array([[300.0, 700.0]])
    adata = AnnData(X=sp.csr_matrix(np.log(tpm + 1).astype(np.float32)),
                    var=pd.DataFrame(index=["MT-CO1", "ACTB"]))
    pp.add_qc_metrics(adata)
    assert abs(float(adata.obs["pct_mito"].iloc[0]) - 30.0) < 0.2
    assert int(adata.obs["n_genes"].iloc[0]) == 2


def test_filter_cells_applies_both_thresholds():
    adata = AnnData(X=sp.csr_matrix(np.ones((4, 3), dtype=np.float32)))
    adata.obs["n_genes"] = [100, 1000, 1000, 100]
    adata.obs["pct_mito"] = [5.0, 5.0, 50.0, 50.0]
    kept, tally = pp.filter_cells(adata, {"min_genes_per_cell": 500, "max_pct_mito": 20.0})
    assert kept.n_obs == 1
    assert tally.loc[tally["step"] == "input", "cells"].iloc[0] == 4


def test_patient_entropy_flags_a_single_patient_cluster():
    obs = pd.DataFrame({
        "leiden": ["0"] * 10 + ["1"] * 10,
        "patient": ["p1"] * 10 + [f"p{i%5}" for i in range(10)],
    })
    adata = AnnData(X=sp.csr_matrix(np.zeros((20, 1), dtype=np.float32)), obs=obs)
    ent = pp.patient_entropy(adata).set_index("cluster")
    assert ent.loc["0", "normalised_entropy"] == 0.0
    assert ent.loc["1", "normalised_entropy"] > 0.9


def test_welch_de_recovers_direction_and_requires_three_per_group():
    rng = np.random.default_rng(0)
    idx = [f"s{i}" for i in range(10)]
    X = rng.normal(size=(10, 20))
    X[5:, 0] += 5.0
    profiles = pd.DataFrame(np.abs(X), index=idx, columns=[f"g{i}" for i in range(20)])
    grp = pd.Series(["ref"] * 5 + ["test"] * 5, index=idx)
    res = DE.welch_de(profiles, grp, "test", "ref").set_index("gene")
    assert res.loc["g0", "log2FC"] > 0
    assert res.loc["g0", "p_adj_BH"] < 0.05
    tiny = DE.welch_de(profiles.iloc[:4], grp.iloc[:4], "test", "ref")
    assert tiny.empty


def test_welch_de_log2fc_matches_the_mean_difference():
    idx = ["a", "b", "c", "d", "e", "f"]
    profiles = pd.DataFrame({"g0": [1.0, 1.1, 0.9, 3.0, 3.1, 2.9]}, index=idx)
    grp = pd.Series(["ref"] * 3 + ["test"] * 3, index=idx)
    res = DE.welch_de(profiles, grp, "test", "ref").set_index("gene")
    expected = (3.0 - 1.0) / np.log(2)
    assert abs(res.loc["g0", "log2FC"] - expected) < 0.05
