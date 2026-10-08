"""Stage 2 - QC, normalisation, integration, clustering, embedding.

The GEO matrix arrives as ``log2(TPM + 1)``. Two consequences drive this module:

1. Library size is already normalised away, so there is no total-count filter
   and no ``normalize_total`` step. Applying one would be a second
   normalisation of already-normalised data.
2. scanpy's fold-change and ``score_genes`` helpers assume *natural* log. Since
   ``ln(x+1) = log2(x+1) * ln 2``, the conversion is a single scalar multiply -
   exact, not an approximation, and it keeps every downstream statistic in the
   convention scanpy documents.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

LOG2_TO_LN = float(np.log(2.0))


def to_natural_log(adata: AnnData) -> AnnData:
    """Rescale ``log2(TPM+1)`` to ``ln(TPM+1)`` in place."""
    if adata.uns.get("log_base") == "e":
        return adata
    adata.X = adata.X.multiply(LOG2_TO_LN).tocsr().astype(np.float32)
    adata.uns["log_base"] = "e"
    return adata


def add_qc_metrics(adata: AnnData) -> AnnData:
    """Per-cell QC on the linear TPM scale recovered from the log values."""
    mito = adata.var_names.str.startswith("MT-")
    adata.var["mito"] = mito
    linear = adata.X.copy()
    linear.data = np.expm1(linear.data)
    total = np.asarray(linear.sum(axis=1)).ravel()
    mito_total = np.asarray(linear[:, mito].sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.diff(adata.X.indptr)
    adata.obs["total_tpm"] = total
    adata.obs["pct_mito"] = np.where(total > 0, 100.0 * mito_total / total, 0.0)
    adata.uns["n_mito_genes"] = int(mito.sum())
    return adata


def filter_cells(adata: AnnData, qc: dict) -> tuple[AnnData, pd.DataFrame]:
    """Apply the QC thresholds, returning the filtered object and a per-step tally."""
    keep_genes = adata.obs["n_genes"] >= qc["min_genes_per_cell"]
    keep_mito = adata.obs["pct_mito"] <= qc["max_pct_mito"]
    tally = pd.DataFrame(
        [
            {"step": "input", "cells": adata.n_obs},
            {"step": f"n_genes >= {qc['min_genes_per_cell']}", "cells": int(keep_genes.sum())},
            {"step": f"pct_mito <= {qc['max_pct_mito']}", "cells": int(keep_mito.sum())},
            {"step": "both", "cells": int((keep_genes & keep_mito).sum())},
        ]
    )
    return adata[keep_genes & keep_mito].copy(), tally


def build_atlas(adata: AnnData, atlas: dict, seed: int = 0) -> AnnData:
    """HVG selection, PCA, Harmony integration over patients, Leiden, UMAP."""
    sc.pp.highly_variable_genes(
        adata, n_top_genes=atlas["n_top_genes"], flavor="seurat", batch_key=atlas["harmony_key"]
    )
    adata.raw = adata
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=atlas["n_pcs"], svd_solver="arpack", random_state=seed)
    sc.external.pp.harmony_integrate(
        hvg, key=atlas["harmony_key"], basis="X_pca", adjusted_basis="X_pca_harmony",
        random_state=seed,
    )
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.obsm["X_pca_harmony"] = hvg.obsm["X_pca_harmony"]
    del hvg

    sc.pp.neighbors(adata, n_neighbors=atlas["neighbors_k"], use_rep="X_pca_harmony",
                    random_state=seed)
    sc.tl.leiden(adata, resolution=atlas["leiden_resolution"], key_added="leiden",
                 random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=seed)
    return adata


def patient_entropy(adata: AnnData, cluster_key: str = "leiden",
                    batch_key: str = "patient") -> pd.DataFrame:
    """Shannon entropy of the patient mix in each cluster, normalised to [0, 1].

    A cluster driven by one patient scores near 0; one drawn evenly from all
    patients scores near 1. This is the check that Harmony actually mixed
    patients rather than leaving per-patient islands.
    """
    tab = pd.crosstab(adata.obs[cluster_key], adata.obs[batch_key])
    p = tab.div(tab.sum(axis=1), axis=0).to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=1)
    max_ent = np.log(tab.shape[1])
    return pd.DataFrame(
        {
            "cluster": tab.index.astype(str),
            "n_cells": tab.sum(axis=1).to_numpy(),
            "n_patients": (tab > 0).sum(axis=1).to_numpy(),
            "top_patient_frac": p.max(axis=1),
            "normalised_entropy": ent / max_ent,
        }
    ).sort_values("normalised_entropy").reset_index(drop=True)
