"""Stage 3 - assign immune identities to Leiden clusters.

Panels are literature-canonical markers for the CD45+ compartment of a solid
tumour, cross-checked against the CELLxGENE CellGuide canonical marker lists
(``data/metadata/cellguide_canonical_markers.json``). Non-immune panels are
included deliberately: a CD45+ sort is never perfectly pure, and an unassigned
melanoma or stromal cluster must be identified rather than mislabelled as a
leucocyte.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

PANELS: dict[str, list[str]] = {
    "T cell (pan)": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC2"],
    "CD8 T": ["CD8A", "CD8B", "GZMK", "GZMH", "CCL5"],
    "CD4 T conv": ["CD4", "IL7R", "CD40LG", "ANXA1"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "TNFRSF4", "IKZF2"],
    "Exhausted/dysfunctional": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13", "ENTPD1"],
    "Memory/stem-like": ["TCF7", "SELL", "CCR7", "LEF1", "IL7R"],
    "Cycling": ["MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "FGFBP2"],
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19", "TNFRSF13C"],
    "Plasma": ["MZB1", "JCHAIN", "DERL3", "XBP1", "TNFRSF17", "IGHG1"],
    "Monocyte/Macrophage": ["LYZ", "CD14", "CD68", "C1QA", "C1QB", "CD163", "AIF1", "TYROBP"],
    "cDC": ["CD1C", "CLEC9A", "FCER1A", "CLEC10A", "XCR1", "BATF3"],
    "pDC": ["LILRA4", "CLEC4C", "IRF7", "IL3RA", "GZMB"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HPGDS"],
    "Melanocytic (non-immune)": ["MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"],
    "Stromal (non-immune)": ["COL1A1", "COL1A2", "DCN", "PECAM1", "VWF"],
}


def score_panels(adata: AnnData, panels: dict[str, list[str]] | None = None,
                 seed: int = 0) -> list[str]:
    """Add one ``score_genes`` column per panel. Returns the column names added."""
    panels = PANELS if panels is None else panels
    added = []
    for name, genes in panels.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        # h5ad keys cannot contain "/"
        col = "score:" + name.replace("/", "-")
        sc.tl.score_genes(adata, present, score_name=col, random_state=seed)
        added.append(col)
    return added


def cluster_panel_matrix(adata: AnnData, score_cols: list[str],
                         cluster_key: str = "leiden") -> pd.DataFrame:
    """Mean panel score per cluster, z-scored down each panel (across clusters)."""
    mat = adata.obs.groupby(cluster_key, observed=True)[score_cols].mean()
    mat.columns = [c.removeprefix("score:") for c in mat.columns]
    return (mat - mat.mean()) / mat.std(ddof=0)


def panel_gene_expression(adata: AnnData, genes: list[str],
                          cluster_key: str = "leiden") -> pd.DataFrame:
    """Mean ln(TPM+1) of individual genes per cluster, for manual inspection."""
    present = [g for g in genes if g in adata.var_names]
    sub = adata[:, present]
    dense = np.asarray(sub.X.todense()) if hasattr(sub.X, "todense") else np.asarray(sub.X)
    df = pd.DataFrame(dense, columns=present, index=adata.obs_names)
    df[cluster_key] = adata.obs[cluster_key].values
    return df.groupby(cluster_key, observed=True).mean()


# Cluster -> identity, decided from the marker evidence in
# results/cluster_annotation_evidence.csv. Leiden cluster "0" was a genuine
# NK/T mixture (only 25% of its cells were CD3+/NK-marker-negative) and was
# sub-clustered at resolution 0.4; "0,3" is the CD8 arm, "0,0"-"0,2" are NK.
CLUSTER_CELL_TYPE: dict[str, str] = {
    "0,0": "NK", "0,1": "NK", "0,2": "NK", "0,3": "CD8 T",
    "1": "CD8 T", "4": "CD8 T", "6": "CD8 T", "7": "CD8 T", "11": "CD8 T", "13": "CD8 T",
    "2": "CD4 T", "3": "CD4 T", "9": "CD4 T",
    "5": "Treg",
    "8": "Cycling T",
    "12": "gdT",
    "10": "pDC",
    "14": "cDC",
    "15": "Monocyte/Macrophage",
    "16": "Plasma", "19": "Plasma",
    "17": "B", "18": "B",
}

CLUSTER_CELL_STATE: dict[str, str] = {
    "0,0": "NK XCL1+ KLRC1-high", "0,1": "NK CD16+", "0,2": "NK CD16+ cytotoxic",
    # GNLY is 5.6 here: the highest of any CD8 state, but the LOWEST of the four
    # sub-clusters of cluster 0, so "GNLY-high" would be ambiguous as a label.
    "0,3": "CD8 effector cytotoxic",
    "1": "CD8 effector-memory GZMK+",
    "4": "CD8 exhausted",
    "6": "CD8 exhausted CXCL13+",
    "7": "CD8 transitional",
    "11": "CD8 exhausted CXCL13-high",
    "13": "CD8 interferon-stimulated",
    "2": "CD4 memory IL7R+ TCF7+",
    "3": "CD4 memory activated",
    "9": "CD4 CXCL13+ (Tfh-like)",
    "5": "Treg",
    "8": "Cycling T",
    "12": "gamma-delta T",
    "10": "pDC",
    "14": "cDC (CD1C+)",
    "15": "Monocyte/Macrophage",
    "16": "Plasma", "19": "Plasma IGHG-high",
    "17": "B (CD79A+ MS4A1-low)", "18": "B naive/memory",
}


def apply_labels(adata: AnnData, mapping: dict[str, str],
                 cluster_key: str = "leiden", key_added: str = "cell_type") -> AnnData:
    """Map cluster ids to cell-type labels, erroring on any unmapped cluster."""
    clusters = set(adata.obs[cluster_key].astype(str).unique())
    missing = clusters - set(mapping)
    if missing:
        raise KeyError(f"clusters without a label: {sorted(missing)}")
    adata.obs[key_added] = (
        adata.obs[cluster_key].astype(str).map(mapping).astype("category")
    )
    return adata
