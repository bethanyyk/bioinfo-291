"""Stage 5 - differential expression within cell types, on pseudobulk profiles.

**Why not DESeq2 here.** The negative-binomial model DESeq2 and edgeR fit is a
model of *counts*: its dispersion-mean trend and size factors are meaningless
applied to values that have already been divided by transcript length and
library size. GSE120575 publishes only a TPM matrix - the raw reads are under
dbGaP control (phs001680) - so there are no counts to give DESeq2. Rounding TPM
and calling it a count matrix would produce confident-looking output from a
misspecified variance model. Instead the pseudobulk profiles are compared on
the log scale they already live on.

**Patient-level pseudobulk.** Profiles are averaged to one per patient per
timepoint before testing, so no patient contributes two correlated observations
and the independence assumption of the two-sample test actually holds. This
costs a little power relative to a mixed model and buys a valid null.

The test is Welch's t-test per gene (unequal variances, which is the realistic
assumption across patients), with Benjamini-Hochberg correction within each
cell type x contrast.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats
from statsmodels.stats.multitest import multipletests

LN2 = float(np.log(2.0))


def pseudobulk(adata: AnnData, cell_type: str, group_keys: list[str],
               cell_type_key: str = "cell_type", min_cells: int = 10
               ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Mean ln(TPM+1) per (patient, group) for one cell type.

    Returns ``(profiles, meta)`` where ``profiles`` is observations x genes.
    Biopsies contributing fewer than ``min_cells`` cells of this type are
    dropped before averaging; a patient with several qualifying biopsies at the
    same timepoint is collapsed to their mean.
    """
    mask = adata.obs[cell_type_key].astype(str) == cell_type
    sub = adata[mask]
    per_sample = sub.obs.groupby("sample", observed=True).size()
    ok = per_sample[per_sample >= min_cells].index
    sub = sub[sub.obs["sample"].astype(str).isin(set(ok.astype(str)))]
    if sub.n_obs == 0:
        return pd.DataFrame(), pd.DataFrame()

    keys = ["patient"] + group_keys
    labels = sub.obs[keys].astype(str).agg("|".join, axis=1)
    frame = pd.DataFrame(
        np.asarray(sub.X.todense()), index=labels.values, columns=sub.var_names
    )
    profiles = frame.groupby(level=0).mean()
    meta = pd.DataFrame(
        [dict(zip(keys, k.split("|"))) for k in profiles.index], index=profiles.index
    )
    return profiles, meta


def welch_de(profiles: pd.DataFrame, group: pd.Series, test: str, reference: str,
             min_detect_frac: float = 0.25) -> pd.DataFrame:
    """Per-gene Welch t-test between two groups of pseudobulk profiles."""
    a = profiles.loc[group[group == test].index]
    b = profiles.loc[group[group == reference].index]
    if len(a) < 3 or len(b) < 3:
        return pd.DataFrame()
    detect = ((profiles > 0).mean(axis=0) >= min_detect_frac)
    a, b = a.loc[:, detect], b.loc[:, detect]
    t, p = stats.ttest_ind(a.to_numpy(), b.to_numpy(), equal_var=False)
    res = pd.DataFrame(
        {
            "gene": a.columns,
            f"mean_ln_{test}": a.mean(axis=0).to_numpy(),
            f"mean_ln_{reference}": b.mean(axis=0).to_numpy(),
            "log2FC": (a.mean(axis=0).to_numpy() - b.mean(axis=0).to_numpy()) / LN2,
            "t": t,
            "p_value": p,
            f"n_{test}": len(a),
            f"n_{reference}": len(b),
        }
    )
    res = res[np.isfinite(res["p_value"])]
    if res.empty:
        return res
    res["p_adj_BH"] = multipletests(res["p_value"], method="fdr_bh")[1]
    return res.sort_values("p_value").reset_index(drop=True)


def run_contrast(adata: AnnData, cell_types: list[str], group_key: str, test: str,
                 reference: str, subset: tuple[str, str] | None = None,
                 cell_type_key: str = "cell_type", min_cells: int = 10,
                 min_samples_per_group: int = 3) -> pd.DataFrame:
    """Run one contrast across every cell type with enough patients."""
    out = []
    for ct in cell_types:
        group_keys = [group_key] if subset is None else [group_key, subset[0]]
        profiles, meta = pseudobulk(adata, ct, group_keys, cell_type_key, min_cells)
        if profiles.empty:
            continue
        if subset is not None:
            keep = meta[subset[0]] == subset[1]
            profiles, meta = profiles.loc[keep.values], meta.loc[keep.values]
        g = meta[group_key]
        if (g == test).sum() < min_samples_per_group or (g == reference).sum() < min_samples_per_group:
            continue
        res = welch_de(profiles, g, test, reference)
        if res.empty:
            continue
        res.insert(0, "cell_type", ct)
        out.append(res)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()
