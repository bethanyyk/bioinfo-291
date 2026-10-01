"""Baseline feature construction for the signature stage.

Kept separate from :mod:`icb_scrna.signature` so the *definition* of a candidate
predictor is testable without running any cross-validation, and so the
pipeline driver and an interactive session build identical features.

Every feature here is computed from **baseline (Pre) biopsies of the unbiased
CD45+ sort only**, and is indexed by patient rather than by biopsy: the
signature question is "can this patient's pre-treatment tumour predict their
outcome", so the patient is the unit.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from .signature import clr

MIN_CD8_CELLS = 10


def _baseline_cells(adata: AnnData, smeta: pd.DataFrame) -> pd.DataFrame:
    pre = smeta.index[smeta["timepoint"].astype(str) == "Pre"]
    obs = adata.obs
    obs = obs[(obs["enrichment"].astype(str) == "CD45")
              & (obs["sample"].astype(str).isin(set(pre)))]
    per_patient = obs.groupby("patient", observed=True)["sample"].nunique()
    if (per_patient > 1).any():
        bad = per_patient[per_patient > 1]
        raise ValueError(f"patients with >1 baseline biopsy, patient-level indexing unsafe: {dict(bad)}")
    return obs


def baseline_features(adata: AnnData, smeta: pd.DataFrame
                      ) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """Return ``({name: feature frame}, y)`` indexed by patient."""
    obs = _baseline_cells(adata, smeta)
    pat_of = obs.groupby("sample", observed=True)["patient"].first().astype(str)
    pre = smeta.index[smeta["timepoint"].astype(str) == "Pre"]
    y = (smeta.loc[pre, "response"].astype(str) == "Responder").astype(int)
    y.index = pat_of.reindex(pre).values

    comp = pd.crosstab(obs["sample"], obs["cell_type"])
    comp.index = pat_of.reindex(comp.index).values
    F_comp = clr(comp.div(comp.sum(axis=1), axis=0)).loc[y.index]

    cd8 = obs[obs["cell_type"].astype(str) == "CD8 T"]
    state = pd.crosstab(cd8["sample"], cd8["cell_state"])
    state.index = pat_of.reindex(state.index).values
    keep = state.sum(axis=1) >= MIN_CD8_CELLS
    F_state = clr(state.div(state.sum(axis=1), axis=0))[keep.values].reindex(y.index).dropna()

    tcf7 = np.asarray(adata[adata.obs_names.isin(cd8.index), "TCF7"].X.todense()).ravel()
    tdf = pd.DataFrame({"sample": cd8["sample"].astype(str).values, "pos": tcf7 > 0})
    frac = tdf.groupby("sample").agg(frac_TCF7pos=("pos", "mean"), n=("pos", "size"))
    frac = frac[frac["n"] >= MIN_CD8_CELLS]
    frac.index = pat_of.reindex(frac.index).values
    F_tcf7 = frac[["frac_TCF7pos"]].reindex(y.index).dropna()

    return (
        {
            "Whole-compartment composition (11 cell types, CLR)": F_comp,
            "CD8 state composition (CLR)": F_state,
            "TCF7+ fraction of CD8 cells [published comparator]": F_tcf7,
        },
        y,
    )


def cd8_pseudobulk_by_patient(adata: AnnData, smeta: pd.DataFrame, y: pd.Series,
                              min_detect_frac: float = 0.25
                              ) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series]:
    """Baseline and on-treatment CD8 pseudobulk, restricted to shared genes."""
    def build(tp: str) -> pd.DataFrame:
        m = ((adata.obs["cell_type"].astype(str) == "CD8 T")
             & (adata.obs["enrichment"].astype(str) == "CD45")
             & (adata.obs["timepoint"].astype(str) == tp))
        sub = adata[m]
        per = sub.obs.groupby("sample", observed=True).size()
        ok = set(per[per >= MIN_CD8_CELLS].index.astype(str))
        sub = sub[sub.obs["sample"].astype(str).isin(ok)]
        df = pd.DataFrame(np.asarray(sub.X.todense()),
                          index=sub.obs["patient"].astype(str).values,
                          columns=sub.var_names)
        prof = df.groupby(level=0).mean()
        return prof.loc[:, (prof > 0).mean(axis=0) >= min_detect_frac]

    pre = build("Pre").reindex(y.index).dropna()
    post = build("Post")
    resp_post = (smeta[smeta["timepoint"].astype(str) == "Post"]
                 .assign(r=lambda d: (d["response"].astype(str) == "Responder").astype(int))
                 .groupby("patient", observed=True)["r"].max())
    post = post.loc[post.index.isin(resp_post.index)]
    y_post = resp_post.reindex(post.index).astype(int)
    shared = pre.columns.intersection(post.columns)
    pat_post = pd.Series(post.index, index=post.index)
    return pre[shared], post[shared], y_post, pat_post
