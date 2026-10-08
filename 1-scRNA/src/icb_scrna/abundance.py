"""Stage 4 - which populations expand or contract.

Cell-type proportions are compositional, overdispersed, and measured on a
handful of biopsies per group, with some patients contributing more than one
biopsy. Three design decisions follow from that:

* **The biopsy is the unit, never the cell.** Treating 14,503 cells as 14,503
  independent observations would turn a 48-biopsy study into a p-value factory.
* **Binomial counts with a cluster-robust covariance.** Each cell type is
  modelled as cells-of-type vs cells-of-other per biopsy with
  ``statsmodels`` GLM(Binomial), and standard errors are clustered on patient.
  The sandwich estimator absorbs both the overdispersion that a plain binomial
  would understate and the correlation between repeat biopsies of one patient.
* **Lineage-enriched sorts are excluded.** Biopsies sorted on a T-cell or
  myeloid marker have proportions fixed by the sort, not by the tumour.

The within-patient paired analysis answers a different question from the
group-level GLM - "did this patient's compartment change" rather than "do the
groups differ" - and is reported alongside it, not as a confirmation of it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests


def proportion_table(obs: pd.DataFrame, cell_type_key: str = "cell_type",
                     sample_key: str = "sample",
                     unbiased_only: bool = True,
                     min_cells: int = 20) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-biopsy counts and proportions over cell types.

    Returns ``(counts, proportions)``, both biopsies x cell types. Biopsies with
    fewer than ``min_cells`` usable cells are dropped.
    """
    df = obs
    if unbiased_only:
        df = df[df["enrichment"].astype(str) == "CD45"]
    counts = pd.crosstab(df[sample_key], df[cell_type_key])
    counts = counts.loc[counts.sum(axis=1) >= min_cells]
    props = counts.div(counts.sum(axis=1), axis=0)
    return counts, props


def _design(meta: pd.DataFrame, group: str, test: str, reference: str,
            covariates: list[str] | None) -> tuple[pd.DataFrame, pd.Series]:
    keep = meta[group].astype(str).isin([test, reference])
    meta = meta.loc[keep]
    X = pd.DataFrame({"Intercept": 1.0}, index=meta.index)
    X[f"{group}[{test}]"] = (meta[group].astype(str) == test).astype(float)
    for cov in covariates or []:
        dummies = pd.get_dummies(meta[cov].astype(str), prefix=cov, drop_first=True)
        X = X.join(dummies.astype(float))
    # Drop constant covariate columns, but never the intercept - without it the
    # group coefficient is an absolute log-odds rather than a contrast.
    keep = [c for c in X.columns if c == "Intercept" or X[c].std(ddof=0) > 0]
    X = X[keep]
    return X, meta.index.to_series()


def test_abundance(counts: pd.DataFrame, sample_meta: pd.DataFrame, group: str,
                   test: str, reference: str, cluster_key: str = "patient",
                   covariates: list[str] | None = None,
                   alpha: float = 0.05) -> pd.DataFrame:
    """Binomial GLM per cell type with patient-clustered robust standard errors.

    ``counts`` is biopsies x cell types; ``sample_meta`` is indexed by biopsy and
    carries ``group`` and ``cluster_key``.
    """
    meta = sample_meta.loc[counts.index]
    X, idx = _design(meta, group, test, reference, covariates)
    counts = counts.loc[idx]
    groups = meta.loc[idx, cluster_key].astype(str)
    total = counts.sum(axis=1)
    coef_name = f"{group}[{test}]"

    rows = []
    for ct in counts.columns:
        n_type = counts[ct].to_numpy(float)
        endog = np.column_stack([n_type, total.to_numpy(float) - n_type])
        try:
            fit = sm.GLM(endog, X.to_numpy(float), family=sm.families.Binomial()).fit(
                cov_type="cluster", cov_kwds={"groups": groups.to_numpy()}
            )
            j = list(X.columns).index(coef_name)
            beta, se, p = fit.params[j], fit.bse[j], fit.pvalues[j]
        except Exception:  # singular design for a near-absent cell type
            beta = se = p = np.nan
        mask_t = X[coef_name].to_numpy() == 1.0
        prop = (counts[ct] / total).to_numpy()
        rows.append(
            {
                "cell_type": ct,
                "n_biopsies_test": int(mask_t.sum()),
                "n_biopsies_ref": int((~mask_t).sum()),
                f"mean_prop_{test}": prop[mask_t].mean(),
                f"mean_prop_{reference}": prop[~mask_t].mean(),
                "log_odds_ratio": beta,
                "se": se,
                "odds_ratio": np.exp(beta) if np.isfinite(beta) else np.nan,
                "p_value": p,
            }
        )
    res = pd.DataFrame([r for r in rows if "cell_type" in r])
    ok = res["p_value"].notna()
    res["p_adj_BH"] = np.nan
    if ok.any():
        res.loc[ok, "p_adj_BH"] = multipletests(res.loc[ok, "p_value"], method="fdr_bh")[1]
    res["significant"] = res["p_adj_BH"] < alpha
    res["direction"] = np.where(res["log_odds_ratio"] > 0, "expands", "contracts")
    res["contrast"] = f"{test} vs {reference}"
    return res.sort_values("p_value").reset_index(drop=True)


def paired_change(props: pd.DataFrame, sample_meta: pd.DataFrame,
                  within: str = "patient", group: str = "timepoint",
                  test: str = "Post", reference: str = "Pre") -> pd.DataFrame:
    """Within-patient Wilcoxon signed-rank on paired biopsies.

    Patients contributing several biopsies at one timepoint are averaged to one
    value per timepoint first, so each patient contributes exactly one pair.
    """
    meta = sample_meta.loc[props.index]
    long = props.join(meta[[within, group]])
    per = long.groupby([within, group], observed=True).mean(numeric_only=True)
    wide = per.unstack(group)
    rows = []
    for ct in props.columns:
        try:
            pair = wide[ct][[reference, test]].dropna()
        except KeyError:
            continue
        if len(pair) < 5:
            continue
        diff = pair[test] - pair[reference]
        stat, p = stats.wilcoxon(pair[test], pair[reference])
        rows.append(
            {
                "cell_type": ct,
                "n_patients_paired": len(pair),
                f"median_prop_{reference}": pair[reference].median(),
                f"median_prop_{test}": pair[test].median(),
                "median_change": diff.median(),
                "n_increased": int((diff > 0).sum()),
                "statistic": stat,
                "p_value": p,
            }
        )
    res = pd.DataFrame(rows)
    if not res.empty:
        res["p_adj_BH"] = multipletests(res["p_value"], method="fdr_bh")[1]
        res = res.sort_values("p_value").reset_index(drop=True)
    return res
