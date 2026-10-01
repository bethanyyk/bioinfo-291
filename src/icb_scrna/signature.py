"""Stage 6 - candidate baseline signatures for responder stratification.

Nineteen baseline biopsies, one per patient, nine of them responders. At that
size the dominant risk is not finding too little - it is reporting an AUC that
describes the selection procedure rather than the biology. Three safeguards:

* **Every fitting step lives inside the cross-validation loop.** Gene selection,
  standardisation and the logistic fit are all refit on the training patients
  of each fold. A signature chosen on all 19 patients and then "validated" by
  cross-validation is scored on genes that already saw the held-out label.
* **Leave-one-patient-out, not leave-one-sample-out.** A patient contributing
  both a baseline and an on-treatment biopsy must be absent from the training
  fold entirely, or an on-treatment-derived program leaks their outcome.
* **A permutation null.** With n = 19 and several candidate signatures, the
  question is not whether AUC exceeds 0.5 but whether it exceeds what label
  shuffling produces under the identical pipeline.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


def clr(props: pd.DataFrame, pseudocount: float = 1e-3) -> pd.DataFrame:
    """Centred log-ratio transform - the standard way to put compositional
    fractions into a model that assumes unbounded, roughly linear features."""
    x = props.to_numpy(float) + pseudocount
    lx = np.log(x)
    return pd.DataFrame(lx - lx.mean(axis=1, keepdims=True),
                        index=props.index, columns=props.columns)


def _fit_predict(Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, seed: int) -> float:
    """Ridge logistic fit on standardised features, returning P(responder).

    Standardisation is done inline rather than through a Pipeline: this is
    called ~19,000 times per permutation test and the scikit-learn Pipeline
    wrapper dominates the runtime at this problem size. The fitted model is
    identical - centre and scale are still estimated on the training fold only,
    never on the held-out patient.
    """
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0)
    sd[sd == 0] = 1.0
    model = LogisticRegression(C=1.0, max_iter=500, random_state=seed)
    model.fit((Xtr - mu) / sd, ytr)
    return float(model.predict_proba((Xte - mu) / sd)[:, 1][0])


def loo_auc(X: pd.DataFrame, y: pd.Series, seed: int = 0) -> tuple[float, np.ndarray]:
    """Leave-one-out CV scores and the resulting AUC for a fixed feature set."""
    scores = np.empty(len(X))
    Xv, yv = X.to_numpy(float), y.to_numpy(int)
    for i in range(len(X)):
        tr = np.ones(len(X), bool); tr[i] = False
        scores[i] = _fit_predict(Xv[tr], yv[tr], Xv[~tr], seed)
    return float(roc_auc_score(yv, scores)), scores


def loo_auc_nested(profiles: pd.DataFrame, y: pd.Series, n_genes: int,
                   selection_profiles: pd.DataFrame | None = None,
                   selection_y: pd.Series | None = None,
                   selection_patient: pd.Series | None = None,
                   patient: pd.Series | None = None,
                   seed: int = 0) -> tuple[float, np.ndarray, list[str]]:
    """LOO-CV where the gene set is reselected from the training fold each time.

    ``profiles`` are the baseline pseudobulk profiles being predicted from.
    If ``selection_profiles`` is given the genes are chosen from that table
    instead (e.g. on-treatment profiles), with the held-out patient removed
    from it as well - otherwise their outcome leaks into the gene choice.
    Returns the AUC, the held-out scores, and the genes selected on the full
    data (reported for interpretation only, never used for scoring).
    """
    sel_X = profiles if selection_profiles is None else selection_profiles
    sel_y = y if selection_y is None else selection_y
    sel_p = (patient if selection_patient is None else selection_patient)
    patient = pd.Series(profiles.index, index=profiles.index) if patient is None else patient

    def pick(mask_sel: np.ndarray) -> list[str]:
        a = sel_X.loc[mask_sel & (sel_y == 1).to_numpy()]
        b = sel_X.loc[mask_sel & (sel_y == 0).to_numpy()]
        if len(a) < 2 or len(b) < 2:
            return []
        t, _ = stats.ttest_ind(a.to_numpy(), b.to_numpy(), equal_var=False)
        t = np.nan_to_num(t)
        order = np.argsort(-np.abs(t))[:n_genes]
        return list(sel_X.columns[order])

    scores = np.empty(len(profiles))
    yv = y.to_numpy(int)
    for i, held in enumerate(profiles.index):
        held_patient = patient.loc[held]
        keep_sel = (sel_p != held_patient).to_numpy()
        genes = pick(keep_sel)
        if not genes:
            scores[i] = 0.5
            continue
        tr = np.ones(len(profiles), bool); tr[i] = False
        Xtr = profiles.loc[:, genes].to_numpy(float)[tr]
        Xte = profiles.loc[:, genes].to_numpy(float)[~tr]
        scores[i] = _fit_predict(Xtr, yv[tr], Xte, seed)
    full_genes = pick(np.ones(len(sel_X), bool))
    return float(roc_auc_score(yv, scores)), scores, full_genes


def bootstrap_auc_ci(y: pd.Series, scores: np.ndarray, n_boot: int = 2000,
                     seed: int = 0) -> tuple[float, float]:
    """Percentile bootstrap CI for AUC, resampling patients."""
    rng = np.random.default_rng(seed)
    yv = y.to_numpy(int)
    out = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(yv), len(yv))
        if len(np.unique(yv[idx])) < 2:
            continue
        out.append(roc_auc_score(yv[idx], scores[idx]))
    if not out:
        return (np.nan, np.nan)
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)))


def permutation_p(auc: float, run_cv, y: pd.Series, n_perm: int = 1000,
                  seed: int = 0) -> float:
    """Fraction of label shuffles whose AUC, through the identical pipeline,
    reaches the observed value. ``run_cv(y_shuffled)`` must return an AUC."""
    rng = np.random.default_rng(seed)
    yv = y.to_numpy(int)
    null = np.empty(n_perm)
    for k in range(n_perm):
        null[k] = run_cv(pd.Series(rng.permutation(yv), index=y.index))
    return float((np.sum(null >= auc) + 1) / (n_perm + 1))
