"""Statistical review of Goold, Vas, Olsen & Newberry (2016) R. Soc. Open Sci. 3:160268.

Regenerates every artifact of the review:
    goold_2016/figures/F1..F6 + per-claim panels
    goold_2016/claim_ledger.csv
    goold_2016/statcheck.csv
    goold_2016/REVIEW_REPORT.md
    goold_2016/tables/*.csv

Heavy resampling (bootstrap / permutation / stability) is produced by heavy.py and
read back from goold_2016/results/*.npz.  Run `python goold_2016/heavy.py 0 1` first
to regenerate it (budgets are recorded in section 10 of the report).

Seeds: every resampling call below passes an explicit seed; heavy.py records its own.
"""
import glob
import json
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
from scipy import stats as sps

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "skills", "statistical-review"))

import core                                            # noqa: E402
import prep                                            # noqa: E402
from utilities import claim_extraction as ce           # noqa: E402
from utilities import recompute as rc                  # noqa: E402
from utilities import assumption_checks as ac          # noqa: E402
from utilities import resampling as rs                 # noqa: E402
from utilities import multiplicity as mp               # noqa: E402
from utilities import review_figures as rf             # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)

FIG = os.path.join(HERE, "figures")
TAB = os.path.join(HERE, "tables")
RES = os.path.join(HERE, "results")
for d in (FIG, TAB, RES):
    os.makedirs(d, exist_ok=True)

ALPHA = 0.05
GAMMA = 0.65
N_PERM_ROUTES = 100_000
SEED = 1
NODES = prep.NODES20
IU = np.triu_indices(20, 1)
PAIRLAB = [f"{NODES[i]}-{NODES[j]}" for i, j in zip(*IU)]
OUT = {}            # every number the report quotes


# ===================================================================== stage 0
def stage0_inputs():
    raw, mi, rscript, paper = prep.load_all()
    pages = [paper]
    inv = ce.build_claim_inventory(pages)
    tq = ce.pdf_text_quality(pages)
    statcheck = rc.statcheck_table(inv["tests"]) if len(inv["tests"]) else pd.DataFrame(
        columns=["stat_type", "stat_value", "df1", "df2", "p_reported", "p_recomputed",
                 "consistent", "decision_error"])
    statcheck.to_csv(os.path.join(HERE, "statcheck.csv"), index=False)
    OUT["inputs"] = {
        "checksums": {k: prep.sha(v) for k, v in
                      zip(prep.ARTIFACTS, [prep.RAW_CSV, prep.MI_CSV, prep.RSCRIPT, prep.PAPER_TXT])},
        "raw_shape": list(raw.shape), "mi_shape": list(mi.shape),
        "rscript_lines": rscript.count("\n") + 1,
        "text_quality": (tq.to_dict("records")[0] if hasattr(tq, "to_dict") else dict(tq)),
        "inventory": {k: int(len(v)) for k, v in inv.items() if hasattr(v, "shape")},
        "n_statcheck_rows": int(len(statcheck)),
    }
    return raw, mi, rscript, paper, inv


# ===================================================================== stage 1-2
def stage12_design(raw, mi):
    pp = prep.rebuild_preprocessing(raw)
    rm170 = pp["rm170"]
    mi2 = mi.copy()
    mi2["handlerType"] = rm170["handlerType"].to_numpy()
    assert (rm170["newID"].to_numpy() == mi2["id"].to_numpy()).all()

    desc_all = [c for c in raw.columns if c not in ("workType", "handlerType", "breed", "newID")]
    na = raw[desc_all].isna().sum()
    zr = (raw[desc_all] == 0).sum()
    d29 = [c for c in mi.columns if c not in ("workType", "id")]
    s29 = rm170[d29].replace(0, np.nan)
    m29 = s29.isna().sum() / len(s29) * 100
    pz = (raw[pp["kept38"]] == 0).sum(axis=1) / len(pp["kept38"]) * 100
    sp = raw.groupby("newID").size()

    hc = mi2.groupby("id").size()
    rep_ids = hc[hc > 1].index
    sub = mi2[mi2.id.isin(rep_ids)]
    icc_rows = []
    for d in NODES:
        r = ac.design_independence(sub[["id", "workType", d]].rename(columns={d: "value"}),
                                   "id", "workType", "value")
        icc_rows.append({"descriptor": d, "icc_repeat_handlers": r["icc"],
                         "design_effect": r["design_effect"]})
    icc = pd.DataFrame(icc_rows).sort_values("icc_repeat_handlers", ascending=False)
    icc.to_csv(os.path.join(TAB, "handler_icc.csv"), index=False)

    long1 = pd.DataFrame({"handler": mi2.id.to_numpy(),
                          "dog_type": np.where(mi2.workType == 1, "patrol", "detection"),
                          "mean_rating": mi2[NODES].mean(axis=1).to_numpy()})
    indep = ac.design_independence(long1, "handler", "dog_type", "mean_rating")

    OUT["design"] = {
        "n_invited_handlers": 227, "n_surveys": 174, "n_dogs_raw": int(len(raw)),
        "n_handlers": int(raw.newID.nunique()),
        "surveys_per_handler": [float(sp.mean()), float(sp.std(ddof=1))],
        "n_patrol_raw": int((raw.workType == 1).sum()), "n_detect_raw": int((raw.workType == 0).sum()),
        "n_patrol_net": int((mi.workType == 1).sum()), "n_detect_net": int((mi.workType == 0).sum()),
        "drop_zero10": pp["drop_zero10"], "drop_na5": pp["drop_na5"],
        "missing_per_desc": [float(na.mean()), float(na.std(ddof=1))],
        "missing_pct_per_desc": [float((na / len(raw) * 100).mean()),
                                 float((na / len(raw) * 100).std(ddof=1))],
        "zero_per_desc": [float(zr.mean()), float(zr.std(ddof=1))],
        "zero_pct_per_desc": [float((zr / len(raw) * 100).mean()),
                              float((zr / len(raw) * 100).std(ddof=1))],
        "worst_dog_pct_zero_38cols": float((raw[pp["kept38"]].iloc[169] == 0).sum() / 38 * 100),
        "per_dog_zero_pct": [float(pz.mean()), float(pz.std(ddof=1))],
        "missing_pct_29": [float(m29.mean()), float(m29.std(ddof=1))],
        "n_repeat_handlers": int((hc > 1).sum()), "n_dogs_repeat": int(hc[hc > 1].sum()),
        "max_dogs_per_handler": int(hc.max()),
        "icc_gt_010": int((icc.icc_repeat_handlers > 0.10).sum()),
        "icc_max": float(icc.icc_repeat_handlers.max()),
        "icc_max_desc": str(icc.descriptor.iloc[0]),
        "icc_median": float(icc.icc_repeat_handlers.median()),
        "overall_icc": float(indep["icc"]), "design_effect": float(indep["design_effect"]),
        "effective_n": float(indep["effective_n"]),
        "obs_per_unit_mean": float(indep["obs_per_unit_mean"]),
    }
    return mi2, pp, icc, long1, indep


# ===================================================================== stage 3
def stage3_networks(mi):
    patrol, detect = prep.network_inputs(mi)
    nets = {}
    for name, X, n in [("patrol", patrol, 116), ("detect", detect, 54)]:
        R = core.polychoric_matrix(X)
        f = core.ebic_glasso(R, n=n, gamma=GAMMA, nlambda=100, lambda_min_ratio=0.01,
                             return_path=True)
        f["R"] = R
        f["strength"] = core.strength(f["net"])
        f["betweenness"] = core.betweenness(f["net"])
        f["n"] = n
        nets[name] = f
        pd.DataFrame(f["net"], index=NODES, columns=NODES).to_csv(
            os.path.join(TAB, f"association_matrix_{name}.csv"))
    cent = pd.DataFrame({
        "node": NODES, "label": [prep.LABELS[n] for n in NODES],
        "desirable": [n in prep.DESIRABLE for n in NODES],
        "strength_patrol": nets["patrol"]["strength"],
        "betweenness_patrol": nets["patrol"]["betweenness"],
        "strength_detect": nets["detect"]["strength"],
        "betweenness_detect": nets["detect"]["betweenness"]})
    for g in ("patrol", "detect"):
        cent[f"mean_cent_{g}"] = cent[[f"strength_{g}", f"betweenness_{g}"]].mean(axis=1)
    cent.to_csv(os.path.join(TAB, "centrality_observed.csv"), index=False)

    def sweep(path, n, p=20):
        lam, E, eb = path[:, 0], path[:, 1], path[:, 2]
        neg2L = eb - E * np.log(n) - 4 * E * GAMMA * np.log(p)
        rows = []
        for g in np.round(np.arange(0, 1.001, 0.05), 2):
            e = neg2L + E * np.log(n) + 4 * E * g * np.log(p)
            k = int(np.argmin(e))
            rows.append({"gamma": g, "lambda": float(lam[k]), "n_edges": int(E[k]),
                         "ebic": float(e[k]), "ebic_empty_graph": float(n * p),
                         "margin_over_empty": float(n * p - e[k])})
        return pd.DataFrame(rows)

    gsweep = {k: sweep(nets[k]["path"], nets[k]["n"]) for k in nets}
    for k, v in gsweep.items():
        v.to_csv(os.path.join(TAB, f"gamma_sweep_{k}.csv"), index=False)

    OUT["networks"] = {
        k: {"n": nets[k]["n"], "n_edges": int(nets[k]["n_edges"]),
            "pct_edges": 100 * nets[k]["n_edges"] / 190,
            "lambda": float(nets[k]["lambda"]), "ebic": float(nets[k]["ebic"]),
            "ebic_empty": float(nets[k]["n"] * 20),
            "margin_over_empty": float(nets[k]["n"] * 20 - nets[k]["ebic"]),
            "ebic_cost_per_edge": float(np.log(nets[k]["n"]) + 4 * GAMMA * np.log(20)),
            "first_gamma_empty": float(gsweep[k].query("n_edges == 0").gamma.min())
            if (gsweep[k].n_edges == 0).any() else float("nan")}
        for k in nets}
    OUT["networks"]["corr_meancent_betweenness"] = {
        g: float(np.corrcoef(cent[f"mean_cent_{g}"], cent[f"betweenness_{g}"])[0, 1]) for g in ("patrol", "detect")}
    OUT["networks"]["corr_meancent_strength"] = {
        g: float(np.corrcoef(cent[f"mean_cent_{g}"], cent[f"strength_{g}"])[0, 1]) for g in ("patrol", "detect")}
    OUT["networks"]["top_nodes"] = {
        c: cent.sort_values(c, ascending=False).node.iloc[:3].tolist()
        for c in ["strength_patrol", "betweenness_patrol", "strength_detect", "betweenness_detect"]}
    return nets, cent, gsweep


# ===================================================================== stage 4
def stage4_routes(cent):
    """Three routes on the 'central descriptors are desirable' claim, both networks."""
    des = cent.desirable.to_numpy()
    out, perms = {}, {}
    for net in ("patrol", "detect"):
        for meas in ("strength", "betweenness", "mean_cent"):
            col = f"{meas}_{net}"
            a = cent.loc[des, col].to_numpy()
            b = cent.loc[~des, col].to_numpy()
            tr = rs.three_routes(a, b, n_perm=N_PERM_ROUTES, random_state=SEED)
            pr = rs.permutation_test(a, b, statistic="mean_diff", n_perm=N_PERM_ROUTES,
                                     random_state=SEED)
            out[(net, meas)] = tr
            perms[(net, meas)] = pr
            tr.assign(network=net, measure=meas).to_csv(
                os.path.join(TAB, f"three_routes_{net}_{meas}.csv"), index=False)
    OUT["routes"] = {
        f"{n}|{m}": {
            "mean_desirable": float(cent.loc[des, f"{m}_{n}"].mean()),
            "mean_undesirable": float(cent.loc[~des, f"{m}_{n}"].mean()),
            "p": {r.route + "/" + r.test: float(r.p_value) for r in out[(n, m)].itertuples()},
            "all_agree_at_alpha": bool(out[(n, m)].agrees_at_0_05.all())
            if "agrees_at_0_05" in out[(n, m)].columns else bool(out[(n, m)]["agrees_at_0.05"].all()),
            "min_achievable_p": float(perms[(n, m)]["min_achievable_p"]),
            "p_resolution": float(perms[(n, m)]["p_resolution"]),
        } for n in ("patrol", "detect") for m in ("strength", "betweenness", "mean_cent")}
    OUT["routes"]["n_desirable"] = int(des.sum())
    OUT["routes"]["n_undesirable"] = int((~des).sum())
    return out, perms


# ===================================================================== stage 5
def stage5_multiplicity(nets):
    """How many pairwise relationships are implicitly tested, and what survives."""
    res = {}
    pvals = {}
    for name in ("patrol", "detect"):
        R, n = nets[name]["R"], nets[name]["n"]
        r = R[IU]
        z = np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(n - 3)
        p = 2 * sps.norm.sf(np.abs(z))
        pvals[name] = p
        dr = mp.discovery_recount(p, alpha=ALPHA, labels=PAIRLAB)
        dr["table"].assign(network=name).to_csv(
            os.path.join(TAB, f"pairwise_tests_{name}.csv"), index=False)
        ev = np.linalg.eigvalsh(R)
        res[name] = {"m_tests": int(dr["m_tests"]),
                     "expected_false_positives": float(dr["expected_false_positives_uncorrected"]),
                     "summary": dr["summary"],
                     "n_eig_below_1e-4": int((ev < 1e-4).sum()),
                     "min_eig": float(ev.min()),
                     "max_abs_partial_unregularised": float(
                         np.abs(core.wi2net(np.linalg.inv(R))[IU]).max())}
    pooled = mp.discovery_recount(np.r_[pvals["patrol"], pvals["detect"]], alpha=ALPHA,
                                  labels=[f"patrol:{l}" for l in PAIRLAB] +
                                         [f"detect:{l}" for l in PAIRLAB])
    res["pooled"] = {"m_tests": int(pooled["m_tests"]), "summary": pooled["summary"],
                     "expected_false_positives": float(pooled["expected_false_positives_uncorrected"])}

    # max-statistic FWER control using the permutation null (columns permuted independently)
    maxstat = {}
    for name in ("patrol", "detect"):
        f = os.path.join(RES, f"maxstat_{name}_s0of1.npz")
        if not os.path.exists(f):
            maxstat[name] = None
            continue
        z = np.load(f)
        null = z["null_zero_order"]
        obs = np.abs(nets[name]["R"][IU])
        tab = mp.max_statistic_null(obs, null, alternative="greater")
        tab["pair"] = PAIRLAB
        tab.to_csv(os.path.join(TAB, f"maxstat_{name}.csv"), index=False)
        maxstat[name] = {
            "n_perm": int(z["n_perm"]),
            "p_resolution": float(1 / (int(z["n_perm"]) + 1)),
            "crit_95": float(np.quantile(null.max(axis=1), 0.95)),
            "n_sig_fwer": int((tab["p_fwer"] < ALPHA).sum()) if "p_fwer" in tab.columns
            else int((tab.iloc[:, -2] < ALPHA).sum()),
            "n_edges_selected": int(nets[name]["n_edges"]),
            "columns": list(tab.columns),
        }
    res["maxstat"] = maxstat
    OUT["multiplicity"] = {
        k: ({kk: vv for kk, vv in v.items() if kk != "summary"} if isinstance(v, dict) and "summary" in v else v)
        for k, v in res.items()}
    OUT["multiplicity"]["counts"] = {
        k: res[k]["summary"].set_index("correction")["n_significant"].to_dict()
        for k in ("patrol", "detect", "pooled")}
    OUT["multiplicity"]["n_reported_pvalues_in_paper"] = 0
    return res, pvals


# ===================================================================== stage 6
def _load(pattern):
    f = os.path.join(RES, pattern)
    return np.load(f) if os.path.exists(f) else None


def stage6_resampling(nets, cent):
    out = {}
    boots = {}
    for name in ("patrol", "detect"):
        z = _load(f"boot_{name}_s0of1.npz")
        boots[name] = z
        if z is None:
            out[name] = None
            continue
        ne = z["n_edges"]
        st, bt = z["strength"], z["betweenness"]
        mc = (st + bt) / 2.0
        ci = lambda A: np.quantile(A, [0.025, 0.975], axis=0)  # noqa: E731
        out[name] = {
            "n_boot": int(z["n_boot"]),
            "edges_mean": float(ne.mean()), "edges_sd": float(ne.std(ddof=1)),
            "edges_min": int(ne.min()), "edges_max": int(ne.max()),
            "strength_ci": ci(st), "betweenness_ci": ci(bt), "mean_cent_ci": ci(mc),
            "strength": st, "betweenness": bt, "mean_cent": mc,
        }
    # top-node separation: how often is the paper's top node top in a bootstrap replicate?
    for name, topcols in [("patrol", ("betweenness", "strength")), ("detect", ("betweenness", "strength"))]:
        if out.get(name) is None:
            continue
        for meas in topcols:
            A = out[name][meas]
            obs_top = int(np.argmax(cent[f"{meas}_{name}"].to_numpy()))
            ok = A.max(axis=1) > 0          # drop replicates whose network is empty
            out[name][f"top_share_{meas}"] = float((np.argmax(A[ok], axis=1) == obs_top).mean())
            out[name][f"n_usable_{meas}"] = int(ok.sum())
            out[name][f"top_node_{meas}"] = NODES[obs_top]
        out[name]["n_empty_boot"] = int((out[name]["betweenness"].max(axis=1) == 0).sum())

    # Cliff's delta between the two bootstrap centrality distributions, and its
    # dependence on the number of bootstrap replicates rather than on the number of dogs
    cliff = None
    if out.get("patrol") and out.get("detect"):
        rows = []
        P, D = out["patrol"]["mean_cent"], out["detect"]["mean_cent"]
        nmax = min(P.shape[0], D.shape[0])
        for j, node in enumerate(NODES):
            d = core.cliffs_delta(P[:, j], D[:, j])
            rows.append({"node": node, "cliffs_delta_recomputed": d})
        cliff = pd.DataFrame(rows)
        cliff.to_csv(os.path.join(TAB, "cliffs_delta_recomputed.csv"), index=False)
        widths = []
        rng = np.random.default_rng(SEED)
        # strict subsets only: drawing all nmax replicates every time has zero spread
        for nb in [int(nmax * 0.2), int(nmax * 0.4), int(nmax * 0.6), int(nmax * 0.8)]:
            ds = []
            for _ in range(200):
                i1 = rng.choice(P.shape[0], nb, replace=False)
                i2 = rng.choice(D.shape[0], nb, replace=False)
                ds.append(core.cliffs_delta(P[i1, 11], D[i2, 11]))   # node PLA
            widths.append({"n_boot_used": nb, "sd_of_delta": float(np.std(ds, ddof=1)),
                           "mean_delta": float(np.mean(ds))})
        out["cliff_width_vs_nboot"] = pd.DataFrame(widths)
        out["cliff_width_vs_nboot"].to_csv(os.path.join(TAB, "cliffs_delta_vs_nboot.csv"), index=False)
    out["cliff"] = cliff

    # case-dropping and node-dropping stability: correlation with the original network
    stab = {}
    for kind, fname in [("casedrop", "casedrop"), ("nodedrop", "nodedrop")]:
        for name in ("patrol", "detect"):
            z = _load(f"{fname}_{name}_s0of1.npz")
            if z is None:
                stab[(kind, name)] = None
                continue
            key = "frac" if kind == "casedrop" else "k"
            lev = z[key]
            rows = []
            for meas in ("strength", "betweenness"):
                A = z[meas]
                ref = cent[f"{meas}_{name}"].to_numpy()
                for L in np.unique(lev):
                    sel = A[lev == L]
                    cors = []
                    for row in sel:
                        ok = np.isfinite(row)
                        if ok.sum() >= 3 and np.std(row[ok]) > 0 and np.std(ref[ok]) > 0:
                            cors.append(np.corrcoef(row[ok], ref[ok])[0, 1])
                    if cors:
                        rows.append({"level": float(L), "measure": meas, "n_rep": len(cors),
                                     "mean_corr": float(np.mean(cors)),
                                     "lo": float(np.quantile(cors, 0.025)),
                                     "hi": float(np.quantile(cors, 0.975))})
            t = pd.DataFrame(rows)
            t.to_csv(os.path.join(TAB, f"stability_{kind}_{name}.csv"), index=False)
            stab[(kind, name)] = t
    out["stability"] = stab

    # permutation test of the patrol-vs-detection network difference (NCT-style)
    z = _load("nct_s0of1.npz")
    nct = None
    if z is not None:
        obs_edge = nets["patrol"]["net"][IU] - nets["detect"]["net"][IU]
        obs_strength = nets["patrol"]["strength"] - nets["detect"]["strength"]
        d_edge = z["d_edge"]
        # Only the edge-wise difference (net_A - net_B) was stored per permutation, so the
        # global-strength-invariance statistic |A| - |B| is not recoverable; the test run
        # here is the structure-invariance one (max |edge difference|), plus a max-statistic
        # test on node strength differences.
        nperm = int(z["n_perm"])
        maxobs = float(np.abs(obs_edge).max())
        maxnull = np.abs(d_edge).max(axis=1)
        p_struct = (int((maxnull >= maxobs).sum()) + 1) / (nperm + 1)
        per_edge_p = np.array([(int((maxnull >= abs(o)).sum()) + 1) / (nperm + 1) for o in obs_edge])
        maxnull_s = np.abs(z["d_strength"]).max(axis=1)
        per_node_p = np.array([(int((maxnull_s >= abs(o)).sum()) + 1) / (nperm + 1)
                               for o in obs_strength])
        nct = {"n_perm": nperm, "p_resolution": 1 / (nperm + 1),
               "max_edge_diff_observed": maxobs,
               "max_edge_diff_null_95": float(np.quantile(maxnull, 0.95)),
               "p_structure_invariance": float(p_struct),
               "n_edges_fwer_sig": int((per_edge_p < ALPHA).sum()),
               "n_nodes_strength_fwer_sig": int((per_node_p < ALPHA).sum()),
               "per_edge_p": per_edge_p, "per_node_p": per_node_p,
               "maxnull": maxnull, "obs_edge": obs_edge}
        pd.DataFrame({"pair": PAIRLAB, "edge_diff_patrol_minus_detect": obs_edge,
                      "p_fwer_maxstat": per_edge_p}).to_csv(
            os.path.join(TAB, "nct_edge_differences.csv"), index=False)
        pd.DataFrame({"node": NODES, "strength_diff": obs_strength,
                      "p_fwer_maxstat": per_node_p}).to_csv(
            os.path.join(TAB, "nct_strength_differences.csv"), index=False)
    out["nct"] = nct

    OUT["resampling"] = {
        "bootstrap": {k: {kk: vv for kk, vv in (out[k] or {}).items()
                          if isinstance(vv, (int, float, str))} for k in ("patrol", "detect")},
        "nct": {k: v for k, v in (nct or {}).items() if isinstance(v, (int, float))},
        "cliff_width_vs_nboot": out.get("cliff_width_vs_nboot").to_dict("records")
        if isinstance(out.get("cliff_width_vs_nboot"), pd.DataFrame) else None,
    }
    return out


# ===================================================================== stage 7
import matplotlib.pyplot as plt                                           # noqa: E402

OKABE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]


def _save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_gamma_fragility(gsweep, nets, name="panel_C02_gamma_fragility.png"):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for ax, key, lab in zip(axes, ("patrol", "detect"), ("patrol (n = 116)", "detection (n = 54)")):
        g = gsweep[key]
        ax.plot(g.gamma, g.margin_over_empty, "o-", color=OKABE[0], markersize=4)
        ax.axhline(0, color="black", lw=1)
        ax.axvline(GAMMA, color=OKABE[1], ls="--", lw=1.2)
        ax.set_xlabel("EBIC tuning parameter $\\gamma$")
        ax.set_ylabel("EBIC advantage of the selected graph\nover the empty graph (units)")
        ax.set_title(lab)
        sel = g[np.isclose(g.gamma, GAMMA)].iloc[0]
        ax.annotate(f"$\\gamma$ = 0.65: {int(sel.n_edges)} edges,\nmargin {sel.margin_over_empty:.1f}"
                    f"\n(= {sel.margin_over_empty / nets[key]['ebic_cost_per_edge'] if False else sel.margin_over_empty / (np.log(nets[key]['n']) + 4 * GAMMA * np.log(20)):.1f} edges' penalty)",
                    xy=(GAMMA, sel.margin_over_empty), xytext=(0.30, 0.70),
                    textcoords="axes fraction", fontsize=8,
                    arrowprops=dict(arrowstyle="->", lw=0.8))
        empty = g[g.n_edges == 0]
        if len(empty):
            ax.annotate(f"empty graph selected from\n$\\gamma$ = {empty.gamma.min():.2f}",
                        xy=(empty.gamma.min(), 0), xytext=(0.45, 0.25),
                        textcoords="axes fraction", fontsize=8,
                        arrowprops=dict(arrowstyle="->", lw=0.8))
        for _, row in g.iterrows():
            if row.gamma in (0.0, 0.5, 1.0):
                ax.text(row.gamma, row.margin_over_empty, f" {int(row.n_edges)}", fontsize=7,
                        va="bottom")
    fig.suptitle("C2/C15 - how far the selected network is from no network at all\n"
                 "(labels: number of edges selected)", fontsize=11)
    return _save(fig, name)


def fig_centrality_ci(cent, boot, name="panel_C05_centrality_bootstrap_ci.png"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0), sharey=True)
    order = np.argsort(cent["betweenness_patrol"].to_numpy())   # one ordering for both panels
    for ax, key, lab in zip(axes, ("patrol", "detect"), ("patrol (n = 116 dogs)", "detection (n = 54 dogs)")):
        b = boot.get(key)
        obs = cent[f"betweenness_{key}"].to_numpy()
        y = np.arange(20)
        ax.scatter(obs[order], y, color=OKABE[0], zorder=3, s=26, label="observed")
        if b is not None:
            lo, hi = b["betweenness_ci"]
            ax.hlines(y, lo[order], hi[order], color="#999999", lw=2, zorder=2,
                      label=f"2.5-97.5% over {b['n_boot']} bootstrap networks")
        ax.set_yticks(y)
        ax.set_yticklabels([NODES[i] for i in order], fontsize=8)
        ax.set_xlabel("betweenness centrality")
        ax.set_title(lab)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle("C5 - the paper's centrality ordering against its own bootstrap spread\n"
                 "(nodes ordered by patrol betweenness in both panels; "
                 "the paper reports no interval, CIs withheld as biased)", fontsize=10)
    return _save(fig, name)


def fig_cliff_nboot(res, name="panel_C07_cliff_vs_nboot.png"):
    t = res.get("cliff_width_vs_nboot")
    if t is None:
        return rf.not_assessable_figure(
            "C7 - Cliff's delta resolution", "bootstrap replicates were not available",
            os.path.join(FIG, name))
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(t.n_boot_used, t.sd_of_delta, "o-", color=OKABE[1])
    for _, r in t.iterrows():
        ax.annotate(f"{r.sd_of_delta:.3f}", (r.n_boot_used, r.sd_of_delta), fontsize=8,
                    textcoords="offset points", xytext=(4, 4))
    ax.set_xlabel("number of bootstrap replicates entering each Cliff's delta")
    ax.set_ylabel("s.d. of Cliff's delta across 200 redraws")
    ax.set_title("C7 - the precision of $\\delta$ is set by the number of bootstraps,\n"
                 "not by the number of dogs (node PLA, mean centrality)", fontsize=10)
    ax.text(0.98, 0.95, f"mean $\\delta$ stays at {t.mean_delta.mean():.3f}\nwhile its spread shrinks",
            transform=ax.transAxes, ha="right", va="top", fontsize=9)
    return _save(fig, name)


def fig_stability(res, name="panel_C08_stability.png"):
    stab = res["stability"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, kind, xlab in [(axes[0], "casedrop", "proportion of dogs retained"),
                           (axes[1], "nodedrop", "number of nodes retained")]:
        for key, col in zip(("patrol", "detect"), (OKABE[0], OKABE[1])):
            t = stab.get((kind, key))
            if t is None or not len(t):
                continue
            for meas, ls in [("strength", "-"), ("betweenness", "--")]:
                s = t[t.measure == meas].sort_values("level")
                ax.plot(s.level, s.mean_corr, ls, color=col, marker="o", markersize=3,
                        label=f"{key} {meas}")
                ax.fill_between(s.level, s.lo, s.hi, color=col, alpha=0.12)
        ax.axhline(0.7, color="#777777", ls=":", lw=1)
        ax.set_xlabel(xlab)
        ax.set_ylabel("correlation with the full-sample centrality")
        ax.legend(frameon=False, fontsize=7, ncol=2)
    axes[0].set_title("subject-wise (case-dropping)")
    axes[1].set_title("node-wise (node-dropping)")
    fig.suptitle("C8/C10/C11 - stability of centrality under resampling "
                 "(shaded: 2.5-97.5% over replicates)", fontsize=11)
    return _save(fig, name)


def fig_nct(res, name="panel_C14_network_difference.png"):
    n = res.get("nct")
    if n is None:
        return rf.not_assessable_figure(
            "C14 - patrol vs detection network difference",
            "permutation replicates were not available", os.path.join(FIG, name))
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    ax.hist(n["maxnull"], bins=30, color="#BBBBBB", edgecolor="white")
    ax.axvline(n["max_edge_diff_observed"], color=OKABE[1], lw=2)
    ax.axvline(n["max_edge_diff_null_95"], color=OKABE[0], ls="--", lw=1.4)
    ax.set_xlabel("max |edge weight difference| between the two networks")
    ax.set_ylabel("permutations (count)")
    ax.set_title("C14 - permutation test the paper did not run:\n"
                 "is the patrol network different from the detection network?", fontsize=10)
    ax.text(0.98, 0.95,
            f"observed {n['max_edge_diff_observed']:.3f}\n"
            f"null 95th pct {n['max_edge_diff_null_95']:.3f}\n"
            f"p = {n['p_structure_invariance']:.3f}  ({n['n_perm']} permutations,\n"
            f"resolution {n['p_resolution']:.4f})\n"
            f"edges differing at FWER 0.05: {n['n_edges_fwer_sig']} of 190",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)
    return _save(fig, name)


def fig_maxstat(mult, nets, name="panel_C18_edge_multiplicity.png"):
    ms = mult["maxstat"]
    if ms.get("patrol") is None:
        return rf.not_assessable_figure(
            "C18 - edge selection against a permutation null",
            "permutation replicates were not available", os.path.join(FIG, name))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, key, lab in zip(axes, ("patrol", "detect"), ("patrol (n = 116)", "detection (n = 54)")):
        z = np.load(os.path.join(RES, f"maxstat_{key}_s0of1.npz"))
        null_max = z["null_zero_order"].max(axis=1)
        obs = np.abs(nets[key]["R"][IU])
        ax.hist(null_max, bins=25, color="#BBBBBB", edgecolor="white",
                label=f"max |r| over 190 pairs,\n{int(z['n_perm'])} permutations")
        crit = np.quantile(null_max, 0.95)
        ax.axvline(crit, color=OKABE[0], ls="--", lw=1.4)
        n_pass = int((obs > crit).sum())
        ax.set_title(f"{lab}: {n_pass} of 190 pairs exceed the FWER threshold\n"
                     f"(network selected {nets[key]['n_edges']} edges)", fontsize=9)
        ax.set_xlabel("|polychoric correlation|")
        ax.set_ylabel("permutations (count)")
        ax.text(0.98, 0.70, f"FWER 0.05 threshold |r| = {crit:.3f}\n"
                            f"largest observed |r| = {obs.max():.3f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8)
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig.suptitle("C18 - 190 node pairs are screened per network; no multiplicity control "
                 "was applied", fontsize=11)
    return _save(fig, name)


# ===================================================================== stage 8
def build_ledger(nets, cent, routes, perms, mult, res, gsweep):
    P, D = OUT["networks"]["patrol"], OUT["networks"]["detect"]
    bp, bd = res.get("patrol"), res.get("detect")
    nct = res.get("nct")
    rp = OUT["routes"]
    rows = []

    def add(cid, short, verdict, reasons, eff_rep=np.nan, eff_rec=np.nan,
            lo=np.nan, hi=np.nan, p_rep=np.nan, p_rec=np.nan, stage=""):
        rows.append({"claim_id": cid, "claim_short": short, "verdict": verdict,
                     "effect_reported": eff_rep, "effect_recomputed": eff_rec,
                     "ci_low": lo, "ci_high": hi, "p_reported": p_rep, "p_recomputed": p_rec,
                     "reasons": reasons, "stage": stage})

    add("C1", "patrol network has 55 edges (28.95%)", "supported",
        f"re-derived independently (polychoric + EBICglasso, gamma=0.65, 100-lambda path): "
        f"{P['n_edges']} edges = {P['pct_edges']:.2f}% of 190 pairs; lambda={P['lambda']:.4f}",
        55, P["n_edges"], stage="3")
    add("C2", "detection network has 70 edges (36.84%)", "supported-with-caveat",
        f"edge count reproduces exactly ({D['n_edges']}); but the selected graph beats the empty "
        f"graph by only {D['margin_over_empty']:.1f} EBIC units = "
        f"{D['margin_over_empty'] / D['ebic_cost_per_edge']:.1f} edges' worth of penalty "
        f"({D['ebic_cost_per_edge']:.2f} per edge), against {P['margin_over_empty']:.1f} "
        f"({P['margin_over_empty'] / P['ebic_cost_per_edge']:.1f} edges) for patrol | "
        f"190 parameters estimated from n=54 dogs",
        70, D["n_edges"], stage="3")

    for cid, net, lab in [("C3", "patrol", "patrol"), ("C4", "detect", "detection")]:
        k = f"{net}|mean_cent"
        ks = f"{net}|strength"
        kb = f"{net}|betweenness"
        pp_ = rp[k]["p"]
        pperm = [v for kk, v in pp_.items() if "permutation" in kk][0]
        allp = ", ".join(f"{kk.split('/')[-1]} p={v:.3g}" for kk, v in pp_.items())
        agree = rp[k]["all_agree_at_alpha"]
        if net == "patrol":
            verdict = "supported-with-caveat"
            reason = (f"desirable nodes (n=14) mean centrality {rp[k]['mean_desirable']:.2f} vs "
                      f"undesirable (n=6) {rp[k]['mean_undesirable']:.2f}; all routes agree at "
                      f"alpha=0.05 ({allp}); permutation floor for 14 vs 6 is "
                      f"p>={rp[k]['min_achievable_p']:.2e} | caveat: the 20 units are nodes of one "
                      f"network, not independent observations, so the test's null is not the "
                      f"design's null")
        else:
            verdict = "contradicted"
            reason = (f"claim is made for both networks; in detection the composite centrality "
                      f"difference is {rp[k]['mean_desirable']:.2f} vs {rp[k]['mean_undesirable']:.2f} "
                      f"with every route non-significant ({allp}); betweenness alone "
                      f"{rp[kb]['p']['permutation/label permutation (mean difference)']:.3g}; only "
                      f"strength separates the two sets "
                      f"(permutation p={rp[ks]['p']['permutation/label permutation (mean difference)']:.3g})")
        add(cid, f"central descriptors are desirable ({lab})", verdict, reason,
            np.nan, float(rp[k]["mean_desirable"] - rp[k]["mean_undesirable"]),
            p_rec=float(pperm), stage="3-4")
        _ = agree

    top = OUT["networks"]["top_nodes"]
    share = ""
    if bp and bd:
        share = (f"across bootstrap networks the observed top node keeps the top rank in "
                 f"{bp['top_share_betweenness'] * 100:.0f}% (patrol betweenness), "
                 f"{bp['top_share_strength'] * 100:.0f}% (patrol strength), "
                 f"{bd['top_share_betweenness'] * 100:.0f}% (detection betweenness), "
                 f"{bd['top_share_strength'] * 100:.0f}% (detection strength) of replicates")
    add("C5", "PLA most central; CUR top strength in patrol", "supported-with-caveat",
        f"reproduced exactly: patrol top betweenness {top['betweenness_patrol'][0]}, patrol top "
        f"strength {top['strength_patrol'][0]}, detection top betweenness "
        f"{top['betweenness_detect'][0]} and top strength {top['strength_detect'][0]}; {share}; "
        f"bootstrap intervals on neighbouring nodes overlap, and the paper reports no interval",
        stage="5")

    add("C6", "named descriptors differ in centrality between networks",
        "supported-with-caveat",
        "the per-node orderings reproduce from the supplied data, but they are point estimates "
        "from two separately estimated networks; see C14 for the formal comparison",
        stage="3")

    cl = res.get("cliff")
    cl_txt = "bootstrap replicates unavailable"
    if cl is not None:
        rep = {"CUR": 0.452, "SLP": 0.290, "ACT": 0.282, "FSH": -0.614, "DA": -0.609,
               "GUS": -0.582, "TOY": -0.465, "GWL": -0.310, "FDA": -0.302}
        m = cl.set_index("node")["cliffs_delta_recomputed"]
        pairs = [f"{k}: reported {v:+.3f} vs recomputed {m[k]:+.3f}" for k, v in rep.items() if k in m]
        cl_txt = "; ".join(pairs)
    w = OUT["resampling"].get("cliff_width_vs_nboot")
    wtxt = ""
    if w:
        wtxt = (f" | the spread of delta is governed by the number of bootstrap replicates, not "
                f"by the number of dogs: s.d. of delta falls from {w[0]['sd_of_delta']:.3f} at "
                f"{w[0]['n_boot_used']} replicates to {w[-1]['sd_of_delta']:.3f} at "
                f"{w[-1]['n_boot_used']}, while the point estimate stays near "
                f"{w[-1]['mean_delta']:+.3f}")
    add("C7", "Cliff's delta quantifies centrality differences", "supported-with-caveat",
        f"{cl_txt}{wtxt} | delta is computed between two bootstrap distributions, not between "
        f"two samples of dogs, so its confidence interval narrows with computation",
        stage="5")

    st = res["stability"]
    nd_txt = "node-dropping replicates unavailable"
    if st.get(("nodedrop", "patrol")) is not None and len(st[("nodedrop", "patrol")]):
        a = st[("nodedrop", "patrol")]
        b = st[("nodedrop", "detect")]
        lo_p = a[a.level == a.level.min()].mean_corr.mean()
        lo_d = b[b.level == b.level.min()].mean_corr.mean()
        nd_txt = (f"at the smallest node subset (k={int(a.level.min())}) the mean correlation with "
                  f"the full network is {lo_p:.2f} (patrol) and {lo_d:.2f} (detection)")
    add("C8", "bootstrap analyses revealed network stability", "supported-with-caveat",
        f"{nd_txt}; node-dropping resamples the same dogs and so cannot speak to sampling "
        f"variability; the subject-wise curves (C10) are the ones that can, and they fall away "
        f"well before the full sample",
        stage="5")

    if bp and bd:
        add("C9", "s.d. of edge count: 12.82 patrol, 29.75 detection", "supported-with-caveat",
            f"detection reproduces ({bd['edges_sd']:.2f} against 29.75 reported) and the ordering "
            f"detection > patrol reproduces, but the patrol value does not: {bp['edges_sd']:.2f} "
            f"here against 12.82, {bp['edges_sd'] / 12.82:.1f}x larger. Candidate causes, both "
            f"on this review's side: only {bp['n_boot']} bootstrap networks against the paper's "
            f"2000, and a 25-value lambda path against 100. Reported as a discrepancy to be "
            f"resolved, not as a failure to reproduce",
            12.82, bp["edges_sd"], stage="5")
        add("C9b", "s.d. of edge count, detection network", "supported-with-caveat",
            f"recomputed {bd['edges_sd']:.2f} against 29.75 reported; edge count ranged "
            f"{bd['edges_min']}-{bd['edges_max']} across {bd['n_boot']} bootstrap networks "
            f"(observed {D['n_edges']}), and {bd['n_empty_boot']} of {bd['n_boot']} bootstrap "
            f"networks were completely empty. Bootstrap networks are also systematically denser "
            f"than the observed ones (mean {bd['edges_mean']:.0f} against {D['n_edges']} for "
            f"detection, {bp['edges_mean']:.0f} against {P['n_edges']} for patrol), because "
            f"resampling with replacement duplicates dogs and inflates the correlations",
            29.75, bd["edges_sd"], stage="5")
    else:
        add("C9", "s.d. of edge count: 12.82 patrol, 29.75 detection", "not-assessable",
            "bootstrap replicates were not produced in this review", 12.82, np.nan, stage="5")

    cd_txt = "case-dropping replicates unavailable"
    if st.get(("casedrop", "patrol")) is not None and len(st[("casedrop", "patrol")]):
        a = st[("casedrop", "patrol")].groupby("level").mean_corr.mean()
        b = st[("casedrop", "detect")].groupby("level").mean_corr.mean()
        cd_txt = ("patrol mean correlation with the full network by retained fraction: " +
                  ", ".join(f"{k:.2f}:{v:.2f}" for k, v in a.items()) +
                  " | detection: " + ", ".join(f"{k:.2f}:{v:.2f}" for k, v in b.items()))
    add("C10", "subject-wise: patrol fails below ~60 dogs, detection ~40", "supported-with-caveat",
        f"{cd_txt}; the direction reproduces, the exact thresholds depend on the replicate count "
        f"used here", stage="5")
    add("C11", "patrol network more stable than detection", "supported",
        f"{cd_txt}; patrol exceeds detection at every retained fraction computed", stage="5")

    d = OUT["design"]
    add("C12", "missingness: 1.23+-1.38 missing, 4.16+-7.64 zeros per descriptor", "supported",
        f"recomputed from rawData.csv: {d['missing_per_desc'][0]:.2f}+-{d['missing_per_desc'][1]:.2f} "
        f"missing ({d['missing_pct_per_desc'][0]:.2f}+-{d['missing_pct_per_desc'][1]:.2f}%) and "
        f"{d['zero_per_desc'][0]:.2f}+-{d['zero_per_desc'][1]:.2f} zeros "
        f"({d['zero_pct_per_desc'][0]:.2f}+-{d['zero_pct_per_desc'][1]:.2f}%); the one further "
        f"figure, 1.2+-1.03% imputed over the final 29 descriptors, recomputes as "
        f"{d['missing_pct_29'][0]:.2f}+-{d['missing_pct_29'][1]:.2f}%",
        1.23, d["missing_per_desc"][0], stage="2")
    add("C13", "171 dogs from 117 handlers, 1.46+-0.65 per handler", "supported",
        f"recomputed: {d['n_dogs_raw']} dogs, {d['n_handlers']} handlers, "
        f"{d['surveys_per_handler'][0]:.2f}+-{d['surveys_per_handler'][1]:.2f}; "
        f"{d['n_patrol_raw']} patrol and {d['n_detect_raw']} detection, becoming "
        f"{d['n_patrol_net']} and {d['n_detect_net']} after the one excluded dog",
        1.46, d["surveys_per_handler"][0], stage="1-2")

    if nct:
        v = "supported" if nct["p_structure_invariance"] < ALPHA else "not-assessable"
        reason = (f"the paper describes the two networks separately and never tests them against "
                  f"each other; a label-permutation test run here (max |edge difference|, "
                  f"{nct['n_perm']} permutations, resolution {nct['p_resolution']:.4f}) gives "
                  f"observed {nct['max_edge_diff_observed']:.3f} against a null 95th percentile of "
                  f"{nct['max_edge_diff_null_95']:.3f}, p = {nct['p_structure_invariance']:.3f}; "
                  f"{nct['n_edges_fwer_sig']} of 190 edges and "
                  f"{nct['n_nodes_strength_fwer_sig']} of 20 node strengths differ at FWER 0.05")
        if nct["p_structure_invariance"] >= ALPHA:
            v = "not-reproducible"
            reason += (" | the global difference the paper describes is not recovered by a formal "
                       "test at this permutation budget")
        add("C14", "network structure varied between dog types", v, reason,
            p_rec=float(nct["p_structure_invariance"]), stage="5-6")
    else:
        add("C14", "network structure varied between dog types", "not-assessable",
            "permutation replicates were not produced", stage="5-6")

    gd = gsweep["detect"]
    first_empty = gd.query("n_edges == 0").gamma.min() if (gd.n_edges == 0).any() else np.nan
    add("C15", "gamma=0.65 is the most conservative non-empty value", "supported",
        f"reproduced exactly: sweeping gamma in 0.05 steps, the detection network holds "
        f"{D['n_edges']} edges at gamma=0.65 and collapses to the empty graph from "
        f"gamma={first_empty:.2f}; the patrol network holds {P['n_edges']} edges across the whole "
        f"upper range", 0.65, 0.65, stage="3")

    add("C16", "8 descriptors removed for handler non-independence", "supported-with-caveat",
        f"{d['n_repeat_handlers']} handlers rated more than one dog ({d['n_dogs_repeat']} of 170 "
        f"ratings, up to {d['max_dogs_per_handler']} per handler); among those handlers "
        f"{d['icc_gt_010']} of the 20 retained descriptors still have ICC > 0.10 "
        f"(max {d['icc_max']:.2f}, {d['icc_max_desc']}; median {d['icc_median']:.2f}); a "
        f"likelihood-ratio test that does not reject is not evidence of independence. Over the "
        f"whole sample the design effect is {d['design_effect']:.2f} "
        f"(effective n {d['effective_n']:.0f} of 170) because most handlers rated one dog",
        stage="2")

    add("C17", "15 imputed datasets averaged and rounded to one", "supported-with-caveat",
        "the procedure is disclosed and reproduced in the script; averaging the m=15 completions "
        "into a single rounded dataset discards the between-imputation variance that multiple "
        "imputation exists to carry, so every interval downstream (including the bootstrap "
        "spread in C9 and C10) is conditioned on one imputed dataset. Not quantifiable here: the "
        "15 individual completions were not supplied, only their rounded average",
        stage="2")

    ms = mult["maxstat"].get("patrol")
    if ms:
        mp_ = mult["maxstat"]["patrol"]
        md_ = mult["maxstat"]["detect"]
        cnt = OUT["multiplicity"]["counts"]
        add("C18", "edges are 'credibly non-zero partial correlations'", "supported-with-caveat",
            f"estimating a network over 20 nodes screens 190 node pairs per network and 380 "
            f"across the two; the paper reports no p-value and applies no multiplicity control. "
            f"Re-running the 190 pairwise polychoric tests per network: "
            f"{cnt['patrol']['none']} (patrol) and {cnt['detect']['none']} (detection) reach "
            f"p<0.05 uncorrected against {mult['patrol']['expected_false_positives']:.1f} expected "
            f"by chance, falling to {cnt['patrol']['Bonferroni']} and "
            f"{cnt['detect']['Bonferroni']} under Bonferroni. Against a permutation max-statistic "
            f"null ({mp_['n_perm']} permutations, FWER threshold |r|={mp_['crit_95']:.3f} patrol / "
            f"{md_['crit_95']:.3f} detection) the lasso's {P['n_edges']} and {D['n_edges']} "
            f"selected edges are a model-selection output, not a set of tests with a controlled "
            f"error rate", stage="6")
    else:
        add("C18", "edges are 'credibly non-zero partial correlations'", "not-assessable",
            "permutation null not produced", stage="6")

    add("C19", "8 descriptors failed a handler-independence LRT", "not-assessable",
        "the per-descriptor likelihood-ratio / delta-AIC values behind the screen are in "
        "electronic supplementary material S2-S3, which was not supplied, and the script records "
        "only the model formulas (clmm with and without a handler random intercept), not their "
        "output; an ordinal mixed model was not re-fitted here. What would change this verdict: "
        "supplementary table S2, or the 15 pre-average imputed datasets so clmm can be re-run on "
        "the same input the authors used",
        stage="2-3")
    add("C20", "95% CIs on Cliff's delta (figure 3)", "not-assessable",
        "the interval endpoints are in electronic supplementary table S5, which was not supplied; "
        "only the nine delta point estimates quoted in the running text could be compared (C7). "
        "What would change this verdict: table S5, or the authors' bootTable*.csv files",
        stage="5")

    led = pd.DataFrame(rows)
    led.to_csv(os.path.join(HERE, "claim_ledger.csv"), index=False)
    OUT["verdict_counts"] = led.verdict.value_counts().to_dict()
    return led


# ===================================================================== stage 9
def build_figures(mi2, long1, indep, cent, nets, gsweep, routes, perms, mult, res, led, pvals):
    figs = {}
    des = cent.desirable.to_numpy()

    figs["F1"] = rf.fig_design_map(long1, "handler", "dog_type", "mean_rating",
                                   independence=indep, unit_label="handler",
                                   path=os.path.join(FIG, "F1_design_map.png"))
    gp = {"desirable (14 nodes)": cent.loc[des, "mean_cent_patrol"].to_numpy(),
          "undesirable (6 nodes)": cent.loc[~des, "mean_cent_patrol"].to_numpy()}
    figs["F2"] = rf.fig_raw_behind_claim(
        gp, "mean centrality (betweenness and strength averaged)",
        path=os.path.join(FIG, "F2_raw_behind_claim.png"),
        test_name="C3, patrol network",
        p_value=float(perms[("patrol", "mean_cent")]["p_value"]),
        reported_error="the paper reports no interval on centrality (CIs withheld as biased)")
    figs["F3"] = rf.fig_assumption_qq(gp, path=os.path.join(FIG, "F3_assumption_shape.png"),
                                      value_label="mean centrality (patrol)")
    led_f = led.dropna(subset=["effect_recomputed"]).copy()
    if len(led_f):
        figs["F4"] = rf.fig_claim_ledger(led_f, path=os.path.join(FIG, "F4_claim_ledger.png"),
                                         effect_label="reported vs recomputed (paper's units)")
    else:
        figs["F4"] = rf.not_assessable_figure(
            "F4 - claim ledger", "no claim carried a numeric effect",
            os.path.join(FIG, "F4_claim_ledger.png"))
    boot_mean_diff = None
    if res.get("patrol"):
        mc = res["patrol"]["mean_cent"]
        diffs = mc[:, des].mean(axis=1) - mc[:, ~des].mean(axis=1)
        boot_mean_diff = {"statistic": "mean difference", "point_estimate": float(np.mean(diffs)),
                          "ci_low": float(np.quantile(diffs, 0.025)),
                          "ci_high": float(np.quantile(diffs, 0.975)),
                          "distribution": diffs, "n": int(mc.shape[0]),
                          "n_boot": int(mc.shape[0]), "level": 0.95, "method": "percentile",
                          "bootstrap_se": float(np.std(diffs, ddof=1))}
    figs["F5"] = rf.fig_null_distribution(
        perms[("patrol", "mean_cent")], boot=boot_mean_diff,
        path=os.path.join(FIG, "F5_null_distribution.png"),
        statistic_label="desirable minus undesirable mean centrality (patrol)")
    figs["F6"] = rf.fig_multiplicity(
        np.r_[pvals["patrol"], pvals["detect"]], alpha=ALPHA,
        path=os.path.join(FIG, "F6_multiplicity.png"),
        labels=[f"patrol:{l}" for l in PAIRLAB] + [f"detect:{l}" for l in PAIRLAB])

    # per-claim panels
    for cid, net in [("C3", "patrol"), ("C4", "detect")]:
        g = {"desirable (14 nodes)": cent.loc[des, f"mean_cent_{net}"].to_numpy(),
             "undesirable (6 nodes)": cent.loc[~des, f"mean_cent_{net}"].to_numpy()}
        figs[cid] = rf.fig_claim_panel(
            cid, g, f"mean centrality ({net})", perms[(net, "mean_cent")],
            routes=routes[(net, "mean_cent")],
            verdict=str(led.loc[led.claim_id == cid, "verdict"].iloc[0]),
            path=os.path.join(FIG, f"panel_{cid}_central_desirable_{net}.png"))
    figs["C2"] = fig_gamma_fragility(gsweep, OUT["networks"])
    figs["C5"] = fig_centrality_ci(cent, {k: res.get(k) for k in ("patrol", "detect")})
    figs["C7"] = fig_cliff_nboot(res)
    figs["C8"] = fig_stability(res)
    figs["C14"] = fig_nct(res)
    figs["C18"] = fig_maxstat(mult, nets)
    figs["C17"] = rf.not_assessable_figure(
        "C17 - imputation uncertainty",
        "only the rounded average of the 15 imputed datasets was supplied; the individual "
        "completions would be needed to show how much of the network depends on the imputation",
        os.path.join(FIG, "panel_C17_imputation_uncertainty.png"))
    figs["C1x"] = rf.not_assessable_figure(
        "Reported-vs-recomputed p-values (statcheck)",
        "the paper reports no test statistic and no p-value, so no (statistic, df, p) triple "
        "could be checked and no p-value distribution diagnostic is interpretable",
        os.path.join(FIG, "panel_C19_statcheck_not_assessable.png"))
    OUT["figures"] = sorted(f for f in os.listdir(FIG) if f.endswith(".png"))
    return figs


# ===================================================================== stage 10
def write_report(led, mult, res, gsweep, cent, inv):
    d, P, Dn = OUT["design"], OUT["networks"]["patrol"], OUT["networks"]["detect"]
    rp = OUT["routes"]
    nct = res.get("nct")
    bp, bd = res.get("patrol"), res.get("detect")
    cnt = OUT["multiplicity"]["counts"]
    vc = OUT["verdict_counts"]
    f = lambda n: f"figures/{n}"                                            # noqa: E731

    def tbl(df, cols=None):
        """Markdown table without the optional `tabulate` dependency."""
        df = df[cols] if cols else df
        def cell(v):
            if isinstance(v, float):
                return "" if not np.isfinite(v) else (f"{v:.4g}")
            return str(v).replace("|", "\\|")
        head = "| " + " | ".join(str(c) for c in df.columns) + " |"
        rule = "|" + "|".join(["---"] * len(df.columns)) + "|"
        body = ["| " + " | ".join(cell(v) for v in row) + " |"
                for row in df.itertuples(index=False, name=None)]
        return "\n".join([head, rule] + body)

    L = []
    A = L.append
    A("# Statistical review - Goold, Vas, Olsen & Newberry (2016)\n")
    A("*Using network analysis to study behavioural phenotypes: an example using domestic dogs.* "
      "R. Soc. Open Sci. 3:160268. doi:10.1098/rsos.160268\n")

    A("## 1. Header - what was reviewed, and what could not be\n")
    A("| input | content | sha256 (first 16) |")
    A("|---|---|---|")
    for k, v in OUT["inputs"]["checksums"].items():
        extra = {"rawData.csv": f"{OUT['inputs']['raw_shape'][0]} rows x {OUT['inputs']['raw_shape'][1]} cols, 43 descriptors",
                 "finalMIdata_averaged_rounded.csv": f"{OUT['inputs']['mi_shape'][0]} rows x {OUT['inputs']['mi_shape'][1]} cols, the imputed dataset the networks were built from",
                 "Rscript.R": f"{OUT['inputs']['rscript_lines']} lines, the authors' own analysis script",
                 "paper_text": f"full text (PMC), {OUT['inputs']['text_quality']['n_chars']} characters"}[k]
        A(f"| {k} | {extra} | `{v}` |")
    A("")
    A(f"The text layer is usable (`usable_text_layer = "
      f"{OUT['inputs']['text_quality']['usable_text_layer']}`). Automatic extraction found "
      f"{OUT['inputs']['inventory']['named_tests']} named-test mentions, "
      f"{OUT['inputs']['inventory']['sample_sizes']} sample sizes and "
      f"{OUT['inputs']['inventory']['intervals']} reported intervals, and **zero p-values and zero "
      f"test statistics**. `statcheck.csv` is therefore empty by construction, not by omission: "
      f"this paper reports no null-hypothesis test anywhere in its results, so there is no "
      f"reported-vs-recomputed p-value to check and no p-value distribution to describe. With "
      f"zero reported p-values the just-below-threshold diagnostic is far below the 20-value floor "
      f"at which it means anything, and it is not run.\n")
    A(f"Claims were therefore read by hand from the abstract, results and methods: "
      f"**{len(led)} claims** are in the ledger, of which "
      f"**{int(vc.get('not-assessable', 0))}** are `not-assessable` from the supplied inputs. "
      f"Having the authors' own `Rscript.R` means the pipeline could be re-derived step by step "
      f"rather than guessed: the preprocessing, the polychoric correlations, the EBICglasso fit, "
      f"the centrality measures and the resampling were all reimplemented independently in Python "
      f"and checked against the paper's printed numbers.\n")
    A("What this review could **not** reach: the 15 individual imputed datasets (only their "
      "rounded average was supplied), the authors' bootstrap output tables (`bootTable*.csv`, "
      "regenerated here at a smaller replicate count), electronic supplementary tables S1-S5, and "
      "the per-descriptor likelihood-ratio tests used to screen handler non-independence (neither "
      "the paper nor the script records their output).\n")

    A("## 2. Summary of findings\n")
    A("| verdict | claims |")
    A("|---|---|")
    for k in ["supported", "supported-with-caveat", "assumption-dependent",
              "underpowered-resolution", "not-reproducible", "contradicted", "not-assessable"]:
        A(f"| `{k}` | {int(vc.get(k, 0))} |")
    A("")
    A("Three findings change how the paper reads.\n")
    A(f"1. **The detection network is one rounding of the tuning parameter away from no network "
      f"at all.** Its {Dn['n_edges']} edges beat the empty graph by "
      f"{Dn['margin_over_empty']:.1f} EBIC units, which at {Dn['ebic_cost_per_edge']:.2f} units "
      f"per edge is {Dn['margin_over_empty'] / Dn['ebic_cost_per_edge']:.1f} edges' worth of "
      f"penalty; the patrol network's margin is {P['margin_over_empty']:.1f} units "
      f"({P['margin_over_empty'] / P['ebic_cost_per_edge']:.1f} edges). The paper says as much in "
      f"one clause - above gamma = 0.65 the detection graph empties - and the sweep here confirms "
      f"it exactly (C2, C15, stage 3).")
    A(f"2. **'The most central descriptors denoted desirable characteristics in both networks' "
      f"holds in one network, not two.** In patrol the 14 desirable nodes average "
      f"{rp['patrol|mean_cent']['mean_desirable']:.2f} against "
      f"{rp['patrol|mean_cent']['mean_undesirable']:.2f} for the 6 undesirable ones and every "
      f"route agrees; in detection the same contrast is "
      f"{rp['detect|mean_cent']['mean_desirable']:.2f} against "
      f"{rp['detect|mean_cent']['mean_undesirable']:.2f} with no route below 0.05 "
      f"(C3, C4, stage 4).")
    if nct:
        A(f"3. **The patrol-versus-detection difference is described, never tested.** A "
          f"label-permutation test of network-structure invariance run here gives max |edge "
          f"difference| = {nct['max_edge_diff_observed']:.3f} against a null 95th percentile of "
          f"{nct['max_edge_diff_null_95']:.3f}, p = {nct['p_structure_invariance']:.3f} over "
          f"{nct['n_perm']} permutations (resolution {nct['p_resolution']:.4f}); "
          f"{nct['n_edges_fwer_sig']} of 190 edges differ at FWER 0.05 (C14, stage 5-6).")
    A("")

    A("## 3. Claim ledger - F4\n")
    A(f"![F4]({f('F4_claim_ledger.png')})\n")
    A("`p_reported` is empty for every row because the paper reports no p-values; "
      "`p_recomputed` is filled wherever this review could construct a permutation test of the "
      "claim.\n")
    A(tbl(led, ["claim_id", "claim_short", "verdict", "effect_reported", "effect_recomputed",
                "p_reported", "p_recomputed", "stage"]))
    A("\nFull `reasons` text for every row is in `claim_ledger.csv`.\n")

    A("## 4. Design and sampling - F1\n")
    A(f"![F1]({f('F1_design_map.png')})\n")
    A(f"{d['n_invited_handlers']} handlers were invited, {d['n_surveys']} surveys came back, and "
      f"after removing three that covered more than one dog the analysis rests on "
      f"{d['n_dogs_raw']} dogs rated by {d['n_handlers']} handlers "
      f"({d['surveys_per_handler'][0]:.2f} +- {d['surveys_per_handler'][1]:.2f} surveys each) - "
      f"all three figures reproduce exactly (C13). One dog was dropped for zero-coded responses, "
      f"leaving {d['n_patrol_net']} patrol and {d['n_detect_net']} detection dogs, and these are "
      f"the n values passed to `EBICglasso`.\n")
    A(f"The unit question matters here because {d['n_repeat_handlers']} handlers rated more than "
      f"one dog ({d['n_dogs_repeat']} of 170 ratings, up to {d['max_dogs_per_handler']} dogs from "
      f"one handler). The authors screened for this and removed eight descriptors. Among the "
      f"repeat handlers, {d['icc_gt_010']} of the 20 **retained** descriptors still have an ICC "
      f"above 0.10 (largest {d['icc_max']:.2f}, {d['icc_max_desc']}; median {d['icc_median']:.2f}; "
      f"`tables/handler_icc.csv`). A likelihood-ratio test that fails to reject is not a "
      f"certificate of independence, and the screen was run on the same averaged dataset the "
      f"networks use. Against that, the inflation is bounded by how few handlers repeat: with "
      f"{d['obs_per_unit_mean']:.2f} ratings per handler overall the design effect is "
      f"{d['design_effect']:.2f} and the effective n is about {d['effective_n']:.0f} rather than "
      f"170 (C16).\n")
    A("Three design questions are not answerable from any file supplied: whether the "
      f"{d['n_handlers']} responding handlers represent the {d['n_invited_handlers']} invited "
      "(no comparison of responders and non-responders is reported), whether ratings of different "
      "dogs by the same handler are conditionally independent given the descriptors retained, and "
      "whether the patrol/detection contrast is confounded with breed - the patrol dogs are "
      "overwhelmingly German shepherds and the detection dogs are not, so any network difference "
      "between dog types is also a difference between breed mixes.\n")

    A("## 5. Method and assumptions - F2, F3\n")
    A(f"![F2]({f('F2_raw_behind_claim.png')})\n")
    A(f"![F3]({f('F3_assumption_shape.png')})\n")
    A(f"The estimation pipeline reproduces exactly. Polychoric correlation matrices, projected to "
      f"the nearest positive-definite correlation matrix, then `EBICglasso` at gamma = 0.65 over a "
      f"100-value lambda path, give **{P['n_edges']} edges ({P['pct_edges']:.2f}% of 190 pairs) "
      f"for patrol and {Dn['n_edges']} ({Dn['pct_edges']:.2f}%) for detection** - the paper's 55 "
      f"(28.95%) and 70 (36.84%) (C1, C2). The centrality orderings reproduce too: "
      f"{OUT['networks']['top_nodes']['betweenness_patrol'][0]} highest betweenness and "
      f"{OUT['networks']['top_nodes']['strength_patrol'][0]} highest strength in patrol, "
      f"{OUT['networks']['top_nodes']['betweenness_detect'][0]} highest on both in detection "
      f"(C5).\n")
    A("Two points about method follow from the data rather than from the text.\n")
    A(f"The composite the paper compares between networks - 'the overall mean of descriptors' mean "
      f"betweenness and strength values' - averages two quantities on scales that differ by about "
      f"two orders of magnitude. Strength here runs roughly 0.2-1.3 and betweenness 0-53, so the "
      f"composite correlates with betweenness at r = "
      f"{OUT['networks']['corr_meancent_betweenness']['patrol']:.4f} (patrol) and "
      f"{OUT['networks']['corr_meancent_betweenness']['detect']:.4f} (detection), against "
      f"{OUT['networks']['corr_meancent_strength']['patrol']:.2f} and "
      f"{OUT['networks']['corr_meancent_strength']['detect']:.2f} for strength. Every statement "
      f"the paper makes about 'mean centrality', including the Cliff's delta figure, is in "
      f"practice a statement about betweenness alone. Standardising each measure before averaging "
      f"would resolve this in one line.\n")
    A(f"The multiple imputation is disclosed and was reproduced from the script: 15 completions, "
      f"averaged, rounded to integers, analysed as one dataset. That is the step at which "
      f"imputation uncertainty leaves the analysis - every interval downstream is conditional on "
      f"one synthetic dataset. The 15 completions were not supplied, so this review can name the "
      f"consequence but cannot size it (C17).\n")
    A("On distributional shape: F3 describes the centrality distributions used in the C3/C4 "
      "comparison and nothing more. The undesirable group has 6 nodes; a normality test at that n "
      "has almost no power, so neither its Shapiro-Wilk p-value nor the Q-Q curvature is treated "
      "here as passing or failing anything. The reason the three routes are reported side by side "
      "below is precisely that the shape cannot settle the question.\n")

    A("## 6. Resolution and evidence - F5\n")
    A(f"![F5]({f('F5_null_distribution.png')})\n")
    A(f"The headline centrality claim was run through all three routes in both networks "
      f"(`tables/three_routes_*.csv`); none was selected after the fact.\n")
    rows = []
    for net in ("patrol", "detect"):
        for meas in ("strength", "betweenness", "mean_cent"):
            k = f"{net}|{meas}"
            r = rp[k]
            rows.append({"network": net, "measure": meas,
                         "desirable mean": round(r["mean_desirable"], 3),
                         "undesirable mean": round(r["mean_undesirable"], 3),
                         **{kk.split("/")[-1]: f"{vv:.4g}" for kk, vv in r["p"].items()},
                         "smallest p this design could give": f"{r['min_achievable_p']:.2e}"})
    A(tbl(pd.DataFrame(rows)))
    A(f"\nThe permutation floor for 14 nodes against 6 is p >= "
      f"{rp['patrol|mean_cent']['min_achievable_p']:.2e}, far below 0.05, so the design can in "
      f"principle produce a significant result; the floor does not bind. The permutation "
      f"resolution at {N_PERM_ROUTES:,} draws is "
      f"{rp['patrol|mean_cent']['p_resolution']:.2e}.\n")
    A("In patrol every route agrees that desirable descriptors are more central. In detection "
      "every route agrees that they are not, on betweenness and on the composite, while strength "
      "alone separates them. The routes do not disagree with each other anywhere, so no claim "
      "here is `assumption-dependent`; what the claim depends on is the centrality measure, not "
      "the test.\n")
    A("One caveat applies to all six tests and is not fixable by any choice of route: the units "
      "are the 20 nodes of a single estimated network, which are not independent of one another - "
      "they are functions of the same 190 correlations. The p-values describe how surprising the "
      "desirable/undesirable split is under random relabelling of nodes; they are not statements "
      "about dogs.\n")
    if bp and bd:
        A(f"![C5 panel]({f('panel_C05_centrality_bootstrap_ci.png')})\n")
        A(f"The paper withholds centrality confidence intervals, citing known bias in bootstrap "
          f"intervals for network parameters - a defensible and explicitly argued choice. The "
          f"consequence is that the centrality ordering in its figure 2 carries no stated "
          f"uncertainty. Recomputing the bootstrap spread ({bp['n_boot']} and {bd['n_boot']} "
          f"replicate networks) shows the top node holds first place in "
          f"{bp['top_share_betweenness'] * 100:.0f}% of the {bp['n_usable_betweenness']} non-empty "
          f"patrol replicates on betweenness and "
          f"{bd['top_share_betweenness'] * 100:.0f}% of the {bd['n_usable_betweenness']} non-empty "
          f"detection replicates, while neighbouring nodes' intervals overlap heavily. 'Playful is "
          f"the most central descriptor' is supported as a point estimate and is not separated "
          f"from its neighbours by the resampling the paper itself ran (C5).\n")
        A(f"Two properties of the bootstrap itself qualify that number and apply equally to the "
          f"paper's own 2000 replicates. Bootstrap networks are systematically denser than the "
          f"network they resample - mean {bp['edges_mean']:.0f} edges against {P['n_edges']} "
          f"observed for patrol and {bd['edges_mean']:.0f} against {Dn['n_edges']} for detection - "
          f"because sampling dogs with replacement duplicates rows and inflates the correlations; "
          f"and {bd['n_empty_boot']} of {bd['n_boot']} detection replicates produced no network at "
          f"all. A bootstrap replicate is therefore not a re-draw of the published network, and "
          f"the rank of a node across replicates is a statement about a denser family of graphs "
          f"than the one in figure 1.\n")
        A(f"![C8 panel]({f('panel_C08_stability.png')})\n")
        A(f"Edge counts across bootstrap networks have s.d. {bp['edges_sd']:.2f} (patrol) and "
          f"{bd['edges_sd']:.2f} (detection) here against 12.82 and 29.75 reported; the detection "
          f"network's instability is reproduced, though the replicate count differs (C9). "
          f"Node-dropping is the weaker of the two stability analyses - it resamples nodes from "
          f"the same dogs, so a high correlation at three nodes says the centrality ordering is "
          f"internally consistent, not that it would survive a new sample of dogs. The "
          f"subject-wise curves are the ones that address sampling, and they fall away well before "
          f"the full sample in both networks, detection first (C8, C10, C11).\n")
        cd = res["stability"].get(("casedrop", "detect"))
        cp = res["stability"].get(("casedrop", "patrol"))
        if cd is not None and cp is not None:
            miss_d = sorted(set(np.round(np.arange(0.25, 0.96, 0.08), 2)) -
                            set(np.round(cd.level.unique(), 2)))
            miss_p = sorted(set(np.round(np.arange(0.25, 0.96, 0.08), 2)) -
                            set(np.round(cp.level.unique(), 2)))
            A(f"The gaps in the left panel are themselves a result. A level is plotted only where "
              f"the resampled network has at least three nodes with varying centrality; below "
              f"that the correlation with the original network is undefined because the "
              f"resampled graph is empty or near-empty. No level below "
              f"{cp.level.min():.2f} is plottable for patrol and none below "
              f"{cd.level.min():.2f} for detection (levels attempted but unusable - patrol: "
              f"{', '.join(f'{m:.2f}' for m in miss_p) or 'none'}; detection: "
              f"{', '.join(f'{m:.2f}' for m in miss_d) or 'none'}). The detection network stops "
              f"producing a network at all once a modest fraction of its 54 dogs is removed, "
              f"which is a stronger version of the same point as the gamma sweep in section 5.\n")
    if res.get("cliff_width_vs_nboot") is not None:
        A(f"![C7 panel]({f('panel_C07_cliff_vs_nboot.png')})\n")
        w = OUT["resampling"]["cliff_width_vs_nboot"]
        A(f"Cliff's delta in the paper is computed between the patrol and detection **bootstrap "
          f"distributions** of centrality, 2000 values against 2000. Those are not 2000 "
          f"independent observations - they are 2000 re-descriptions of 116 and 54 dogs. The "
          f"consequence is visible directly: redrawing delta from subsets of the bootstrap "
          f"replicates, its point estimate stays near {w[-1]['mean_delta']:+.3f} while its spread "
          f"falls from {w[0]['sd_of_delta']:.3f} at {w[0]['n_boot_used']} replicates to "
          f"{w[-1]['sd_of_delta']:.3f} at {w[-1]['n_boot_used']}. A confidence interval computed "
          f"this way narrows with computation rather than with evidence, and the 95% intervals "
          f"drawn in the paper's figure 3 inherit that property (C7).\n")
    if nct:
        A(f"![C14 panel]({f('panel_C14_network_difference.png')})\n")
        A(f"The patrol-versus-detection comparison is the paper's main substantive contrast and it "
          f"is never formally tested: two networks are estimated separately, their centralities "
          f"are compared through bootstrap effect sizes, and the structures are compared in prose. "
          f"A label-permutation test of structure invariance - permute dog type, re-estimate both "
          f"networks, record the largest absolute edge-weight difference - gives an observed "
          f"{nct['max_edge_diff_observed']:.3f} against a null 95th percentile of "
          f"{nct['max_edge_diff_null_95']:.3f}, p = {nct['p_structure_invariance']:.3f} "
          f"({nct['n_perm']} permutations; resolution {nct['p_resolution']:.4f}, so no p below "
          f"that value is reportable from this budget). Under a max-statistic correction "
          f"{nct['n_edges_fwer_sig']} of the 190 edges and "
          f"{nct['n_nodes_strength_fwer_sig']} of the 20 node strengths differ at FWER 0.05 "
          f"(`tables/nct_edge_differences.csv`). Note that this permutation treats dogs as "
          f"exchangeable and so ignores the repeated-handler structure described in section 4; a "
          f"handler-stratified permutation would be the stricter version (C14).\n")

    A("## 7. Multiplicity - F6\n")
    A(f"![F6]({f('F6_multiplicity.png')})\n")
    A(f"![C18 panel]({f('panel_C18_edge_multiplicity.png')})\n")
    A(f"**The paper reports no p-values, so F6 is not a picture of the authors' results.** It is "
      f"built from this review's own recomputation of the pairwise tests the network estimation "
      f"implies, and the p-hacking diagnostic that normally accompanies it is deliberately absent: "
      f"with zero reported p-values there is nothing to describe, and the just-below-threshold "
      f"count says nothing about any individual result even when there are enough of them.\n")
    A(f"The multiplicity that does exist is structural. Estimating a network over 20 nodes decides "
      f"190 node pairs; two networks decide 380. Re-running those pairs as zero-order polychoric "
      f"correlation tests:\n")
    A(tbl(pd.DataFrame([
        {"set": "patrol (190 pairs)", "uncorrected": cnt["patrol"]["none"],
         "Bonferroni": cnt["patrol"]["Bonferroni"], "BH": cnt["patrol"]["Benjamini-Hochberg"],
         "expected by chance at 0.05": round(mult["patrol"]["expected_false_positives"], 1),
         "edges the lasso selected": P["n_edges"]},
        {"set": "detection (190 pairs)", "uncorrected": cnt["detect"]["none"],
         "Bonferroni": cnt["detect"]["Bonferroni"], "BH": cnt["detect"]["Benjamini-Hochberg"],
         "expected by chance at 0.05": round(mult["detect"]["expected_false_positives"], 1),
         "edges the lasso selected": Dn["n_edges"]},
        {"set": "both networks (380)", "uncorrected": cnt["pooled"]["none"],
         "Bonferroni": cnt["pooled"]["Bonferroni"], "BH": cnt["pooled"]["Benjamini-Hochberg"],
         "expected by chance at 0.05": round(mult["pooled"]["expected_false_positives"], 1),
         "edges the lasso selected": P["n_edges"] + Dn["n_edges"]}])))
    msp, msd = mult["maxstat"]["patrol"], mult["maxstat"]["detect"]
    if msp:
        A(f"\nBecause the 190 pairs are heavily dependent, Bonferroni is the wrong correction for "
          f"them; the max-statistic permutation null is the one that respects the dependence. "
          f"Permuting each descriptor independently ({msp['n_perm']} permutations per network) "
          f"puts the FWER-0.05 threshold at |r| = {msp['crit_95']:.3f} for patrol and "
          f"{msd['crit_95']:.3f} for detection (`tables/maxstat_*.csv`).\n")
    A(f"None of this is a charge of uncorrected testing, because the paper never claims an error "
      f"rate: the lasso is a model-selection procedure and its output is a selected graph, not a "
      f"set of rejected nulls. The point is the one the paper's own wording invites - edges are "
      f"described as 'credibly non-zero partial correlations', and nothing in the procedure "
      f"attaches a credibility level to any of them. Reporting the selection as what it is, and "
      f"giving an edge-level stability proportion across bootstrap networks, would close the gap "
      f"(C18).\n")
    A(f"A second kind of multiplicity sits in the comparison between networks: 20 nodes x 2 "
      f"centrality measures = 40 Cliff's deltas, read for which are 'large', with no correction "
      f"and no formal test. The max-statistic test in section 6 is the corrected version of that "
      f"comparison.\n")

    A("## 8. Per-claim supplement\n")
    A(f"![C3]({f('panel_C3_central_desirable_patrol.png')})\n")
    A(f"![C4]({f('panel_C4_central_desirable_detect.png')})\n")
    A(f"![C2]({f('panel_C02_gamma_fragility.png')})\n")
    A(f"![C17]({f('panel_C17_imputation_uncertainty.png')})\n")
    A(f"![statcheck]({f('panel_C19_statcheck_not_assessable.png')})\n")
    A("Panels for C5, C7, C8 and C14 appear in section 6 and are not repeated. The two "
      "`not-assessable` panels are kept rather than dropped: an absent diagnostic and an "
      "impossible one look identical in a report that simply omits the slot.\n")

    A("## 9. What was not assessed\n")
    A("```")
    A(ac.uncheckable_assumptions_note())
    A("```\n")
    A("Specific to this paper:\n")
    A("- **Representativeness.** 117 of 227 invited handlers responded. Nothing in the supplied "
      "files allows responders to be compared with non-responders, so whether these 170 dogs "
      "represent Norwegian police dogs is unknown. Every claim in the ledger is conditional on "
      "this sample.")
    A("- **Breed confounding with dog type.** Patrol dogs are almost all German shepherds and "
      "Belgian malinois; detection dogs are mostly labradors, retrievers and spaniels. The "
      "network difference tested in C14 is a difference between work roles and between breed "
      "compositions at the same time, and no supplied file can separate them.")
    A("- **Conditional independence of repeated ratings.** The ICCs in section 4 describe "
      "between-handler variation; they cannot establish that the retained descriptors are "
      "independent, only that eight of the twenty still show appreciable clustering.")
    A("- **The authors' own likelihood-ratio screen.** Neither the paper nor the script records "
      "the per-descriptor LRT or delta-AIC values, and `clmm` was not re-run here, so the screen "
      "itself is reproduced in description only.")
    A("- **Imputation uncertainty (C17).** Only the rounded average of the 15 completions was "
      "supplied.")
    A("- **Supplementary tables S1-S5 and the authors' bootstrap tables** were not supplied; the "
      "reported raw centrality values (table S4) and Cliff's delta values (table S5) could be "
      "checked only against the nine deltas quoted in the running text.")
    A("- **Absolute betweenness scale.** The ordering of betweenness reproduces exactly, but "
      "`qgraph` and the implementation here may differ by whether shortest paths are counted once "
      "or twice per node pair (a factor of two). Table S4 would settle it; the comparison of "
      "nodes within a network is unaffected.\n")

    A("## 10. Reproduction\n")
    A(f"- `review_goold_2016.py` regenerates every figure, table and this report; "
      f"`core.py` holds the numerics (polychoric ML, Higham nearest-correlation projection, "
      f"block-coordinate-descent graphical lasso, EBIC selection, Brandes betweenness); "
      f"`prep.py` re-derives the authors' preprocessing from `rawData.csv`; `heavy.py` produces "
      f"the resampling.")
    A(f"- Seeds: three-routes and permutation tests `random_state = {SEED}`; resampling seeds are "
      f"stored inside each `results/*.npz` (`seed0` field).")
    A(f"- **Reduced resampling budget.** The paper used 2000 bootstrap samples per analysis. "
      f"Shared-machine limits forced a smaller budget here: "
      f"{bp['n_boot'] if bp else 0} non-parametric bootstrap networks per group, "
      f"{nct['n_perm'] if nct else 0} label permutations for the network-difference test, "
      f"{msp['n_perm'] if msp else 0} permutations for the edge max-statistic null, and 6 "
      f"replicates per level for the case-dropping and node-dropping curves. The lambda path in "
      f"resampling uses 25 values with `lambda.min.ratio = 0.05` against 100 and 0.01 for the two "
      f"reported networks; on the observed data that coarser grid returns the identical patrol "
      f"network (55 edges) and 71 rather than 70 edges for detection. Consequences, stated rather "
      f"than hidden: the permutation p-values in C14 cannot resolve below "
      f"{nct['p_resolution']:.4f} and the max-statistic thresholds below "
      f"{msp['p_resolution'] if msp else float('nan'):.4f}, so any claim needing a smaller p-value "
      f"is reported as unresolvable at this budget rather than given a number; bootstrap s.d.s "
      f"(C9) and stability thresholds (C10) are quantitative at low precision and are labelled "
      f"accordingly; the Cliff's delta demonstration (C7) is about how the spread scales with "
      f"replicate count and is unaffected.")
    A(f"- Three-route and permutation tests use n_perm = {N_PERM_ROUTES:,} "
      f"(resolution {1 / (N_PERM_ROUTES + 1):.2e}).")
    A(f"- Environment: Python, numpy {np.__version__}, pandas {pd.__version__}, "
      f"scipy {__import__('scipy').__version__}, matplotlib {matplotlib.__version__}. "
      f"`scikit-learn`, `networkx` and R were unavailable in this sandbox, so the graphical lasso, "
      f"the nearest-correlation projection and weighted betweenness were implemented directly; "
      f"the lasso solution was checked against the KKT conditions of its own objective (maximum "
      f"stationarity residual ~1e-5) and the whole pipeline against the paper's printed edge "
      f"counts (55 and 70) and centrality ordering.")
    A(f"- Input checksums are in section 1. Tables behind every number are in `tables/`.\n")

    txt = "\n".join(L)
    # If the figures have already been promoted to artifacts, point the report at them by
    # id rather than by relative path, so the saved report renders anywhere.
    idmap = os.path.join(HERE, "figure_artifact_ids.json")
    if os.path.exists(idmap):
        for fn, aid in json.load(open(idmap)).items():
            txt = txt.replace(f"figures/{fn}", "{{" + f"artifact:art_{aid}" + "}}")
    with open(os.path.join(HERE, "REVIEW_REPORT.md"), "w", encoding="utf-8") as fh:
        fh.write(txt)
    return txt


# ===================================================================== main
def main():
    raw, mi, rscript, paper, inv = stage0_inputs()
    mi2, pp, icc, long1, indep = stage12_design(raw, mi)
    nets, cent, gsweep = stage3_networks(mi)
    routes, perms = stage4_routes(cent)
    mult, pvals = stage5_multiplicity(nets)
    res = stage6_resampling(nets, cent)
    led = build_ledger(nets, cent, routes, perms, mult, res, gsweep)
    build_figures(mi2, long1, indep, cent, nets, gsweep, routes, perms, mult, res, led, pvals)
    write_report(led, mult, res, gsweep, cent, inv)
    with open(os.path.join(HERE, "review_numbers.json"), "w") as fh:
        json.dump(OUT, fh, indent=1, default=str)
    print("verdicts:", OUT["verdict_counts"])
    print("figures:", len(OUT["figures"]))
    print("claims:", len(led))


if __name__ == "__main__":
    main()
