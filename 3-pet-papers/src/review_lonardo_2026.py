"""
Statistical review of:

  Lonardo L., Volter C.J. & Huber L. (2026) Pupil size changes reveal that dogs
  are sensitive to the social conventions behind human gestures.
  Cognition 274:106586.  doi:10.1016/j.cognition.2026.106586

Regenerates every artifact under lonardo_2026/:
  figures/F1..F6 + per-claim panels + F7 pupil time course
  claim_ledger.csv, statcheck.csv, REVIEW_REPORT.md

The paper's full text was NOT available (Elsevier paywall); the claims reviewed
are the abstract's claims. The inputs are the authors' deposited data and their
own analysis code.

Run from the workspace root:  python lonardo_2026/review_lonardo_2026.py
"""

import hashlib
import os
import sys
import json
import datetime

import numpy as np
import pandas as pd
import scipy
import scipy.stats as st
import pyarrow
import pyarrow.parquet as pq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SKILL = os.path.join("skills", "statistical-review")
sys.path.insert(0, SKILL)

import kernel as K                                     # noqa: E402
from utilities import assumption_checks as ac          # noqa: E402
from utilities import claim_extraction as ce           # noqa: E402
from utilities import multiplicity as mult             # noqa: E402
from utilities import recompute as rc                  # noqa: E402
from utilities import resampling as rs                 # noqa: E402
from utilities import review_figures as rf             # noqa: E402

# ----------------------------------------------------------------------------
# 0. Fixed inputs, seeds, paths
# ----------------------------------------------------------------------------
SEED = 20260101
N_PERM = 10000
N_BOOT = 10000
ALPHA = 0.05

DATA = os.path.join("data_in", "dog_non_verbal_communication_ET")
OUT = "lonardo_2026"
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

P_SAMPLE = os.path.join(DATA, "data", "dog_nonverbal_communication_samplereport_ds.parquet")
P_ALL = os.path.join(DATA, "data", "dog_nonverbal_communication_samplereport_ds_ALL_DATA.parquet")
P_IA = os.path.join(DATA, "data", "IA_report_IP_from_reaction.txt")
P_IA579 = os.path.join(DATA, "data", "IA_report_IP_from_reaction_frame579.txt")
P_IANEW = os.path.join(DATA, "data", "IA_report_new_video_IP.txt")
P_CB = os.path.join(DATA, "data", "Non_verbal_communication_counterbalancing.csv")
P_RMD_PUP = os.path.join(DATA, "non_verbal_communication_pupil_size_analysis.Rmd")
P_RMD_DWL = os.path.join(DATA, "non_verbal_communication_dwell_time_analysis.rmd")


def sha256(path, n=16):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()[:n]


# ----------------------------------------------------------------------------
# 1. The abstract (the only paper text available) -> claim inventory
# ----------------------------------------------------------------------------
ABSTRACT = """While many animals display behaviours resembling social conventions when \
interacting with conspecifics, it is unknown whether animals who live in close contact with \
humans develop sensitivity to our social conventions. With an eye-tracking study (N = 39), we \
asked whether dogs' (Canis familiaris) eye movements and pupil sizes would show a reaction when \
observing a violation in the customary response to human conventional communicative gestures \
('come here' and 'go away'). In the video stimuli, a human actress either approached or walked \
away from a second human actress. Both these actions could be expected or surprising \
(within-subject), depending on whether they happened in response to the 'come here' or the 'go \
away' gesture of the second actress. While dogs' looking times and number of gaze transitions \
between the two actresses did not differ significantly between congruent and incongruent \
responses, their pupil size changes did: dogs' pupils dilated more when watching the incongruent \
responses compared to the congruent ones. We discuss the possibility that pet dogs, during their \
ontogeny, implicitly form expectations about the usual response to some conventional human \
gestures. These expectations are not self-centred, as they were evidenced when dogs observed a \
third-party interaction between two humans, and they are context-specific, as the same responses \
elicited surprise following one gesture but not the other."""

with open(os.path.join(OUT, "paper_abstract.txt"), "w", encoding="utf-8") as fh:
    fh.write(ABSTRACT + "\n")

# SECOND PASS: the published PDF was supplied after the first version of this
# review. Everything below now reads the paper's own text; the first pass had
# only the abstract above.
import pypdfium2 as pdfium                                   # noqa: E402

PDF = ("C:/Users/sungk/.claude-science/orgs/5885b251-f542-45c9-b932-d3bdb36b54c4/"
       "artifacts/proj_cca0f303f25d/243eff84-837a-4b59-9483-8329ea64eeb9/"
       "vec3193b9_1-s2.0-S0010027726001538-main.pdf")
_doc = pdfium.PdfDocument(PDF)
PAGES = [p.get_textpage().get_text_range() for p in _doc]
with open(os.path.join(OUT, "paper_fulltext.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n\f\n".join(PAGES))

inventory = ce.build_claim_inventory(PAGES)
quality = ce.pdf_text_quality(PAGES)
inv_summary = ce.inventory_summary(inventory)
N_AUTO_N = int(len(inventory.get("sample_sizes", [])))
N_AUTO_P = int(len(inventory.get("p_values", [])))
N_AUTO_T = int(len(inventory.get("tests", [])))

# The automatic extractor reads t and F but not the Greek chi-square glyph, so
# the six chi-square triples are transcribed by hand from the sentences the
# p-value extractor did capture. Source sentence recorded with each.
HAND_TESTS = pd.DataFrame([
    dict(claim_id="HND001", page=4, stat_type="chi2", stat_value=0.05, df1=1, df2=np.nan,
         p_relation="=", p_reported=0.818, tail_hint=2,
         sentence="dwell proportion GLMM: trial type x agent type interaction"),
    dict(claim_id="HND002", page=4, stat_type="chi2", stat_value=0.12, df1=1, df2=np.nan,
         p_relation="=", p_reported=0.731, tail_hint=2,
         sentence="dwell proportion GLMM: trial type (reduced model)"),
    dict(claim_id="HND003", page=4, stat_type="chi2", stat_value=13.85, df1=1, df2=np.nan,
         p_relation="<", p_reported=0.001, tail_hint=2,
         sentence="dwell proportion GLMM: agent type"),
    dict(claim_id="HND004", page=5, stat_type="chi2", stat_value=10.32, df1=6, df2=np.nan,
         p_relation="=", p_reported=0.112, tail_hint=2,
         sentence="dwell proportion GLMM: condition + condition x agent type (full vs reduced)"),
    # The paper prints "chi2 < 0.01, df = 1, p = .951": the '<' bounds the STATISTIC,
    # while p is reported exactly. p_relation therefore tracks the operator before p
    # ('='), and the statistic's bound is carried in STAT_BOUND below. Labelling this
    # row '<' would send it to statcheck's one-sided branch, which only asks whether
    # the recomputed p falls below the reported one -- a near-automatic pass that
    # would hide a real gap instead of testing for one.
    dict(claim_id="HND005", page=5, stat_type="chi2", stat_value=0.01, df1=1, df2=np.nan,
         p_relation="=", p_reported=0.951, tail_hint=2,
         sentence="transition count GLMM: trial type (statistic reported as chi2 < 0.01)"),
    dict(claim_id="HND006", page=5, stat_type="chi2", stat_value=2.95, df1=1, df2=np.nan,
         p_relation="=", p_reported=0.086, tail_hint=2,
         sentence="transition count GLMM: session number"),
])
ALL_TESTS = pd.concat([inventory["tests"], HAND_TESTS], ignore_index=True)
inventory_checked = dict(inventory)
inventory_checked["tests"] = ALL_TESTS
audit = K.audit_reported_statistics(inventory_checked, alpha=ALPHA)
STATCHECK = audit["statcheck"]
STATCHECK.to_csv(os.path.join(OUT, "statcheck.csv"), index=False)

# statcheck recomputes p from the statistic AS PRINTED. Where a statistic is
# printed to two decimals its own rounding moves the implied p; propagating that
# window decides whether a flagged row is an arithmetic error or a rounding
# artefact. Computed, not asserted.
# Rows where the printed statistic is an upper BOUND rather than a rounded point:
# the implied window is [0, bound], not +/- half a printed digit.
STAT_BOUND = {"HND005": (0.0, 0.01)}


def _rounding_resolved(row, decimals=2):
    """Window of p implied by the precision with which the STATISTIC was printed."""
    if row["stat_type"] != "chi2" or not np.isfinite(row.get("df1", np.nan)):
        return None
    if row["claim_id"] in STAT_BOUND:
        s_lo, s_hi = STAT_BOUND[row["claim_id"]]
        kind = "bound"
    else:
        half = 0.5 * 10 ** (-decimals)
        s_lo, s_hi = row["stat_value"] - half, row["stat_value"] + half
        kind = "rounding"
    lo = float(st.chi2.sf(s_hi, row["df1"]))
    hi = float(st.chi2.sf(max(s_lo, 0.0), row["df1"]))
    return bool(lo <= row["p_reported"] <= hi), lo, hi, kind


ROUNDING = {}
for _, row in STATCHECK.iterrows():
    if row.get("consistent") is False:
        res = _rounding_resolved(row)
        if res is not None:
            ROUNDING[row["claim_id"]] = dict(stat=row["stat_value"], df1=row["df1"],
                                             p_reported=row["p_reported"],
                                             p_point=row["p_recomputed"],
                                             resolved=res[0], lo=res[1], hi=res[2],
                                             kind=res[3], sentence=row.get("sentence", ""))
N_SC_CHECKED = int(STATCHECK["consistent"].notna().sum()) if len(STATCHECK) else 0
N_SC_FLAGGED = int((STATCHECK["consistent"] == False).sum()) if len(STATCHECK) else 0  # noqa: E712
N_SC_ROUNDING = int(sum(1 for v in ROUNDING.values() if v["resolved"]))
N_SC_NOTCHECK = int(len(STATCHECK) - N_SC_CHECKED)

# The paper's own numbers, transcribed from the sentences quoted in section 1.
REPORTED = dict(
    pupil_gamm_t=2.04, pupil_gamm_p=0.041,
    pupil_wilcox_Tplus=313.0, pupil_wilcox_n=28, pupil_wilcox_p=0.011,
    pupil_session_t=-0.01, pupil_session_p=0.990,
    pupil_n39_t=0.54, pupil_n39_p=0.592,
    curve_lo=7697.0, curve_hi=9111.0,
    dwell_trialtype_chi2=0.12, dwell_trialtype_p=0.731,
    dwell_interaction_chi2=0.05, dwell_interaction_p=0.818,
    dwell_agenttype_chi2=13.85, dwell_agenttype_p=0.001,
    dwell_condition_chi2=10.32, dwell_condition_df=6, dwell_condition_p=0.112,
    trans_trialtype_p=0.951, trans_session_chi2=2.95, trans_session_p=0.086,
    reacting_mean=5105.0, reacting_sd=1924.0, communicating_mean=3218.0, communicating_sd=2065.0,
    n_total=39, n_pupil=28, n_female=18, mean_age_months=59.3,
    n_excluded_trials=12, n_excluded_pupil_lost=4, n_excluded_offscreen=8,
)

# ----------------------------------------------------------------------------
# 2. Load the data
# ----------------------------------------------------------------------------
pup_ds = pq.read_table(P_SAMPLE).to_pandas()        # authors' "dat"  (downsampled, QC'd)
pup_all = pq.read_table(P_ALL).to_pandas()          # authors' "dat2" (ALL_DATA, N = 39)
cb = pd.read_csv(P_CB)
ia = pd.read_csv(P_IA, sep="\t", na_values=".", low_memory=False)
ia579 = pd.read_csv(P_IA579, sep="\t", na_values=".", low_memory=False)
ianew = pd.read_csv(P_IANEW, sep="\t", na_values=".", low_memory=False)

# The 11 dogs the authors exclude from the primary pupil analysis
EXCL = ["Ace", "Filou", "Floki", "George", "Ignaz", "Jimi", "Lenny",
        "Sahibu", "Sheila", "Timo", "Zserbo"]

n_cond = pup_ds.groupby("subject", observed=True)["condition"].nunique()
incomplete = sorted(n_cond[n_cond < 4].index.astype(str))
assert incomplete == sorted(EXCL), "authors' exclusion list no longer matches the data"

pup_28 = pup_ds[~pup_ds["subject"].astype(str).isin(EXCL)].copy()

# ----------------------------------------------------------------------------
# 3. Subject bookkeeping (stage 1)
# ----------------------------------------------------------------------------
cb["EDF3f"] = cb["EDF3"].replace({"Georgia3": "Georgia_3", "Mathild3": "Mathilda_3"})
cb["subject_from_EDF3"] = cb["EDF3f"].astype(str).str.split("_").str[0]
named = cb[cb["Name"].notna()]
subj_parquet = sorted(pup_ds["subject"].astype(str).unique())
join_fail = sorted(set(subj_parquet) - set(named["subject_from_EDF3"]))
# dogs sharing a Dog_ID share the counterbalancing `order` string -> Dog_ID is an
# order slot, not a dog identifier. Verified, not assumed:
order_per_dogid = named.groupby("Dog_ID")["order"].nunique()
dogid_is_order = bool((order_per_dogid == 1).all())
dogid_shared = order_per_dogid[named.groupby("Dog_ID").size() > 1]

NBOOK = dict(
    n_rows_cb=int(len(cb)),
    n_named=int(len(named)),
    n_unique_dogid=int(named["Dog_ID"].nunique()),
    n_blank_rows=int(cb["Name"].isna().sum()),
    n_subject_parquet_ds=int(pup_ds["subject"].nunique()),
    n_subject_parquet_all=int(pup_all["subject"].nunique()),
    n_primary_pupil=int(pup_28["subject"].nunique()),
    dogid_is_order_slot=dogid_is_order,
    n_shared_dogid=int(len(dogid_shared)),
    demo_join_failures=join_fail,
    n_female=int((named["sex"] == "f").sum()),
    n_breeds=int(named["Breed"].nunique()),
)

# ----------------------------------------------------------------------------
# 4. Replicate the authors' dwell-time pipeline (dwell .rmd, verbatim logic)
# ----------------------------------------------------------------------------
RECODE = {"Georgia1": "Georgia_1", "Georgia2": "Georgia_2", "Georgia3": "Georgia_3",
          "Georgia4": "Georgia_4", "Mathild1": "Mathilda_1", "Mathild2": "Mathilda_2",
          "Mathild3": "Mathilda_3", "Mathild4": "Mathilda_4"}
DROP_SESSIONS = ["Nala_2", "Sahibu_2", "Melody_1", "Zserbo_1"]

rip = pd.concat([ia[~ia["RECORDING_SESSION_LABEL"].isin(DROP_SESSIONS)],
                 ia579[ia579["RECORDING_SESSION_LABEL"] == "Zserbo_1"]], ignore_index=True)
rip["RECORDING_SESSION_LABEL"] = rip["RECORDING_SESSION_LABEL"].replace(RECODE)
rip[["subject", "session_num"]] = rip["RECORDING_SESSION_LABEL"].str.split("_", n=1, expand=True)
rip = rip[rip["session_number"].notna()].copy()

dw = rip[rip["IA_LABEL"] != "FIXATION_TRIGGER"].copy()
dw["agent_type"] = np.where(dw["IA_LABEL"].str.contains("inviting|pointing", regex=True),
                            "communicating", "reacting")
dw["prop"] = dw["IA_DWELL_TIME"] / dw["TRIAL_DWELL_TIME"]

# reproduction check against the two values the authors wrote into their own script
_agent = dw.groupby(["subject", "agent_type"])["IA_DWELL_TIME"].mean().unstack().mean()
REPRO_AGENT = {"reacting": float(_agent["reacting"]),
               "communicating": float(_agent["communicating"])}
assert abs(REPRO_AGENT["reacting"] - 5105.282) < 0.01
assert abs(REPRO_AGENT["communicating"] - 3217.737) < 0.01


def count_transitions(seq):
    """The authors' count_transitions(), transcribed from the dwell-time .rmd."""
    if not isinstance(seq, str):
        return np.nan
    s = seq.replace("[", "").replace("]", "").replace(" ", "")
    vals = []
    for part in s.split(","):
        try:
            vals.append(float(part))
        except ValueError:
            pass
    tr, prev = 0, None
    for v in vals:
        if v in (21.0, 22.0, 23.0):
            if prev is not None and prev != 20.0 and prev != v:
                tr += 1
            prev = v
    return tr


rip["transition_count"] = rip["INTEREST_AREA_FIXATION_SEQUENCE"].map(count_transitions)
trd = (rip[rip["IA_LABEL"] != "FIXATION_TRIGGER"]
       .groupby(["subject", "trial_number", "condition", "trial_type", "TRIAL_DWELL_TIME"],
                as_index=False)
       .agg(transition_count=("transition_count", "first")))

# ----------------------------------------------------------------------------
# 5. Subject-level paired tables. Subject is the unit of analysis throughout.
# ----------------------------------------------------------------------------
def paired(df, value, by="trial_type"):
    return df.groupby(["subject", by], observed=True)[value].mean().unstack(by).dropna()


W = {
    "pupil, primary complete cases (dat3)": paired(pup_28, "pupil_base"),
    "pupil, ds file, no exclusion": paired(pup_ds, "pupil_base"),
    "pupil, ALL_DATA (abstract's N = 39)": paired(pup_all, "pupil_base"),
    "looking time to agents (ms)": dw.groupby(["subject", "trial_type"])["IA_DWELL_TIME"].mean().unstack().dropna(),
    "looking time, proportion of trial": dw.groupby(["subject", "trial_type"])["prop"].mean().unstack().dropna(),
    "gaze transitions per trial": trd.groupby(["subject", "trial_type"])["transition_count"].mean().unstack().dropna(),
}
KEY28, KEYDS, KEYALL = list(W)[0], list(W)[1], list(W)[2]
KEYDWL, KEYPROP, KEYTR = list(W)[3], list(W)[4], list(W)[5]


def routes_block(w, label, seed=SEED):
    """All available routes + interval + floors for one paired contrast.

    Effect is signed incongruent - congruent, so a positive effect is the
    direction the abstract claims for pupil size.
    """
    x = w["incong"].to_numpy(float)
    y = w["cong"].to_numpy(float)
    d = x - y
    r = rs.three_routes(x, y, paired=True, n_perm=N_PERM, random_state=seed)
    p = rs.permutation_test(x, y, paired=True, n_perm=N_PERM, random_state=seed)
    b = rs.bootstrap_ci(d, np.mean, n_boot=N_BOOT, random_state=seed)
    return dict(label=label, n=int(len(w)), effect=float(d.mean()), sd=float(d.std(ddof=1)),
                mean_cong=float(y.mean()), mean_incong=float(x.mean()),
                ci=(float(b["ci_low"]), float(b["ci_high"])), boot=b, perm=p, routes=r,
                sign=rs.sign_test(d), shape=ac.paired_difference_check(x, y),
                p_t=float(r.loc[0, "p_value"]), p_rank=float(r.loc[1, "p_value"]),
                p_sign=float(r.loc[2, "p_value"]), p_perm=float(p["p_value"]),
                min_p=float(p["min_achievable_p"]), sign_floor=float(rs.min_sign_test_p(len(w))),
                p_resolution=float(p["p_resolution"]),
                routes_agree=bool(r["agrees_at_0.05"].nunique() == 1),
                n_pos=int((d > 0).sum()))


R = {k: routes_block(w, k) for k, w in W.items()}

# --- direct reproduction of the paper's own reported statistics ---------------
# The paper reports the subject-level pupil test as T+ = 313, N = 28, p = .011,
# using the authors' better.wilc() convention T+ = max(W, N(N+1)/2 - W).
_d28 = (W[KEY28]["incong"] - W[KEY28]["cong"]).to_numpy(float)
_n28 = int(len(_d28))
_W = float(st.wilcoxon(_d28).statistic)
REPRO = dict(
    wilcox_W=_W,
    wilcox_Tplus=float(max(_W, _n28 * (_n28 + 1) / 2 - _W)),
    wilcox_n=_n28,
    wilcox_p=float(st.wilcoxon(_d28).pvalue),
    n_ties=int((_d28 == 0).sum()),
)
REPRO["Tplus_matches"] = bool(REPRO["wilcox_Tplus"] == REPORTED["pupil_wilcox_Tplus"])
REPRO["p_matches_3dp"] = bool(round(REPRO["wilcox_p"], 3) == REPORTED["pupil_wilcox_p"])

# Trial accounting against the methods section, as a containment chain.
_set_all = set(map(tuple, pup_all[["subject", "condition"]].astype(str).drop_duplicates().values))
_set_ds = set(map(tuple, pup_ds[["subject", "condition"]].astype(str).drop_duplicates().values))
_trialmean_all = pup_all.groupby(["subject", "condition"], observed=True)["pupil_base"].mean()
_empty = sorted(f"{s}/{c}" for s, c in _trialmean_all[_trialmean_all.isna()].index)
_dropped = sorted(f"{s}/{c}" for s, c in (_set_all - _set_ds))
TRIALS = dict(
    intended=39 * 4,
    alldata=len(_set_all),
    alldata_with_pupil=len(_set_all) - len(_empty),
    ds=len(_set_ds),
    ds_is_subset=bool(_set_ds <= _set_all),
    empty_trials=_empty,
    dropped=len(_dropped),
    dropped_list=_dropped,
    text_says_dropped=REPORTED["n_excluded_trials"],
)
TRIALS["missing_sessions"] = TRIALS["intended"] - TRIALS["alldata"]
TRIALS["dropped_with_pupil"] = TRIALS["dropped"] - len(_empty)
TRIALS["gap_vs_text"] = TRIALS["text_says_dropped"] - TRIALS["dropped"]
TRIALS["gap_vs_text_strict"] = TRIALS["text_says_dropped"] - TRIALS["dropped_with_pupil"]

# The paper prints "5105 +/- 1924 ms" without naming the error. Now checkable.
_g_agent = dw.groupby(["subject", "agent_type"])["IA_DWELL_TIME"].mean().unstack()
ERRBAR = {a: rc.error_bar_identity(e, float(_g_agent[a].std(ddof=1)), int(len(_g_agent)))
          for a, e in (("reacting", REPORTED["reacting_sd"]),
                       ("communicating", REPORTED["communicating_sd"]))}

# Demographics against the Subjects paragraph.
_age = pd.to_numeric(named["age"], errors="coerce")
DEMO = dict(n_named=int(len(named)), n_female=int((named["sex"] == "f").sum()),
            mean_age=float(_age.mean()), n_age=int(_age.notna().sum()))
DEMO["female_matches"] = bool(DEMO["n_female"] == REPORTED["n_female"])
DEMO["age_matches_1dp"] = bool(round(DEMO["mean_age"], 1) == REPORTED["mean_age_months"])

# ----------------------------------------------------------------------------
# 6. Time-resolved cluster-based permutation (the estimand the authors' GAMM
#    targets: a difference localised in time, not averaged over the window)
# ----------------------------------------------------------------------------
def cluster_perm(df, drop=None, n_perm=N_PERM, seed=SEED):
    d = df if drop is None else df[~df["subject"].astype(str).isin(drop)]
    m = (d.groupby(["subject", "trial_type", "time"], observed=True)["pupil_base"]
           .mean().unstack("trial_type").dropna())
    diff = (m["incong"] - m["cong"]).unstack("time").dropna(axis=0, how="any")
    X = diff.to_numpy(float)
    times = diff.columns.to_numpy(float)
    n = X.shape[0]
    thr = float(st.t.ppf(1 - ALPHA / 2, n - 1))

    def tstat(M):
        return M.mean(0) / (M.std(0, ddof=1) / np.sqrt(n))

    def clusters(tv):
        out, cur = [], None
        for i, v in enumerate(tv):
            if abs(v) > thr:
                cur = [i, i, v] if cur is None else [cur[0], i, cur[2] + v]
            elif cur is not None:
                out.append(tuple(cur)); cur = None
        if cur is not None:
            out.append(tuple(cur))
        return out

    tobs = tstat(X)
    obs = clusters(tobs)
    mass = max([abs(c[2]) for c in obs], default=0.0)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        s = rng.choice([-1.0, 1.0], size=n)[:, None]
        null[i] = max([abs(c[2]) for c in clusters(tstat(X * s))], default=0.0)
    win = [(float(times[c[0]]), float(times[c[1]]), float(c[2])) for c in obs]
    big = max(win, key=lambda w: abs(w[2])) if win else (np.nan, np.nan, 0.0)
    return dict(n=n, times=times, t=tobs, thr=thr, X=X, null=null, obs_mass=float(mass),
                windows=win, biggest=big,
                p_cluster=float((1 + int((null >= mass).sum())) / (n_perm + 1)),
                p_resolution=1.0 / (n_perm + 1))


CP28 = cluster_perm(pup_ds, drop=EXCL)
CPALL = cluster_perm(pup_all)
CW0, CW1 = CP28["biggest"][0] + 5000.0, CP28["biggest"][1] + 5000.0

# ----------------------------------------------------------------------------
# 7. Context-specificity: same response, different preceding gesture
#    approach response: GA (incongruent) vs CA (congruent)
#    leaving response:  CG (incongruent) vs GG (congruent)
# ----------------------------------------------------------------------------
def cond_pair(df, c_incong, c_cong, drop=EXCL):
    d = df[~df["subject"].astype(str).isin(drop)]
    d = d[d["condition"].astype(str).isin([c_incong, c_cong])]
    w = d.groupby(["subject", "condition"], observed=True)["pupil_base"].mean().unstack("condition").dropna()
    return w.rename(columns={c_incong: "incong", c_cong: "cong"})


KEY_APP = "pupil, approach response (GA vs CA)"
KEY_LEA = "pupil, leaving response (CG vs GG)"
CTX = {
    KEY_APP: routes_block(cond_pair(pup_ds, "GA", "CA"), KEY_APP),
    KEY_LEA: routes_block(cond_pair(pup_ds, "CG", "GG"), KEY_LEA),
}
_wA = cond_pair(pup_ds, "GA", "CA")
_wL = cond_pair(pup_ds, "CG", "GG")
_common = _wA.index.intersection(_wL.index)
_dA = (_wA.loc[_common, "incong"] - _wA.loc[_common, "cong"]).to_numpy(float)
_dL = (_wL.loc[_common, "incong"] - _wL.loc[_common, "cong"]).to_numpy(float)
INTER = dict(
    n=int(len(_common)),
    effect=float((_dA - _dL).mean()),
    routes=rs.three_routes(_dA, _dL, paired=True, n_perm=N_PERM, random_state=SEED),
    perm=rs.permutation_test(_dA, _dL, paired=True, n_perm=N_PERM, random_state=SEED),
    boot=rs.bootstrap_ci(_dA - _dL, np.mean, n_boot=N_BOOT, random_state=SEED),
)

# ----------------------------------------------------------------------------
# 8. Claims. review_claim() for every contrast the data can address; hand rows
#    for the claims they cannot.
# ----------------------------------------------------------------------------
def long_frame(w, value_name, incong_label="incongruent", cong_label="congruent"):
    """Long format with the incongruent level first, so review_claim's
    levels[0] - levels[1] is incongruent - congruent."""
    rows = ([{"subject": s, "trial_type": incong_label, value_name: float(v)}
             for s, v in w["incong"].items()]
            + [{"subject": s, "trial_type": cong_label, value_name: float(v)}
               for s, v in w["cong"].items()])
    return pd.DataFrame(rows)


PUP_LAB = "pupil size change (baseline-corrected, a.u.)"
DWL_LAB = "mean dwell time to agents (ms)"
TR_LAB = "gaze transitions per trial"

c1 = K.review_claim(long_frame(W[KEY28], PUP_LAB), value_col=PUP_LAB,
                    group_col="trial_type", pair_col="subject", paired=True,
                    claim_id="C1", claim_short="pupils dilated more to incongruent responses",
                    n_perm=N_PERM, n_boot=N_BOOT, alpha=ALPHA, random_state=SEED,
                    figure_path=os.path.join(FIG, "panel_C1.png"))

PROP_LAB = "proportion of trial dwell time on the agents"
c2 = K.review_claim(long_frame(W[KEYPROP], PROP_LAB), value_col=PROP_LAB, group_col="trial_type",
                    pair_col="subject", paired=True, claim_id="C2",
                    claim_short="looking times did not differ cong vs incong",
                    n_perm=N_PERM, n_boot=N_BOOT, alpha=ALPHA, random_state=SEED,
                    figure_path=os.path.join(FIG, "panel_C2.png"))

c3 = K.review_claim(long_frame(W[KEYTR], TR_LAB), value_col=TR_LAB, group_col="trial_type",
                    pair_col="subject", paired=True, claim_id="C3",
                    claim_short="gaze transitions did not differ cong vs incong",
                    n_perm=N_PERM, n_boot=N_BOOT, alpha=ALPHA, random_state=SEED,
                    figure_path=os.path.join(FIG, "panel_C3.png"))

c5 = K.review_claim(long_frame(cond_pair(pup_ds, "CG", "GG"), PUP_LAB,
                               "incongruent (CG)", "congruent (GG)"),
                    value_col=PUP_LAB, group_col="trial_type", pair_col="subject", paired=True,
                    claim_id="C5", claim_short="effect is context-specific (one gesture only)",
                    n_perm=N_PERM, n_boot=N_BOOT, alpha=ALPHA, random_state=SEED,
                    figure_path=os.path.join(FIG, "panel_C5.png"))

# -- reviewer verdicts. propose_verdict() cannot see the methods or the second
# -- sample, so each override is recorded with the numbers that drove it.
proposed = {c["claim_id"]: c["verdict"] for c in (c1, c2, c3, c5)}
A, B, C = R[KEY28], R[KEYALL], R[KEYDS]
DWL, TR = R[KEYDWL], R[KEYTR]
LEA, APP = CTX[KEY_LEA], CTX[KEY_APP]

c1["verdict"] = "supported"
c1["reasons"] = [
    f"the paper's own subject-level test reproduces exactly: it reports T+ = "
    f"{REPORTED['pupil_wilcox_Tplus']:.0f}, N = {REPORTED['pupil_wilcox_n']}, p = "
    f"{REPORTED['pupil_wilcox_p']:.3f}; recomputing the Wilcoxon signed-rank on the same "
    f"{REPRO['wilcox_n']} dogs from the raw parquet gives W = {REPRO['wilcox_W']:.0f}, hence T+ = "
    f"max(W, N(N+1)/2 - W) = {REPRO['wilcox_Tplus']:.0f} (match: {REPRO['Tplus_matches']}) and "
    f"p = {REPRO['wilcox_p']:.4f} (match to 3 dp: {REPRO['p_matches_3dp']})",
    f"every route agrees at alpha on those {A['n']} dogs: paired t p = {A['p_t']:.4f}, Wilcoxon "
    f"p = {A['p_rank']:.4f}, sign test p = {A['p_sign']:.4f}, sign-flip permutation p = "
    f"{A['p_perm']:.4f}; incongruent minus congruent = {A['effect']:.1f} a.u. (95% bootstrap CI "
    f"{A['ci'][0]:.1f} to {A['ci'][1]:.1f})",
    f"the time-resolved cluster permutation run here, blind to the paper, places the difference at "
    f"{CW0:.0f}-{CW1:.0f} ms (cluster p = {CP28['p_cluster']:.4f}), inside the "
    f"{REPORTED['curve_lo']:.0f}-{REPORTED['curve_hi']:.0f} ms window the paper's GAMM difference "
    f"curve reports",
    f"the paper reports the N = 39 analysis itself and reports it as null (t = "
    f"{REPORTED['pupil_n39_t']}, p = {REPORTED['pupil_n39_p']:.3f}); the subject-level contrast "
    f"recomputed here on that sample agrees in conclusion though not in test ({B['effect']:.1f} "
    f"a.u., 95% CI {B['ci'][0]:.1f} to {B['ci'][1]:.1f}, permutation p = {B['p_perm']:.3f})",
    f"the first pass of this review, which had only the abstract, recorded the 28-vs-39 split as an "
    f"undisclosed sample-definition dependence and returned supported-with-caveat; the methods and "
    f"results sections state N for every analysis, so that caveat was the reviewer's and is withdrawn",
]
PROP = R[KEYPROP]
c2["verdict"] = "supported-with-caveat"
c2["reasons"] = [
    f"the non-difference reproduces on the paper's own DV: it reports chi2 = "
    f"{REPORTED['dwell_trialtype_chi2']}, df = 1, p = {REPORTED['dwell_trialtype_p']:.3f} for trial "
    f"type on the proportion of dwell time, and the subject-level contrast recomputed here is "
    f"{PROP['effect']:+.4f} (95% bootstrap CI {PROP['ci'][0]:.4f} to {PROP['ci'][1]:.4f}), all "
    f"routes p >= {min(PROP['p_t'], PROP['p_rank'], PROP['p_perm']):.2f}",
    f"in the interpretable unit the paper uses for its descriptives, the same contrast is "
    f"{DWL['effect']:.1f} ms on a congruent mean of {DWL['mean_cong']:.0f} ms (95% CI "
    f"{DWL['ci'][0]:.0f} to {DWL['ci'][1]:.0f} ms, permutation p = {DWL['p_perm']:.3f})",
    f"CAVEAT: the interval admits differences up to about "
    f"{max(abs(DWL['ci'][0]), abs(DWL['ci'][1])) / DWL['mean_cong'] * 100:.0f}% of the congruent "
    f"mean, so the data are consistent with a null but do not establish one; the paper states no "
    f"smallest effect of interest and reports no equivalence test, and its bootstrap CIs are in "
    f"Table S2, which was not supplied",
    f"automatic proposal: {proposed['C2']}",
]
c3["verdict"] = "supported-with-caveat"
c3["reasons"] = [
    f"the non-difference reproduces: incongruent minus congruent = {TR['effect']:.3f} transitions on a "
    f"congruent mean of {TR['mean_cong']:.2f}; all routes p >= "
    f"{min(TR['p_t'], TR['p_rank'], TR['p_perm']):.2f}",
    f"the 95% bootstrap interval runs {TR['ci'][0]:.2f} to {TR['ci'][1]:.2f} transitions, about "
    f"+/-{max(abs(TR['ci'][0]), abs(TR['ci'][1])) / TR['mean_cong'] * 100:.0f}% of the congruent mean; "
    f"consistent with a null, not an equivalence demonstration",
    f"{int((trd['transition_count'] == 0).sum())} of {len(trd)} trials contain zero transitions, so the "
    f"subject-level mean summarises a count with a hard floor at 0",
    f"automatic proposal: {proposed['C3']}",
]
c5["verdict"] = "supported-with-caveat"
c5["reasons"] = [
    "the paper argues context-relativity from the design, not from a per-gesture test: because the "
    "same physical action is congruent after one gesture and incongruent after the other, the "
    "action cancels in the trial-type contrast, and the paper's evidence is the main effect of "
    f"trial type (t = {REPORTED['pupil_gamm_t']}, p = {REPORTED['pupil_gamm_p']:.3f}; Wilcoxon "
    f"p = {REPORTED['pupil_wilcox_p']:.3f})",
    f"the design premise verifies exactly in the data: all {NBOOK['n_primary_pupil']} dogs in the "
    f"pupil sample contribute one trial in each of CA, CG, GA and GG, so approaching appears once as "
    f"congruent and once as incongruent within every dog, and likewise for leaving - the contrast is "
    f"action-balanced within subject, not merely on average",
    f"the gesture-by-congruency interaction, which the paper does not test, is null here: "
    f"{INTER['effect']:.1f} a.u. (95% CI {INTER['boot']['ci_low']:.1f} to "
    f"{INTER['boot']['ci_high']:.1f}), permutation p = {INTER['perm']['p_value']:.3f}; that is "
    f"consistent with an expectation that operates in both gesture contexts and is not evidence "
    f"against the claim",
    f"caveat: nothing in the design separates a context-specific expectation from a general response "
    f"to incongruity, and the separate per-gesture contrasts are not interpretable on their own - "
    f"leaving CG vs GG is {LEA['effect']:.1f} a.u. with routes disagreeing at alpha (t "
    f"{LEA['p_t']:.4f}, permutation {LEA['p_perm']:.4f}, Wilcoxon {LEA['p_rank']:.4f}, sign "
    f"{LEA['p_sign']:.4f}) while approach GA vs CA is {APP['effect']:.1f} a.u. (all routes p >= "
    f"{min(APP['p_t'], APP['p_perm']):.2f}); reading those two side by side would be a difference of "
    f"significance, and the paper does not do so",
    "the first pass of this review, which had only the abstract, read 'the same responses elicited "
    "surprise following one gesture but not the other' as a per-gesture significance split and "
    "returned assumption-dependent; the full text shows it restates the within-subject design, so "
    "that verdict was the reviewer's misreading and is withdrawn",
]
c3["reasons"].insert(2, (
    f"the paper reports chi2 < 0.01, df = 1, p = {REPORTED['trans_trialtype_p']:.3f} for trial type "
    f"and chi2 = {REPORTED['trans_session_chi2']}, df = 1, p = {REPORTED['trans_session_p']:.3f} for "
    f"session number; the subject-level recomputation agrees in conclusion for trial type"))

claims = [c1, c2, c3, c5]
P_REPORTED = {"C1": REPORTED["pupil_wilcox_p"], "C2": REPORTED["dwell_trialtype_p"],
              "C3": REPORTED["trans_trialtype_p"], "C5": np.nan}
for c in claims:
    row = c["ledger_row"]
    row.loc[0, "verdict"] = c["verdict"]
    row.loc[0, "reasons"] = " | ".join(c["reasons"])
    row.loc[0, "p_reported"] = P_REPORTED[c["claim_id"]]

# C5's ledger row carries the interaction - the quantity the claim asserts and the
# one the paper never tests - not the single per-gesture contrast review_claim was
# built on. p_reported stays blank because no test of it is reported anywhere.
c5["ledger_row"].loc[0, "effect_recomputed"] = INTER["effect"]
c5["ledger_row"].loc[0, "ci_low"] = INTER["boot"]["ci_low"]
c5["ledger_row"].loc[0, "ci_high"] = INTER["boot"]["ci_high"]
c5["ledger_row"].loc[0, "p_recomputed"] = float(INTER["perm"]["p_value"])
c5["ledger_row"].loc[0, "claim_short"] = "effect is context-specific (action-balanced design)"

hand = pd.DataFrame([
    {"claim_id": "C4", "claim_short": "eye-tracking study with N = 39 dogs",
     "verdict": "supported-with-caveat",
     "effect_reported": 39.0, "effect_recomputed": float(NBOOK["n_subject_parquet_ds"]),
     "ci_low": np.nan, "ci_high": np.nan, "p_reported": np.nan, "p_recomputed": np.nan,
     "reasons": (
         f"N = 39 reproduces: 39 distinct subjects in both pupil parquet files and 39 named dogs in "
         f"the counterbalancing file | the sample description reproduces too - the paper reports "
         f"{REPORTED['n_female']} females and a mean age of {REPORTED['mean_age_months']} months, "
         f"the file gives {DEMO['n_female']} (match: {DEMO['female_matches']}) and "
         f"{DEMO['mean_age']:.2f} months (match to 1 dp: {DEMO['age_matches_1dp']}) | the paper "
         f"states N separately for every analysis (N = {REPORTED['n_pupil']} for pupil size, "
         f"N = {REPORTED['n_total']} for the all-data re-analysis), so the differing n per analysis "
         f"is disclosed, not a gap | the counterbalancing file's {NBOOK['n_rows_cb']} rows and "
         f"{NBOOK['n_unique_dogid']} unique Dog_ID values are likewise not a subject-count "
         f"discrepancy: {NBOOK['n_blank_rows']} rows are blank and every repeated Dog_ID maps to "
         f"exactly one counterbalancing order string (order unique within Dog_ID: "
         f"{NBOOK['dogid_is_order_slot']}, over {NBOOK['n_shared_dogid']} shared ids), so Dog_ID "
         f"indexes the order slot and Name is the subject key | ONE ARITHMETIC GAP REMAINS: the "
         f"methods state {REPORTED['n_excluded_trials']} trials were excluded from the pupil "
         f"analyses ({REPORTED['n_excluded_pupil_lost']} for pupil detected under 70% of the IP, "
         f"{REPORTED['n_excluded_offscreen']} for gaze on-screen under 70%), but the deposited files "
         f"differ by {TRIALS['dropped']}: ALL_DATA holds {TRIALS['alldata']} trials and the analysis "
         f"file holds {TRIALS['ds']}, and only {TRIALS['dropped_with_pupil']} of the removed trials "
         f"carried pupil data to exclude ({', '.join(TRIALS['empty_trials'])} is empty across all 41 "
         f"bins), so {TRIALS['gap_vs_text']} trial is unaccounted for on the loosest reading and "
         f"{TRIALS['gap_vs_text_strict']} on the strictest; no verdict depends on it because all "
         f"{TRIALS['dropped']} belong to dogs already outside the 28-dog sample | the "
         f"{TRIALS['intended'] - TRIALS['alldata']} trials missing from the "
         f"{TRIALS['intended']} intended match the methods exactly (one dog 3 missing sessions, one "
         f"dog 1)")},
    {"claim_id": "C6", "claim_short": "expectations are not self-centred (third party)",
     "verdict": "not-assessable",
     "effect_reported": np.nan, "effect_recomputed": np.nan, "ci_low": np.nan, "ci_high": np.nan,
     "p_reported": np.nan, "p_recomputed": np.nan,
     "reasons": ("the full text confirms the first pass's reading: the paper argues this from the "
                 "stimulus design and not from any test - 'dogs' pupil sizes reacted to the "
                 "incongruent actions despite the dogs never being addressed by the depicted "
                 "agents' - so there is no statistic to audit | the data confirm the structural "
                 "premise - every trial carries two agent interest "
                 "areas, one communicating (inviting/pointing) and one reacting (approaching/leaving), "
                 "so the dog is not the addressee - but whether the resulting expectation is "
                 "'not self-centred' is an inference from the stimulus design, not a statistic | a "
                 "first-person control condition, in which the gesture is addressed to the dog, "
                 "would change this verdict")},
    {"claim_id": "C7", "claim_short": "expectations form implicitly during ontogeny",
     "verdict": "not-assessable",
     "effect_reported": np.nan, "effect_recomputed": np.nan, "ci_low": np.nan, "ci_high": np.nan,
     "p_reported": np.nan, "p_recomputed": np.nan,
     "reasons": ("the full text confirms the first pass's reading: the paper frames this as a "
                 "hypothesis ('we hypothesised that, due to the experience they gain during their "
                 "ontogeny in a human household') and as a possibility for discussion, never as a "
                 "tested result | the data contain no developmental or experience measure | age in months is "
                 "recorded in the counterbalancing file but no claim links it to the pupil effect; "
                 "a rearing-history or longitudinal measure would change this verdict")},
])

ledger = pd.concat([K.review_summary_table(claims), hand], ignore_index=True)
ledger = ledger.sort_values("claim_id").reset_index(drop=True)
ledger.to_csv(os.path.join(OUT, "claim_ledger.csv"), index=False)
counts = ledger["verdict"].value_counts().reindex(rf.VERDICT_ORDER, fill_value=0)

# statcheck.csv was written above, from the paper's own reported triples.

# ----------------------------------------------------------------------------
# 9. Multiplicity: the test inventory is counted from the authors' code, because
#    no reported p-values exist to form a distribution from.
# ----------------------------------------------------------------------------
TEST_INVENTORY = [
    ("pupil GAMM, complete cases (dat3): trial_type term + s(time, by=trial_type)", 2),
    ("pupil GAMM, ALL_DATA (dat2): trial_type term + s(time, by=trial_type)", 2),
    ("pupil, subject-aggregated paired t-test", 1),
    ("pupil, Shapiro-Wilk on the difference scores", 1),
    ("pupil, exact Wilcoxon signed-rank", 1),
    ("dwell proportion beta GLMM, trial_type x agent_type interaction (drop1)", 1),
    ("dwell proportion beta GLMM, main effects (drop1: trial_type, agent_type, z.trial)", 3),
    ("dwell time, reacting vs communicating agent, paired t-test", 1),
    ("dwell proportion beta GLMM, condition x agent_type vs reduced model (LR test)", 1),
    ("transition count Poisson GLMM (drop1: trial_type, z.trial)", 2),
    ("transition count, zero-inflation simulation", 1),
    ("overdispersion tests (two models)", 2),
]
M_TESTS = int(sum(n for _, n in TEST_INVENTORY))
recomputed_p = {
    "C1 pupil, 28 dogs (permutation)": A["p_perm"],
    "C1 pupil, 28 dogs (cluster permutation)": CP28["p_cluster"],
    "C1 pupil, 38 dogs ALL_DATA (permutation)": B["p_perm"],
    "C2 looking time, ms (permutation)": DWL["p_perm"],
    "C2 looking time, proportion (permutation)": R[KEYPROP]["p_perm"],
    "C3 gaze transitions (permutation)": TR["p_perm"],
    "C5 leaving contrast (permutation)": LEA["p_perm"],
    "C5 approach contrast (permutation)": APP["p_perm"],
    "C5 interaction (permutation)": float(INTER["perm"]["p_value"]),
}
# The paper's own reported p-values; bounds (p < .001) are recorded at the bound,
# which is the most favourable reading for the paper and is stated as such.
REPORTED_P = audit["pvalues"]
recount_paper = mult.discovery_recount(REPORTED_P, alpha=ALPHA)
rp_sum = recount_paper["summary"]
recount_paper_md = "\n".join(
    f"| {r['correction']} | {r['controls']} | {int(r['n_significant'])} | {r['threshold_on_raw_p']:.4f} |"
    for _, r in rp_sum.iterrows())
JUSTBELOW = audit["just_below"]

recount = mult.discovery_recount(list(recomputed_p.values()), alpha=ALPHA,
                                 labels=list(recomputed_p.keys()))
rc_tab, rc_sum = recount["table"], recount["summary"]
SURV_BH = rc_tab.loc[rc_tab["significant_bh"], "label"].tolist()
SURV_BONF = rc_tab.loc[rc_tab["significant_bonferroni"], "label"].tolist()
SURV_RAW = rc_tab.loc[rc_tab["significant_uncorrected"], "label"].tolist()
recount_md = "\n".join(
    f"| {r['correction']} | {r['controls']} | {int(r['n_significant'])} | {r['threshold_on_raw_p']:.4f} |"
    for _, r in rc_sum.iterrows())

# ----------------------------------------------------------------------------
# 10. Figures
# ----------------------------------------------------------------------------
PUP_SHORT = "pupil change (a.u.)"
trial_level = (pup_all.groupby(["subject", "condition", "trial_type"], observed=True)["pupil_base"]
               .mean().reset_index().rename(columns={"pupil_base": PUP_SHORT}))
trial_level["subject"] = trial_level["subject"].astype(str)
trial_level["trial_type"] = trial_level["trial_type"].astype(str)
indep = ac.design_independence(trial_level, "subject", "trial_type", PUP_SHORT)

figs = {}
rf.fig_design_map(trial_level, "subject", "trial_type", PUP_SHORT, independence=indep,
                  path=os.path.join(FIG, "F1_design_map.png"), unit_label="dog")
figs["F1_design_map"] = os.path.join(FIG, "F1_design_map.png")

rf.fig_raw_behind_claim(c1["groups"], PUP_LAB, path=os.path.join(FIG, "F2_raw_behind_claim.png"),
                        test_name="sign-flip permutation (mean of paired differences)",
                        p_value=A["p_perm"],
                        reported_error=("SD (the paper's '+/-' resolves to SD, ratio "
                                        f"{ERRBAR['reacting']['ratio_to_candidates']['SD']:.3f})"))
figs["F2_raw_behind_claim"] = os.path.join(FIG, "F2_raw_behind_claim.png")

rf.fig_assumption_qq(c1["groups"], path=os.path.join(FIG, "F3_assumption_shape.png"),
                     paired_differences=c1["differences"], value_label=PUP_LAB)
figs["F3_assumption_shape"] = os.path.join(FIG, "F3_assumption_shape.png")

# C4's "effect" is a sample size, not a contrast in the pupil's units; plotting 39
# on the same axis as an a.u. difference would misread the axis, so it is blanked
# in the plotted copy only. The real values stay in claim_ledger.csv.
ledger_plot = ledger.copy()
ledger_plot.loc[ledger_plot["claim_id"] == "C4",
                ["effect_reported", "effect_recomputed", "ci_low", "ci_high"]] = np.nan
rf.fig_claim_ledger(ledger_plot, path=os.path.join(FIG, "F4_claim_ledger.png"),
                    effect_label="effect (incongruent - congruent, paper's units)")
figs["F4_claim_ledger"] = os.path.join(FIG, "F4_claim_ledger.png")

rf.fig_null_distribution(c1["permutation"], c1["bootstrap"],
                         path=os.path.join(FIG, "F5_null_distribution.png"),
                         statistic_label="mean paired difference in pupil size (a.u.)")
figs["F5_null_distribution"] = os.path.join(FIG, "F5_null_distribution.png")

# F6 is assessable on this pass: the paper's own reported p-values are extractable.
rf.fig_multiplicity(REPORTED_P, alpha=ALPHA, path=os.path.join(FIG, "F6_multiplicity.png"))
figs["F6_multiplicity"] = os.path.join(FIG, "F6_multiplicity.png")

# re-render each per-claim panel with the reviewer's final verdict: review_claim
# stamped the automatic proposal, which was overridden in section 8.
for c in claims:
    if c["claim_id"] == "C5":
        # C5's panel is the interaction its ledger row carries: each dog's two
        # gesture-wise congruency contrasts, and the null for their difference.
        rf.fig_claim_panel("C5", {"after 'go away' (GA - CA)": _dA,
                                  "after 'come here' (CG - GG)": _dL},
                           "congruency contrast in pupil size (a.u.)",
                           INTER["perm"], INTER["boot"], INTER["routes"],
                           verdict=c["verdict"], path=c["figure_path"])
    else:
        rf.fig_claim_panel(c["claim_id"], c["groups"],
                           "pupil size change (a.u.)" if c["claim_id"] == "C1"
                           else (PROP_LAB if c["claim_id"] == "C2" else TR_LAB),
                           c["permutation"], c["bootstrap"], c["routes"],
                           verdict=c["verdict"], path=c["figure_path"])
    figs[f"panel_{c['claim_id']}"] = c["figure_path"]


def fig_timecourse(cp28, cpall, path):
    """F7: the time-resolved contrast, the estimand the authors' GAMM targets."""
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    for ax, cp, title in ((axes[0], cp28, "authors' primary sample (complete cases)"),
                          (axes[1], cpall, "ALL_DATA (the abstract's N = 39)")):
        tv = cp["times"] + 5000.0
        mean = cp["X"].mean(0)
        sem = cp["X"].std(0, ddof=1) / np.sqrt(cp["n"])
        for (t0, t1, _m) in cp["windows"]:
            ax.axvspan(t0 + 5000, t1 + 5000, color="#E69F00", alpha=0.30, linewidth=0)
        ax.axhline(0, color="#444444", lw=0.8)
        ax.fill_between(tv, mean - 1.96 * sem, mean + 1.96 * sem, color="#0072B2", alpha=0.20,
                        linewidth=0)
        ax.plot(tv, mean, color="#0072B2", lw=1.8)
        ax.set_title(f"{title}\nn = {cp['n']} dogs  ·  cluster p = {cp['p_cluster']:.4f}"
                     f"  (resolution {cp['p_resolution']:.0e})", fontsize=9.5)
        ax.set_xlabel("video time (ms)")
        ax.tick_params(direction="out")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("incongruent - congruent pupil size\n(a.u., mean of dogs +/- 95% CI)")
    axes[0].text(0.02, 0.04, f"shaded: bins with |t| > {cp28['thr']:.2f}\n"
                             f"{len(cp28['windows'])} candidate cluster(s); largest "
                             f"{CW0:.0f}-{CW1:.0f} ms",
                 transform=axes[0].transAxes, fontsize=7.5, color="#555555", va="bottom")
    axes[1].text(0.02, 0.04,
                 (f"no bin reaches |t| > {cpall['thr']:.2f}" if not cpall["windows"]
                  else f"shaded: bins with |t| > {cpall['thr']:.2f}"),
                 transform=axes[1].transAxes, fontsize=7.5, color="#555555", va="bottom")
    fig.suptitle("F7 - Time-resolved pupil contrast, cluster-based sign-flip permutation "
                 "(subject is the unit)", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=rf.FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    return path


figs["F7_pupil_timecourse"] = fig_timecourse(CP28, CPALL,
                                             os.path.join(FIG, "F7_pupil_timecourse.png"))

# ----------------------------------------------------------------------------
# 11. Report
# ----------------------------------------------------------------------------
def routes_md(block):
    out = []
    for _, r in block["routes"].iterrows():
        out.append(f"| {r['route']} | {r['test']} | {r['p_value']:.4f} | {r['assumes']} | "
                   f"{r['min_achievable_p']:.2e} | "
                   f"{'yes' if bool(r['reachable_at_0.05']) else 'no'} |")
    return "\n".join(out)


def ledger_md():
    out = ["| id | claim | verdict | effect rep. | effect recomp. | 95% CI | p rep. | p recomp. |",
           "|----|-------|---------|-------------|----------------|--------|--------|-----------|"]
    for r in ledger.itertuples():
        ci = "-" if not np.isfinite(r.ci_low) else f"{r.ci_low:.3g} to {r.ci_high:.3g}"
        er = "-" if not np.isfinite(r.effect_reported) else f"{r.effect_reported:.4g}"
        ec = "-" if not np.isfinite(r.effect_recomputed) else f"{r.effect_recomputed:.4g}"
        pr = "-" if not np.isfinite(r.p_reported) else f"{r.p_reported:.4g}"
        pc = "-" if not np.isfinite(r.p_recomputed) else f"{r.p_recomputed:.4g}"
        out.append(f"| {r.claim_id} | {r.claim_short} | `{r.verdict}` | {er} | {ec} | {ci} | {pr} | {pc} |")
    return "\n".join(out)


ALL_CONTRASTS = list(R.items()) + list(CTX.items())
contrast_rows = "\n".join(
    f"| {k} | {v['n']} | {v['effect']:.4g} | {v['ci'][0]:.4g} to {v['ci'][1]:.4g} | {v['p_t']:.4f} | "
    f"{v['p_rank']:.4f} | {v['p_sign']:.4f} | {v['p_perm']:.4f} | {'yes' if v['routes_agree'] else 'NO'} |"
    for k, v in ALL_CONTRASTS)
floor_rows = "\n".join(
    f"| {k} | {v['n']} | {v['p_perm']:.4f} | {v['min_p']:.2e} | {v['sign_floor']:.2e} | "
    f"{v['p_resolution']:.0e} | {'yes' if v['perm']['floor_binds'] else 'no'} |"
    for k, v in ALL_CONTRASTS)

checks = {
    "pupil sample report (ds)": (P_SAMPLE, pup_ds.shape),
    "pupil sample report (ALL_DATA)": (P_ALL, pup_all.shape),
    "IA report, reaction IP": (P_IA, ia.shape),
    "IA report, reaction IP from frame 579": (P_IA579, ia579.shape),
    "IA report, whole-video IP (unused by both scripts)": (P_IANEW, ianew.shape),
    "counterbalancing / demographics": (P_CB, cb.shape),
    "pupil analysis code": (P_RMD_PUP, None),
    "dwell-time analysis code": (P_RMD_DWL, None),
}
inputs_md = "\n".join(
    f"| {k} | `{os.path.basename(p)}` | {('%d x %d' % s) if s else 'R Markdown'} | `{sha256(p)}` |"
    for k, (p, s) in checks.items())

report = f"""# Statistical Review — *Lonardo, Volter & Huber (2026), Cognition 274:106586*

*Reviewer: Statistical Reviewer (automated audit) · Date: {datetime.date.today().isoformat()} · Protocol: statistical-review v1*

**Paper.** Lonardo L., Volter C.J. & Huber L. (2026) Pupil size changes reveal that dogs are
sensitive to the social conventions behind human gestures. *Cognition* **274**:106586.
doi:10.1016/j.cognition.2026.106586

---

## 1. Inputs and coverage

**This review was done in two passes.** The first pass had only the abstract (Elsevier paywall);
the published PDF was supplied afterwards, and this version is the revision against it. The first
pass's limitations are gone: statcheck now runs, `p_reported` is populated across the ledger, and
F4 and F6 have both axes.

**What the full text confirmed, and what it overturned.** Two of the five abstract-based readings
did not survive contact with the paper, and both errors were the reviewer's:

| First-pass reading | Full text | Outcome |
|--------------------|-----------|---------|
| **C1** — the 28-vs-39 sample split is an undisclosed dependence; `supported-with-caveat` | Methods: "we only retained dogs with four valid test trials, one per condition (N = 28)"; Results report the N = 39 re-analysis and its null (t = {REPORTED['pupil_n39_t']}, p = {REPORTED['pupil_n39_p']:.3f}) | **Overturned.** No disclosure gap. C1 is now a reproduction result and the verdict is `supported` |
| **C5** — "context-specific" rests on a difference of significance between the two gestures; `assumption-dependent` | The paper never splits the pupil analysis by gesture. "Context-relative" is an argument from the design: the same action is congruent after one gesture and incongruent after the other, so the action cancels in the trial-type contrast | **Overturned.** The reviewer read a per-gesture significance split into a sentence that restates the design. Verdict is now `supported-with-caveat` |
| **C2, C3** — the nulls reproduce but the intervals do not establish equivalence | Reported as chi2 = {REPORTED['dwell_trialtype_chi2']}, df = 1, p = {REPORTED['dwell_trialtype_p']:.3f} and chi2 < 0.01, df = 1, p = {REPORTED['trans_trialtype_p']:.3f}; no equivalence test, no smallest effect of interest | **Confirmed**, and now evidenced from the text rather than inferred |
| **C6, C7** — argued from design and offered as discussion, with no test to audit | Confirmed verbatim in the discussion and hypothesis sections | **Confirmed**; both stay `not-assessable` |
| **C4** — N = 39 reproduces; the per-analysis n differ | The paper states n separately for every analysis, and its demographics reproduce exactly | **Mostly confirmed**; one trial-accounting gap remains (section 4) |

**How well the blind recomputation did.** The first pass derived its numbers from the raw data
without seeing any of the paper's. The subject-level Wilcoxon reproduces the paper's reported
statistic exactly — reported T+ = {REPORTED['pupil_wilcox_Tplus']:.0f}, N = {REPORTED['pupil_wilcox_n']},
p = {REPORTED['pupil_wilcox_p']:.3f} against recomputed T+ = {REPRO['wilcox_Tplus']:.0f},
N = {REPRO['wilcox_n']}, p = {REPRO['wilcox_p']:.4f} — and the cluster-based permutation test,
run as a stand-in for a GAMM that could not be refitted, placed the difference at
{CW0:.0f}–{CW1:.0f} ms against the paper's reported significant window of
{REPORTED['curve_lo']:.0f}–{REPORTED['curve_hi']:.0f} ms. The reviewer's errors were in reading the
prose, not in the arithmetic.

| Input | File | Size | sha256 (first 16) |
|-------|------|------|-------------------|
{inputs_md}

| Input | File | Size | sha256 (first 16) |
|-------|------|------|-------------------|
| published article | `1-s2.0-S0010027726001538-main.pdf` | {len(PAGES)} pages, text layer usable ({quality['usable_text_layer']}) | `{sha256(PDF)}` |

| | |
|--|--|
| Claims reviewed | 7 (C1–C7) |
| Extracted automatically from the PDF | {N_AUTO_T} (statistic, p) pairs, {N_AUTO_P} p-values, {N_AUTO_N} sample-size statements, {len(inventory['intervals'])} intervals |
| Added by hand | 6 chi-square triples the extractor's regex misses (Greek glyph), transcribed from the sentences its p-value pass did capture; and the 7 claims themselves, since the paper states them in prose |
| Claims the data cannot address | C6, C7 — `not-assessable` |
| Seeds | `random_state = {SEED}`, `n_perm = {N_PERM}`, `n_boot = {N_BOOT}`, alpha = {ALPHA} |

**Coverage, stated honestly.** This reviews the seven claims the abstract makes, now checked
against the full text and the deposited data. {M_TESTS} inferential tests are visible in the
authors' code and {len(REPORTED_P)} p-values are reported in the main text; this review re-derives
five of them from raw data and arithmetic-checks {N_SC_CHECKED} of the reported triples. What
remains outside it: **Tables S1–S4 and the supplementary figures were not supplied**, so the GAMM
and GLMM coefficient tables, their degrees of freedom, and the 1000-sample bootstrap CIs the
authors computed are unseen — that is why {N_SC_NOTCHECK} of the extracted test statistics could
not be arithmetic-checked (section 5). The authors' three model classes could not be refitted
here, so the published fits themselves are not reproduced; this review audits the contrasts those
fits estimate.

**Which data correspond to which claim.**

| Claim | Supporting file |
|-------|-----------------|
| C1, C5 (pupil size) | `..._samplereport_ds.parquet` ({pup_ds.shape[0]} rows) and `..._ALL_DATA.parquet` ({pup_all.shape[0]} rows) — 41 time bins per trial over 1000–5000 ms of the baseline+reaction interest period (video time 6000–10000 ms) |
| C2 (looking time) | `IA_report_IP_from_reaction.txt` plus `..._frame579.txt` for Zserbo session 1 — {len(dw)} agent interest-area rows from {dw['RECORDING_SESSION_LABEL'].nunique()} sessions |
| C3 (gaze transitions) | the same IA reports, via `INTEREST_AREA_FIXATION_SEQUENCE` — {len(trd)} trial-level transition counts |
| C4 (N = 39) | both parquet files plus `Non_verbal_communication_counterbalancing.csv` |
| C6, C7 | no corresponding data |

`IA_report_new_video_IP.txt` covers the whole-video interest period and is used by neither of the
authors' scripts; it is not used here either.

**The models the authors actually specified**, read from the deposited `.Rmd` files. The task
brief described the analysis as "GLMMs with bootstrapped CIs"; that is accurate for looking time
and gaze transitions but not for the headline pupil analysis:

- **Pupil size (headline): a GAMM, not a GLMM** — `bam(pupil_base ~ trial_type + session_number +
  s(time, k=20) + s(time, by=trial_type, k=20) + s(Xgaze_RIGHT, Ygaze_RIGHT) +
  s(time, subject, bs='fs', m=1) + s(time, Event, bs='fs', m=1), method="ML")`, fitted twice: on
  complete cases (`dat3`, {NBOOK['n_primary_pupil']} dogs) and on `ALL_DATA`
  ({NBOOK['n_subject_parquet_all']} dogs), with inference read off the `plot_diff` difference
  curve. The same script additionally runs a subject-aggregated paired t-test, a Shapiro-Wilk test
  on the difference scores, and an exact Wilcoxon signed-rank test on
  {NBOOK['n_primary_pupil']} dogs.
- **Looking time:** `glmmTMB(prop.dwell.time.agents.scaled ~ trial_type + agent_type + z.trial +
  (1 + trial_type.c + agent_type.c + z.trial || subject), family = beta_family)`, evaluated by
  `drop1(test="Chisq")`, with 1000 parametric bootstrap CIs.
- **Gaze transitions:** `glmer(transition_count ~ trial_type + z.trial +
  offset(log(TRIAL_DWELL_TIME)) + (1 + z.trial | subject), family = poisson)`, `drop1(test="Chisq")`,
  1000 bootstrap CIs.

None of these three model classes is available in this environment, so **this review does not
reproduce the authors' models.** It re-derives the same contrasts with dog as the unit of analysis
and, for pupil size, adds a time-resolved cluster-based permutation test, which targets the same
estimand as the GAMM difference curve — a difference localised in time — without assuming a
parametric smooth.

---

## 2. Summary of findings

| Verdict | Claims |
|---------|--------|
| supported | {int(counts['supported'])} |
| supported-with-caveat | {int(counts['supported-with-caveat'])} |
| assumption-dependent | {int(counts['assumption-dependent'])} |
| underpowered-resolution | {int(counts['underpowered-resolution'])} |
| not-reproducible | {int(counts['not-reproducible'])} |
| contradicted | {int(counts['contradicted'])} |
| not-assessable | {int(counts['not-assessable'])} |

1. **C1 reproduces exactly, and the first pass's caveat against it was the reviewer's error.** The
   paper reports T+ = {REPORTED['pupil_wilcox_Tplus']:.0f}, N = {REPORTED['pupil_wilcox_n']},
   p = {REPORTED['pupil_wilcox_p']:.3f}; recomputing from the raw parquet gives
   T+ = {REPRO['wilcox_Tplus']:.0f}, N = {REPRO['wilcox_n']}, p = {REPRO['wilcox_p']:.4f}, with all
   four routes agreeing (paired t {A['p_t']:.4f}, Wilcoxon {A['p_rank']:.4f}, sign
   {A['p_sign']:.4f}, permutation {A['p_perm']:.4f}) and a blind cluster permutation landing at
   {CW0:.0f}–{CW1:.0f} ms inside the paper's reported {REPORTED['curve_lo']:.0f}–{REPORTED['curve_hi']:.0f} ms
   window. The 28-vs-39 split the first pass flagged is stated in the methods and the N = 39 null
   is reported by the authors themselves (t = {REPORTED['pupil_n39_t']},
   p = {REPORTED['pupil_n39_p']:.3f}); verdict moves to `supported` (stages 4–5).
2. **One to two excluded trials are unaccounted for between the methods and the deposit.** The
   methods state {REPORTED['n_excluded_trials']} trials were excluded from the pupil analyses
   ({REPORTED['n_excluded_pupil_lost']} for pupil detection below 70% of the interest period,
   {REPORTED['n_excluded_offscreen']} for on-screen gaze below 70%); the deposited files differ by
   {TRIALS['dropped']} trials — `ALL_DATA` holds {TRIALS['alldata']} and the analysis file
   {TRIALS['ds']} — and only {TRIALS['dropped_with_pupil']} of those carried pupil data to exclude,
   since {', '.join(TRIALS['empty_trials'])} is empty across all 41 bins. No verdict depends on it:
   all {TRIALS['dropped']} belong to dogs already outside the 28-dog sample. Everything else in the
   accounting reconciles to the digit: 39 subjects,
   {DEMO['n_female']} females against {REPORTED['n_female']} reported, mean age
   {DEMO['mean_age']:.2f} against {REPORTED['mean_age_months']} months, and the
   {TRIALS['intended'] - TRIALS['alldata']} missing sessions exactly as described (stage 1).
3. **Both nulls reproduce; neither interval establishes equivalence, and no correction is applied
   across {len(REPORTED_P)} reported p-values.** Looking time differs by {DWL['effect']:.1f} ms
   (95% CI {DWL['ci'][0]:.0f} to {DWL['ci'][1]:.0f}, on a congruent mean of
   {DWL['mean_cong']:.0f}) and transitions by {TR['effect']:.3f}
   (95% CI {TR['ci'][0]:.2f} to {TR['ci'][1]:.2f}, mean {TR['mean_cong']:.2f}) — roughly
   ±{max(abs(DWL['ci'][0]), abs(DWL['ci'][1])) / DWL['mean_cong'] * 100:.0f}% and
   ±{max(abs(TR['ci'][0]), abs(TR['ci'][1])) / TR['mean_cong'] * 100:.0f}% of their means remain
   compatible with the data, and the paper states no smallest effect of interest to test against
   (stages 5–6).

---

## 3. Claim ledger — F4

![F4 — claim ledger](figures/F4_claim_ledger.png)

*Recomputed effect with its 95% bootstrap CI (left) and reported against recomputed p-values on
log axes with alpha marked (right). Effects are signed incongruent − congruent, so a positive
pupil effect is the direction the paper claims. C1, C2 and C3 now carry a reported p-value; C5 does
not, because the paper reports no test of the quantity that claim asserts (section 5), and C4's
"effect" is a sample size, which is held out of the effect axis and verified in section 4 instead.*

{ledger_md()}

The full `reasons` text for every row is in `claim_ledger.csv`.

---

## 4. Design and sampling — stages 1–2

![F1 — design map](figures/F1_design_map.png)

*Trials per dog and dogs per trial type, from `ALL_DATA`. The n available to a test is dogs — not
trials, and not the {pup_all.shape[0]} time-bin rows. The panel counts
{TRIALS['alldata_with_pupil']} trials rather than the {TRIALS['alldata']} present in the file
because one trial ({', '.join(TRIALS['empty_trials'])}) has no pupil value in any of its 41 bins
and drops out of a figure keyed on that value; section 4(d) follows that trial through the
accounting.*

![F2 — raw data behind the claim](figures/F2_raw_behind_claim.png)

*Each dog's window-mean pupil size by trial type, with SD and SEM shown on one axis. The paper's
own error-bar definition could not be checked, because its text was not supplied.*

**Population.** The abstract claims a property of pet dogs in general ("pet dogs, during their
ontogeny"). The deposited data describe {NBOOK['n_named']} dogs tested in a laboratory eye-tracker,
{NBOOK['n_female']} of them female, carrying {NBOOK['n_breeds']} distinct breed labels. Whether
that sample represents pet dogs as a population is a methods question this review cannot answer:
recruitment frame, screening and the chin-rest training requirement are described only in the
paper. Chin-rest eye-tracking selects for trainable, cooperative dogs, and nothing in the data
bounds that selection.

**Unit of replication.** The pupil files hold {pup_ds.shape[0]} and {pup_all.shape[0]} rows, but
these are 41 time bins × up to 4 trials × dog, so the row count is not an n. Two separate things
determine what n a test may use, and they point the same way here for different reasons.

The first is clustering, which turns out to be slight: across the {trial_level.shape[0]} trial-level
pupil means from {NBOOK['n_subject_parquet_all']} dogs, the ICC within dog is {indep['icc']:.3f}
with a design effect of {indep['design_effect']:.2f}, giving an effective n of
{indep['effective_n']:.1f} — essentially the row count. A dog's four trial means are, on this
summary, about as informative as four independent observations, so a variance-inflation argument
against pooling trials is weak. (Within a trial, the 41 time bins are a different matter: they are
a smoothed time series, which is why the authors model them with `s(time, Event)` rather than
counting them.)

The second is the design, and it is the binding one: congruency was assigned **within** dog, so the
comparison each claim makes is a paired one and the number of independent contrasts is the number of
dogs — {NBOOK['n_subject_parquet_all']} at most, fewer wherever a dog is missing a condition. A near-zero
ICC does not license treating the trials as independent replicates of the contrast, because it is
the pairing, not the correlation, that defines the unit. Every test in this review therefore takes
dog as the unit by construction; the authors' GAMM reaches the same place through its
`s(time, subject)` and `s(time, Event)` random smooths.

**Who actually contributes.** The review's first substantive finding concerns a count, not a
p-value:

| Analysis | Dogs contributing |
|----------|-------------------|
| Abstract's stated N | 39 |
| Distinct subjects in either parquet file | {NBOOK['n_subject_parquet_ds']} |
| Named dogs in the counterbalancing file | {NBOOK['n_named']} (of {NBOOK['n_rows_cb']} rows; {NBOOK['n_blank_rows']} blank) |
| Authors' primary pupil GAMM and aggregated tests (`dat3`) | {NBOOK['n_primary_pupil']} |
| Pupil contrast available at subject level, `ds` file | {C['n']} |
| Pupil contrast available at subject level, `ALL_DATA` | {B['n']} |
| Looking time / transitions | {int(trd['subject'].nunique())} dogs, {len(trd)} trials; {W[KEYDWL].shape[0]} with a paired contrast |

Three things follow.

**(a) N = 39 reproduces exactly** — 39 distinct subjects in both parquet files, 39 named dogs in
the counterbalancing file.

**(b) The counterbalancing file's {NBOOK['n_rows_cb']} rows and {NBOOK['n_unique_dogid']} unique
`Dog_ID` values are not a subject-count discrepancy**, and this review checked rather than assumed
it. {NBOOK['n_blank_rows']} rows are blank, and every `Dog_ID` appearing more than once maps to
exactly one counterbalancing `order` string (`order` unique within `Dog_ID`:
{NBOOK['dogid_is_order_slot']}, over {NBOOK['n_shared_dogid']} shared ids; `Dog_ID` 18, for
instance, is shared by three dogs that all ran `GG-GA-CA-CG`). `Dog_ID` indexes the
counterbalancing slot; `Name` is the subject key. The apparent 45/30 mismatch is a column-semantics
question, not a missing-dogs question, and it is resolved by that one sentence.

**(c) The analyses do not all run on 39.** Eleven dogs ({', '.join(EXCL)}) lack at least one of the
four conditions. The authors exclude exactly those eleven from the primary pupil analysis, leaving
{NBOOK['n_primary_pupil']}; this review confirms the hard-coded exclusion list in their script is
precisely the set of dogs with fewer than four conditions in the data, so the rule is principled
and reproducible. **The paper states this**: "To ensure a within-subject comparison, we only
retained dogs with four valid test trials, one per condition (N = 28)", and it reports the N = 39
re-analysis and its null result as well. The abstract's "N = 39" is the recruited and eye-tracked
sample and the pupil result rests on {A['n']} dogs — a reader of the abstract alone could conflate
them, as the first pass of this review did, but the paper itself does not.

**(d) Trial accounting, which is where the one remaining gap is.** The two pupil files nest
exactly — every trial in the analysis file is also in `ALL_DATA` (`subset = {TRIALS['ds_is_subset']}`)
— so the chain can be followed step by step:

| step | trials | source |
|------|--------|--------|
| intended: 39 dogs × 4 conditions | {TRIALS['intended']} | design |
| less sessions never recorded | −{TRIALS['missing_sessions']} | methods: "one dog had 3 missing sessions – CC, CG, GC; the other dog had 1 missing session - GC" |
| present in `ALL_DATA` | {TRIALS['alldata']} | {TRIALS['alldata']} × 41 time bins = {TRIALS['alldata'] * 41} rows, its exact row count |
| of which carry no pupil data at all | {len(TRIALS['empty_trials'])} ({', '.join(TRIALS['empty_trials'])}: all 41 bins missing) | data |
| present in the analysis file | {TRIALS['ds']} | {TRIALS['ds']} × 41 = {TRIALS['ds'] * 41} rows |
| removed between the two files | {TRIALS['dropped']} | set difference |
| of those, carrying usable pupil data | {TRIALS['dropped_with_pupil']} | data |
| excluded for quality per the methods | {REPORTED['n_excluded_trials']} | {REPORTED['n_excluded_pupil_lost']} pupil detected under 70% of the IP + {REPORTED['n_excluded_offscreen']} on-screen gaze under 70% |

Everything above the last line reconciles. The last line does not: the methods describe
{REPORTED['n_excluded_trials']} quality exclusions, while {TRIALS['dropped']} trials are absent
from the analysis file, of which only {TRIALS['dropped_with_pupil']} carried any pupil data to
exclude. **So {TRIALS['gap_vs_text']} trial is unaccounted for on the loosest reading and
{TRIALS['gap_vs_text_strict']} on the strictest**, depending on whether the empty trial is counted
among the twelve. The removed trials are
{', '.join(TRIALS['dropped_list'])}.

This touches no verdict: all {TRIALS['dropped']} belong to dogs already outside the 28-dog
complete-case sample, which is defined by missing conditions rather than by these exclusions. It
is the kind of discrepancy an author resolves in one sentence, and it is recorded here only
because a reader reconstructing the sample from the deposit will hit it.

A smaller point, now resolved rather than raised: the demographics join in the deposited
dwell-time script recovers the subject name by splitting a truncated EyeLink filename (`EDF3`) and
fails for {len(NBOOK['demo_join_failures'])} dogs ({', '.join(NBOOK['demo_join_failures'])}). The
first pass flagged this as a possible contaminant of the reported age and sex summaries. It is
not: computed over all {DEMO['n_named']} named dogs the file gives {DEMO['n_female']} females and a
mean age of {DEMO['mean_age']:.2f} months, against {REPORTED['n_female']} and
{REPORTED['mean_age_months']} months in the paper — both reproduce — so the published demographics
were evidently not taken through that join.

One naming note for anyone matching the paper to the deposit: the paper labels the four conditions
CC, CG, GC, GG while the data files use CA, CG, GA, GG (A for "approach" where the paper writes C
for "come closer"). The mapping is one-to-one and the congruency coding agrees — CA/CC and GG are
congruent, CG and GA/GC incongruent — but the strings differ.

**Uncertainty reporting — now checkable, and it checks out.** The paper prints
"{REPORTED['reacting_mean']:.0f} ± {REPORTED['reacting_sd']:.0f} ms" and
"{REPORTED['communicating_mean']:.0f} ± {REPORTED['communicating_sd']:.0f} ms" without naming what
the ± is. Recomputing the three candidates from the subject means settles it: for the reacting
agent the SD across dogs is {ERRBAR['reacting']['candidates']['SD']:.1f}, the SEM is
{ERRBAR['reacting']['candidates']['SEM']:.1f} and the 95% CI half-width is
{ERRBAR['reacting']['candidates']['CI95_half_width']:.1f}, so the printed value matches **SD**
(ratio {ERRBAR['reacting']['ratio_to_candidates']['SD']:.4f}, unambiguous:
{ERRBAR['reacting']['unambiguous']}); the communicating agent matches SD on the same test
(ratio {ERRBAR['communicating']['ratio_to_candidates']['SD']:.4f}). SD is the right choice for a
descriptive statement about how dogs vary, and it is the conservative one — SEM would have been
{ERRBAR['reacting']['sd_to_sem_factor']:.2f}× narrower. The only reporting improvement left here
is to name it in the text. F2 draws SD and SEM on one axis so the distinction stays visible.

---

## 5. Method and assumptions — stages 3–4

![F3 — distribution shape](figures/F3_assumption_shape.png)

*Q-Q of the paired differences behind C1 with n, skewness, excess kurtosis and Shapiro-Wilk. This
is descriptive. Shape is not a gate, and nothing on this panel licenses calling any test invalid.*

**statcheck — the paper's own arithmetic.** {len(ALL_TESTS)} test statistics were recovered from
the text ({N_AUTO_T} by the automatic extractor, 6 chi-square triples transcribed by hand because
the extractor's regex does not match the Greek glyph). Of these, **{N_SC_CHECKED} carry the
degrees of freedom needed to recompute p and were checked; {N_SC_NOTCHECK} could not be** — the
GAMM's t and F statistics are reported in the text without their df, which live in Tables S2–S4,
and the supplementary file was not supplied.

Of the {N_SC_CHECKED} checkable triples, statcheck flags **{N_SC_FLAGGED}** and finds
**{int((STATCHECK['decision_error'] == True).sum())} decision errors** (a decision error being a
reported p on the other side of alpha from the recomputed one). Every flagged row resolves once the
precision of the *statistic* is propagated, which the review computed rather than assumed: statcheck
recomputes p from the statistic as printed, so a chi-square given to two decimals pins p far more
tightly than the paper actually did —

""" + "\n".join(
    f"- {v['sentence']}: chi2 = {v['stat']}, df = {v['df1']:.0f}, reported p = {v['p_reported']:.3f}; "
    f"p at the printed value is {v['p_point']:.4f}, but across the statistic's "
    f"{'stated bound [0, ' + format(v['stat'], 'g') + ']' if v['kind'] == 'bound' else 'two-decimal rounding window'}"
    f" the implied p runs {v['lo']:.4f} to {v['hi']:.4f}, which contains the reported value "
    f"({v['resolved']})" for v in ROUNDING.values()) + f"""

So **{N_SC_ROUNDING} of the {N_SC_FLAGGED} flags are artefacts of the tool reading printed
precision as exact, and no genuine arithmetic inconsistency was found** in the statistics this
review could check. The {N_SC_NOTCHECK} unchecked statistics are a gap in coverage, not a finding
about them: supplying Tables S2–S4 would close it.

One transcription decision is worth stating, because getting it wrong would have manufactured a
clean bill of health. For the gaze-transition trial-type test the paper prints "chi2 < 0.01,
df = 1, p = .951": the inequality bounds the *statistic*, while p is reported exactly. Coding that
row's `p_relation` as `<` would route it into statcheck's one-sided branch, which only asks whether
the recomputed p falls below the reported one — a near-automatic pass that would have reported the
row as consistent without ever comparing {ROUNDING.get('HND005', {}).get('p_point', float('nan')):.4f}
to 0.951. It is coded `=` here, which flags it, and the flag is then resolved on the statistic's
stated bound rather than waved through: any chi-square in [0, 0.01] implies p in
[{ROUNDING.get('HND005', {}).get('lo', float('nan')):.4f},
{ROUNDING.get('HND005', {}).get('hi', float('nan')):.4f}], and .951 corresponds to chi2 ≈ 0.0038,
comfortably inside the bound the paper gave.

**Test selection.** On every axis checkable without the methods section the tests match the design:
the comparison is within-dog, every test used is paired or carries a per-subject random effect, and
all are two-sided. The substantive point concerns what a null means rather than which test was
chosen. The authors' pupil script runs a Shapiro-Wilk test on the difference scores, annotates the
result "Assumptions violated? --> Yes", and switches to an exact Wilcoxon signed-rank test. This
review reports the same diagnostic descriptively — Shapiro-Wilk p = {c1['shape']['shapiro_p']:.4f}
on n = {A['n']} paired differences, skewness {c1['shape']['skewness']:.3f}, excess kurtosis
{c1['shape']['excess_kurtosis']:.3f}, median difference {c1['shape']['median_difference']:.1f} a.u.
against a mean of {c1['shape']['mean_difference']:.1f} — and does not treat it as a pass/fail gate.
The shape is right-skewed with heavy tails relative to a normal; that is a description, not a
verdict on any test. At this n the normality test has limited power, so a significant result need
not mean a consequential departure and a non-significant one would not have certified normality.
The reason to report all routes is not that a diagnostic licensed one of them, but that agreement
across them shows the conclusion does not depend on the choice. Note also that the Wilcoxon
signed-rank null is symmetry of the differences about zero rather than equality of means, so the
switch changes the hypothesis slightly; here it changes neither the direction nor the decision.
The skill's own note on this n is carried verbatim rather than paraphrased:
*"{c1['shape']['interpretation_note']}"*

**Three routes — C1, pupil size, {A['n']} dogs** (incongruent − congruent = {A['effect']:.1f} a.u.):

| route | test | p | assumes | smallest achievable p | reachable at 0.05 |
|-------|------|---|---------|-----------------------|-------------------|
{routes_md(A)}

All routes agree at alpha = {ALPHA}. The same four routes across every contrast the data support:

| contrast | n dogs | effect (incong − cong) | 95% bootstrap CI | paired t | Wilcoxon | sign | permutation | routes agree |
|----------|--------|------------------------|------------------|----------|----------|------|-------------|--------------|
{contrast_rows}

Two readings of that table matter.

First, the pupil p-value moves from {A['p_perm']:.4f} ({A['n']} dogs, complete cases) to
{C['p_perm']:.3f} ({C['n']} dogs, same file without the exclusion) to {B['p_perm']:.3f}
({B['n']} dogs, `ALL_DATA`), while the four routes agree with one another *within* every sample.
This is not a modelling-route dependence — the four routes agree with one another *within* every
sample — it is a sample-definition dependence, and **the paper owns it explicitly**: the methods
state the rule and its rationale ("to ensure a within-subject comparison, we only retained dogs
with four valid test trials, one per condition (N = 28)"), and the results report the all-data
re-analysis and its null (t = {REPORTED['pupil_n39_t']}, p = {REPORTED['pupil_n39_p']:.3f}). What
this review adds is the magnitude of the dependence at the subject level, which the paper's own
GAMM t-statistics do not convey: the same contrast moves from p = {A['p_perm']:.4f} to
p = {B['p_perm']:.3f}, and its interval widens from [{A['ci'][0]:.0f}, {A['ci'][1]:.0f}] to
[{B['ci'][0]:.0f}, {B['ci'][1]:.0f}] a.u. Whether the rule was fixed before the outcomes were
examined is still not readable from either the paper or the code — it is the kind of thing a
preregistration answers, and none is cited — but the rule itself is principled, stated, and
reproduces exactly against the data (section 4).

Second, the leaving-response contrast (CG vs GG) is the one place in this review where the routes
**do** disagree at alpha: paired t p = {LEA['p_t']:.4f} and permutation p = {LEA['p_perm']:.4f}
fall below 0.05 while Wilcoxon p = {LEA['p_rank']:.4f} and the sign test p = {LEA['p_sign']:.4f} do
not. On the first pass, when this contrast was wrongly taken to carry C5, that disagreement drove
an `assumption-dependent` verdict. **The paper does not report this contrast**, so the
disagreement is a property of a contrast only this review computed, and it bears on no claim in the
ledger; it is reported because the route table above is reported in full, including the rows that
disagree, and because it is the reason the per-gesture contrasts should not be read individually.
No route was selected after the fact in any contrast here.

**What the review could not test.** Whether the GAMM's smoothing parameters, its `k`, or the
`s(Xgaze_RIGHT, Ygaze_RIGHT)` term change the `trial_type` effect; whether the beta model's
zero/one rescaling affects the looking-time null; whether `log(TRIAL_DWELL_TIME)` is the right
exposure for the transition counts. All three require refitting the authors' model classes, which
this environment cannot do.

---

## 6. Resolution and evidence — stage 5

![F5 — null distribution](figures/F5_null_distribution.png)

*Sign-flip permutation null for C1 with the observed mean difference marked, and the bootstrap
distribution of the estimate with its 95% interval.*

![F7 — time-resolved pupil contrast](figures/F7_pupil_timecourse.png)

*The same contrast resolved over the 41 bins of the analysis window, in both samples. Left: the
authors' primary sample, bins exceeding |t| > {CP28['thr']:.2f} shaded, largest cluster
{CW0:.0f}–{CW1:.0f} ms of video time, cluster-mass permutation p = {CP28['p_cluster']:.4f}. Right:
`ALL_DATA`, where no bin reaches the threshold at all.*

**Resolution — every p-value next to the smallest its design could have produced.** No floor binds
anywhere in this review:

| contrast | n pairs | permutation p | smallest achievable p (sign-flip, 2^-n) | sign-test floor | resolution at n_perm = {N_PERM} | floor binds |
|----------|---------|---------------|------------------------------------------|-----------------|---------------------------------|-------------|
{floor_rows}
| C5 interaction | {INTER['n']} | {INTER['perm']['p_value']:.4f} | {INTER['perm']['min_achievable_p']:.2e} | {rs.min_sign_test_p(INTER['n']):.2e} | {INTER['perm']['p_resolution']:.0e} | {'yes' if INTER['perm']['floor_binds'] else 'no'} |
| C1 cluster permutation | {CP28['n']} | {CP28['p_cluster']:.4f} | {2.0 ** (-CP28['n']):.2e} | — | {CP28['p_resolution']:.0e} | no |

With {A['n']} within-dog pairs the sign-flip design admits p-values down to {A['min_p']:.2e}, so
every value reported here is one the design could have produced and none sits against a floor. The
binding limit is the **resolution** of the permutation itself: {N_PERM} permutations resolve no
finer than p = {A['p_resolution']:.0e}, and no number in this report should be read below that.
The place a floor would have mattered is the sign test, whose p = {A['p_sign']:.4f} on C1 reflects
that {A['n'] - A['n_pos']} of {A['n']} dogs moved against the claimed direction and
{A['n_pos']} with it — a count-based statement that is weaker than the mean-based one and agrees
with it.

**How the permutation was done here.** Signs were flipped **within dog**, which encodes the null
that congruent and incongruent are exchangeable within a dog while preserving between-dog
variation. Shuffling trial labels across the whole dataset would encode a different and weaker
null. How the paper's own 1000 bootstraps were resampled cannot be stated: `boot.glmmTMB` and
`boot.glmm.pred` resample at the model level, and whether the published intervals are those
parametric-bootstrap intervals could not be verified without the paper.

**Intervals.** C1, {A['n']} dogs: {A['effect']:.1f} a.u., 95% BCa bootstrap CI {A['ci'][0]:.1f} to
{A['ci'][1]:.1f}, n_boot = {N_BOOT}, bootstrap SE {A['boot']['bootstrap_se']:.1f}, bias
{A['boot']['bootstrap_bias']:.2f}. For the two nulls the interval *is* the finding: looking time
{DWL['effect']:.1f} ms, CI {DWL['ci'][0]:.0f} to {DWL['ci'][1]:.0f} ms; transitions
{TR['effect']:.3f}, CI {TR['ci'][0]:.2f} to {TR['ci'][1]:.2f}. Both intervals contain zero, and
both also contain differences large enough to matter behaviourally — roughly a
{max(abs(DWL['ci'][0]), abs(DWL['ci'][1])) / DWL['mean_cong'] * 100:.0f}% change in looking time
and a {max(abs(TR['ci'][0]), abs(TR['ci'][1])) / TR['mean_cong'] * 100:.0f}% change in transition
rate. "Did not differ significantly" is a correct description of these data; "did not differ" is
not supported by them, and no equivalence bound was pre-specified in anything this review could
see. This is the single change that would most strengthen the two null claims: state a smallest
effect of interest and test against it, rather than reading non-significance as absence.

---

## 7. Multiplicity — stage 6

![F6 — multiplicity](figures/F6_multiplicity.png)

*Assessable on this pass. The {len(REPORTED_P)} p-values reported in the paper's main text, against
the uniform expectation, with hits surviving each correction on the right. Bounds reported as
p < .001 are plotted at .001, which is the reading most favourable to the paper.*

**Counting what was run.** Two counts, which differ and should: {len(REPORTED_P)} p-values appear
in the paper's main text, while the authors' deposited code contains {M_TESTS} inferential tests.
The code is the larger number because it includes diagnostics and model-comparison steps the paper
does not print:

| analysis | inferential tests |
|----------|-------------------|
""" + "\n".join(f"| {name} | {n} |" for name, n in TEST_INVENTORY) + f"""
| **total visible in deposited code** | **{M_TESTS}** |

At alpha = {ALPHA} with m = {M_TESTS}, {ALPHA * M_TESTS:.2f} significant results are expected by
chance alone; with m = {len(REPORTED_P)} reported, {ALPHA * len(REPORTED_P):.2f}. **The words
"Bonferroni", "false discovery" and "correct for multiple" do not appear in the paper**, and no
correction is described in its analysis section. Recounting the reported p-values under two
standard corrections:

| correction | controls | n significant at {ALPHA} | threshold on raw p |
|------------|----------|--------------------------|--------------------|
{recount_paper_md}

Two things qualify that table before it is read as a criticism. First, the {len(REPORTED_P)}
reported p-values are not one family: they span three outcomes (pupil size, dwell proportion, gaze
transitions) with largely pre-specified hypotheses, and several are control predictors (session
number, trial number) or model-comparison steps rather than tests of the paper's claim. A blanket
correction across all {len(REPORTED_P)} is the most severe reading available, not the correct one.
Second, the claim-bearing results are not marginal: the pupil Wilcoxon at
p = {REPORTED['pupil_wilcox_p']:.3f} and the agent-type effect at p < .001 are the two that carry
conclusions, and the two null results are unaffected by any correction. What the review can say
plainly is that the family over which the paper's alpha is controlled is not stated, and stating
it would let a reader do this arithmetic themselves.

**Correction applied to this review's own p-values.** {len(recomputed_p)} permutation p-values were
computed here; correcting *them* as a family is a statement about this review, not about the paper:

| correction | controls | n significant at {ALPHA} | threshold on raw p |
|------------|----------|--------------------------|--------------------|
{recount_md}

Expected by chance across these {len(recomputed_p)}:
{recount['expected_false_positives_uncorrected']:.2f}. Uncorrected, {len(SURV_RAW)} are significant
({'; '.join(SURV_RAW) if SURV_RAW else 'none'}). Under Benjamini-Hochberg
{len(SURV_BH)} survive ({'; '.join(SURV_BH) if SURV_BH else 'none'}); under Bonferroni
{len(SURV_BONF)} survive ({'; '.join(SURV_BONF) if SURV_BONF else 'none'}). The leaving-response
contrast of C5 (p = {LEA['p_perm']:.4f}) is
{'among the survivors' if KEY_LEA in SURV_BH or 'C5 leaving contrast (permutation)' in SURV_BH else 'not among the survivors of either correction'},
which — independently of the route disagreement — is a second reason not to read it as establishing
context-specificity.

**Distribution of p-values.** The diagnostic now has inputs but still has no resolution, and that
limit is the point rather than a preamble to it. A just-below-threshold diagnostic is a property of
a *collection* of p-values: it describes a set, says nothing about any individual result in it, and
says nothing about why any result looks the way it does. **{len(REPORTED_P)} reported p-values is
below the roughly 20 needed for the diagnostic to distinguish anything from chance**
(`interpretable = {JUSTBELOW['interpretable']}`), so while the counts can be printed —
{int(JUSTBELOW['n_in_suspect_window'])} in the window {JUSTBELOW['suspect_window']} against
{int(JUSTBELOW['n_in_reference_window'])} in {JUSTBELOW['reference_window']}, binomial p for an
excess = {JUSTBELOW['binomial_p_excess']:.3f}, {int(JUSTBELOW['exactly_at_threshold'])} exactly at
threshold — they are not interpretable, and this review does not interpret them. The skill's own
note is carried rather than paraphrased: *"{JUSTBELOW['note']}"* No inference about analytic flexibility,
selective reporting or intent is drawn here; none could be drawn from {len(REPORTED_P)} values, and
none would be drawn about an individual result from any number of them.

---

## 8. Per-claim supplement

![C1](figures/panel_C1.png)

*C1 — pupil size change, {A['n']} dogs. Each dog's two values, the sign-flip permutation null, and
the four routes side by side. All agree at alpha.*

![C2](figures/panel_C2.png)

*C2 — mean dwell time to the agents, {DWL['n']} dogs. The interval, not the p-value, is the result.*

![C3](figures/panel_C3.png)

*C3 — gaze transitions per trial, {TR['n']} dogs.
{int((trd['transition_count'] == 0).sum())} of {len(trd)} trials contain no transition at all, so
the per-dog mean summarises a count with a floor.*

![C5](figures/panel_C5.png)

*C5 — the gesture-by-congruency interaction, {INTER['n']} dogs. Left: each dog's congruency
contrast within each gesture context. Centre and right: the permutation null and the routes for
their difference, {INTER['effect']:.1f} a.u. (95% CI {INTER['boot']['ci_low']:.1f} to
{INTER['boot']['ci_high']:.1f}, permutation p = {INTER['perm']['p_value']:.3f}). The paper runs no
such test, and does not need to: its context-relativity argument is the action-balanced design,
not a contrast between these two columns. The null here is consistent with that argument.*

C4, C6 and C7 carry no per-claim panel: C4 is a count, verified in section 4, and C6 and C7 are
`not-assessable` with no contrast to draw.

---

## 9. What was not assessed

{ac.uncheckable_assumptions_note()}

Additionally, specific to this review:

- **The supplementary material.** Tables S1–S4 and the supplementary figures were not supplied.
  They hold the GAMM and GLMM coefficient tables with their degrees of freedom and the authors'
  1000-sample bootstrap CIs, which is why {N_SC_NOTCHECK} of the {len(ALL_TESTS)} extracted test
  statistics could not be arithmetic-checked and why no reported interval appears anywhere in this
  report. Supplying them would close both gaps.
- **Reviewer error on the first pass, recorded rather than silently corrected.** With only the
  abstract, this review returned `supported-with-caveat` on C1 against an undisclosed sample split
  and `assumption-dependent` on C5 against a difference of significance. The paper discloses the
  split and makes no per-gesture claim; both verdicts were withdrawn on reading it. The arithmetic
  in both cases was correct — the reading of the prose was not, which is the failure mode of
  reviewing an abstract as though it were a paper.
- **The authors' models.** Neither the `mgcv::bam` factor-smooth GAMM, the `glmmTMB` beta model,
  nor the `glmer` Poisson model with offset was refitted — none is available here, and refitting
  them would be redoing the analysis rather than auditing it. The contrasts were re-derived at the
  subject level, and a cluster-based permutation test was added for pupil size because the window
  mean and the GAMM difference curve are different estimands. **A subject-level null on `ALL_DATA`
  does not demonstrate that the authors' GAMM is wrong;** it shows the effect is not visible in the
  dogs-as-units summary of the same data, and the discrepancy is one the authors can address in a
  sentence by stating which sample the published pupil result comes from.
- **Claims with no supporting data:** C6 (not self-centred) and C7 (ontogeny), both
  `not-assessable`. C6 would be changed by a first-person control condition; C7 by any
  developmental, rearing-history or longitudinal measure.
- **Inclusion and exclusion as written.** The eleven-dog exclusion was verified against the data
  (it is exactly the set lacking a condition) and the four excluded sessions (`Nala_2`, `Sahibu_2`,
  `Melody_1`, `Zserbo_1`) were reproduced from the code's own comments, but whether these criteria
  were fixed before the outcomes were examined is a methods question, readable only from the paper,
  and is not assessed.
- **Nothing was subsampled, filtered or approximated** beyond the authors' own documented filters.
  Every number here comes from the full deposited files. Where a dog drops out of a contrast it is
  because one of its two trial types is absent, and the n is stated on every panel and in every
  table.
- **Not computed, by choice:** residual diagnostics for the authors' models, power curves, and
  equivalence tests against a bound the paper did not state.

---

## 10. Reproduction

```bash
python lonardo_2026/review_lonardo_2026.py
```

| | |
|--|--|
| Seed | `random_state = {SEED}` — sign-flip permutations, bootstraps and the cluster permutation |
| Permutations / resamples | `n_perm = {N_PERM}`, `n_boot = {N_BOOT}`, cluster permutations = {N_PERM} |
| Cluster-forming threshold | two-sided t at alpha = {ALPHA}: {CP28['thr']:.3f} ({CP28['n']} dogs), {CPALL['thr']:.3f} ({CPALL['n']} dogs) |
| Bootstrap method | BCa percentile for paired differences (`resampling.bootstrap_ci`) |
| Package versions | numpy {np.__version__}, scipy {scipy.__version__}, pandas {pd.__version__}, pyarrow {pyarrow.__version__} |
| Data checksums | section 1 (sha256, first 16 hex) |
| Skill version | statistical-review v1 |
| Pipeline check | the authors' own commented values for mean dwell time to the reacting ({REPRO_AGENT['reacting']:.3f} ms) and communicating ({REPRO_AGENT['communicating']:.3f} ms) agent were reproduced to three decimals from the raw IA reports, confirming this review's dwell-time pipeline matches theirs |

---

*Findings in this report are arithmetic statements about recomputed values and about the inputs
named above. No reported value was available to compare them against, and where they differ from
the paper the cause may be an input this review could not see; the inputs used are named so that
any discrepancy can be resolved. Nothing here attributes intent.*
"""

with open(os.path.join(OUT, "REVIEW_REPORT.md"), "w", encoding="utf-8") as fh:
    fh.write(report)

summary = dict(verdict_counts={k: int(v) for k, v in counts.items()},
               figures=sorted(os.path.basename(p) for p in figs.values()),
               n_claims=int(len(ledger)), nbook=NBOOK,
               cluster_p_28=CP28["p_cluster"], cluster_p_all=CPALL["p_cluster"],
               cluster_window=[CW0, CW1],
               key_effects={k: dict(n=v["n"], effect=v["effect"], ci=v["ci"], p_perm=v["p_perm"])
                            for k, v in ALL_CONTRASTS},
               interaction=dict(n=INTER["n"], effect=INTER["effect"],
                                ci=[INTER["boot"]["ci_low"], INTER["boot"]["ci_high"]],
                                p_perm=float(INTER["perm"]["p_value"])),
               recount={k: (int(v) if isinstance(v, (int, np.integer)) else v)
                        for k, v in recount.items() if not isinstance(v, (pd.DataFrame, np.ndarray))},
               repro_agent=REPRO_AGENT, m_tests=M_TESTS)
with open(os.path.join(OUT, "review_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2, default=str)

print(json.dumps({k: summary[k] for k in
                  ("verdict_counts", "figures", "cluster_p_28", "cluster_p_all", "m_tests")},
                 indent=2, default=str))
