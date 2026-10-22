"""
Statistical review of:

  Wang A. & Horowitz A. (2026) Dogs (Canis familiaris) distinguish conspecific
  emotional chemosignals. Sci Rep 16:11176.  doi:10.1038/s41598-026-41426-1

Inputs: the paper text plus the authors' Dryad deposit (10.5061/dryad.bzkh189pt):
Part_1.xlsx, Part_2.xlsx, statistical_analyses.R, README.md.

Method: the six-stage `statistical-review` protocol; figures on `data-visualization`.

Two stages, in order:
  1. Rscript handoff/recompute_authors_tests.R   -- reproduces the authors' own R
     pipeline (wilcox.test / kruskal.test / aov / shapiro.test / glm) on their data,
     plus the same comparisons at the dog level and with a count family. Writes
     handoff/r_results/*.csv. R is used for this step because the paper's statistics
     are R outputs: V, the tie-corrected H, and R's exact-vs-normal-approximation
     switch are properties of those routines.
  2. python review_wang_horowitz_2026.py          -- this script: the three routes,
     design floors, bootstrap intervals, multiplicity, the figure spine, the claim
     ledger and the report.

Seeds fixed (SEED = 0; n_perm = 10000, n_boot = 10000).
"""

import os, re, sys, json, math, hashlib, itertools
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED, N_PERM, N_BOOT, ALPHA = 0, 10000, 10000, 0.05
np.random.seed(SEED)
SLUG = "wang_horowitz_2026"
sys.path.insert(0, "skills/statistical-review")
from utilities import claim_extraction as ce
from utilities import recompute as rc
from utilities import resampling as rs
from utilities import multiplicity as mp
from utilities import assumption_checks as ac
from utilities import review_figures as rf

DRY = os.path.join(SLUG, "inputs", "dryad")
RRES = os.path.join("handoff", "r_results")
FIGDIR = os.path.join(SLUG, "figures")
os.makedirs(FIGDIR, exist_ok=True)
OKABE = {"blue": "#0072B2", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
         "vermillion": "#D55E00", "purple": "#CC79A7", "grey": "#999999", "black": "#000000"}
PTS = [f"{i}s" for i in range(5, 65, 5)]


def _style():
    try:
        apply_figure_style(sizes=(9, 8, 7))                      # noqa: F821
    except NameError:
        plt.rcParams.update({"figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
                             "font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                             "legend.fontsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                             "axes.spines.top": False, "axes.spines.right": False,
                             "legend.frameon": False})


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


# ===========================================================================
# STAGE 0 — what was supplied
# ===========================================================================
PAPER = os.path.join(SLUG, "inputs", "paper_fulltext.txt")
raw_text = open(PAPER, encoding="utf-8", errors="replace").read()
INPUTS = {"paper": dict(path=os.path.basename(PAPER), bytes=os.path.getsize(PAPER),
                        sha256=sha(PAPER),
                        note="PMC-extracted plain text, not the publisher PDF; figure images absent"),
          "Part_1.xlsx": dict(bytes=os.path.getsize(f"{DRY}/Part_1.xlsx"), sha256=sha(f"{DRY}/Part_1.xlsx")),
          "Part_2.xlsx": dict(bytes=os.path.getsize(f"{DRY}/Part_2.xlsx"), sha256=sha(f"{DRY}/Part_2.xlsx")),
          "statistical_analyses.R": dict(bytes=os.path.getsize(f"{DRY}/statistical_analyses.R"),
                                         sha256=sha(f"{DRY}/statistical_analyses.R")),
          "README.md": dict(bytes=os.path.getsize(f"{DRY}/README.md"), sha256=sha(f"{DRY}/README.md"))}

clean = raw_text.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
clean = ce.normalize_text(clean)
clean = re.sub(r"\b([A-Za-z])\s+\((\d)\)", r"\1(\2)", clean)
clean = re.sub(r"\(\s+", "(", re.sub(r"\s+\)", ")", clean))
inventory = ce.build_claim_inventory([clean])
auto_p = len(inventory["p_values"])

RSCRIPT = open(f"{DRY}/statistical_analyses.R", encoding="utf-8").read()
SCRIPT_CALLS = {name: len(re.findall(pat, RSCRIPT))
                for name, pat in [("wilcox.test", r"wilcox\.test\("), ("lm", r"\blm\("),
                                  ("shapiro.test", r"shapiro\.test\("), ("kruskal.test", r"kruskal\.test\("),
                                  ("aov", r"\baov\("), ("dunn.test", r"dunn\.test\("),
                                  ("glm summary", r"summary\(glm\(")]}
TYPO = "Frequenc_MFy" in RSCRIPT

# ===========================================================================
# STAGE 1-2 — the data, and what the unit of replication actually is
# ===========================================================================
x1 = pd.ExcelFile(f"{DRY}/Part_1.xlsx", engine="openpyxl")
x2 = pd.ExcelFile(f"{DRY}/Part_2.xlsx", engine="openpyxl")
num = lambda s: pd.to_numeric(s, errors="coerce")
allsheet = x1.parse("all").apply(num)
pairings = x1.parse("emotionalpairings")
orders = x1.parse("specificorders")

PAIR_COLS = {
    "JB": dict(h1=["Jhabit1_vB", "Bhabit1_vJ"], h3=["Jhabit3_vB", "Bhabit3_vJ"],
               h4=["Jhabit4_vB", "Bhabit4_vJ"], nov=["Bnovel_vJ", "Jnovel_vB"], hab=("J", "B")),
    "SB": dict(h1=["Shabit1_vB", "Bhabit_vS"], h3=["Shabit3_vB", "Bhabit3_vS"],
               h4=["Shabit4_vB", "Bhabit4_vS"], nov=["Bnovel_vS", "Snovel_vB"], hab=("S", "B")),
    "SJ": dict(h1=["Shabit1_vJ", "Jhabit1_vS"], h3=["Shabit3_vJ", "Jhabit3_vS"],
               h4=["Shabit4_vJ", "Jhabit4_vS"], nov=["Jnovel_vS", "Snovel_vJ"], hab=("S", "J"))}
rows = []
for pair, cols in PAIR_COLS.items():
    for _, r in orders.iterrows():
        v, which = {}, None
        for key in ("h1", "h3", "h4", "nov"):
            ca, cb = cols[key]
            a, b = num(pd.Series([r[ca]])).iloc[0], num(pd.Series([r[cb]])).iloc[0]
            v[key] = a if pd.notna(a) else b
            if key == "h3":
                which = cols["hab"][0] if pd.notna(a) else cols["hab"][1]
        rows.append(dict(Subject=r["Subject"], Season=r["Season"], Sex=r["Sex"], pairing=pair,
                         habituated=which, habit1=v["h1"], habit3=v["h3"], habit4=v["h4"], novel=v["nov"]))
long = pd.DataFrame(rows)
long = long[long[["habit1", "habit3", "habit4", "novel"]].notna().any(axis=1)].reset_index(drop=True)
# the reconstruction must reproduce the authors' pooled sheet exactly, or the dog
# labels attached to those 122 rows are not theirs
RECON_OK = {c: sorted(long[c].dropna().round(6)) == sorted(allsheet[c].dropna().round(6))
            for c in ["habit1", "habit3", "habit4", "novel"]}
assert all(RECON_OK.values()), RECON_OK
assert len(long) == len(allsheet) == 122
dogmeans = long.groupby("Subject")[["habit1", "habit3", "habit4", "novel"]].mean().reset_index()
TRIALS_PER_DOG = long.groupby("Subject").size().value_counts().to_dict()
N_DOGS_P1 = long["Subject"].nunique()

long["discrim_diff"] = long["novel"] - long["habit4"]
long["dishab_diff"] = long["novel"] - long["habit3"]
INDEP = {k: ac.design_independence(long.dropna(subset=[k]), "Subject", "pairing", k)
         for k in ("discrim_diff", "dishab_diff", "novel")}

# --- the authors' own R output ---------------------------------------------
RW = pd.read_csv(f"{RRES}/part1_wilcoxon.csv").set_index("test_id")
RROUTES = pd.read_csv(f"{RRES}/part1_routes.csv")
RLM = pd.read_csv(f"{RRES}/part1_sidebias_lm.csv")
ROMNI = pd.read_csv(f"{RRES}/part2_omnibus.csv").set_index("sheet")
RSHAP = pd.read_csv(f"{RRES}/part2_shapiro.csv").set_index("sheet")
RDOG = pd.read_csv(f"{RRES}/part2_doglevel.csv").set_index("sheet")
RDOGQ = pd.read_csv(f"{RRES}/part2_doglevel_quasipois.csv").set_index("sheet")
RDUNN = pd.read_csv(f"{RRES}/part2_dunn.csv")
RGLM = pd.read_csv(f"{RRES}/part2_glm_gaussian.csv")
RSIMP = pd.read_csv(f"{RRES}/part2_glm_simple_effects.csv")
RSEXDOG = pd.read_csv(f"{RRES}/part2_sex_doglevel.csv")
RPW = pd.read_csv(f"{RRES}/part2_pairwise_doglevel.csv")
RGLMQ = pd.read_csv(f"{RRES}/part2_glm_dogcount.csv")

# --- Part 1: paired vectors behind each reported comparison ----------------
def pair_vec(df, ca, cb):
    a, b = num(df[ca]).values, num(df[cb]).values
    m = np.isfinite(a) & np.isfinite(b)
    return a[m], b[m]


P1_SPEC = {  # claim_id: (frame, habituated column, novel/compared column, scope, description)
    "W1":  (allsheet, "habit1", "habit3", "pooled trial", "1st vs 3rd habituation presentation"),
    "W2":  (allsheet, "habit4", "novel",  "pooled trial", "novel vs 4th habituated (discrimination)"),
    "W3":  (allsheet, "habit3", "novel",  "pooled trial", "novel vs 3rd habituated (dishabituation)"),
    "W4":  (allsheet, "habit3", "habit4", "pooled trial", "3rd vs 4th habituated (format control)"),
    "W5":  (pairings, "JBhabit4", "JBnovel", "dog (pairing)", "joy-baseline: novel vs 4th"),
    "W6":  (pairings, "JBhabit3", "JBnovel", "dog (pairing)", "joy-baseline: novel vs 3rd"),
    "W11": (pairings, "SBhabit4", "SBnovel", "dog (pairing)", "stress-baseline: novel vs 4th"),
    "W12": (pairings, "SJhabit4", "SJnovel", "dog (pairing)", "stress-joy: novel vs 4th"),
    "W14": (pairings, "SJhabit3", "SJnovel", "dog (pairing)", "stress-joy: novel vs 3rd"),
    "W7":  (orders, "Jhabit4_vB", "Bnovel_vJ", "dog (order)", "joy habituated / baseline novel: novel vs 4th"),
    "W8":  (orders, "Bhabit4_vJ", "Jnovel_vB", "dog (order)", "baseline habituated / joy novel: novel vs 4th"),
    "W9":  (orders, "Jhabit3_vB", "Bnovel_vJ", "dog (order)", "joy habituated / baseline novel: novel vs 3rd"),
    "W10": (orders, "Bhabit3_vJ", "Jnovel_vB", "dog (order)", "baseline habituated / joy novel: novel vs 3rd"),
    "W13": (orders, "Shabit4_vJ", "Jnovel_vS", "dog (order)", "stress habituated / joy novel: novel vs 4th"),
    "W15": (orders, "Shabit3_vJ", "Jnovel_vS", "dog (order)", "stress habituated / joy novel: novel vs 3rd"),
    "W16": (orders, "Jhabit3_vS", "Snovel_vJ", "dog (order)", "joy habituated / stress novel: novel vs 3rd"),
}
P1_REPORTED = {"W1": ("<", 0.001, 7108.0), "W2": ("=", 0.008, 2704.0), "W3": ("<", 0.001, 2061.5),
               "W4": ("=", 0.087, 2977.0), "W5": ("=", 0.013, 226.0), "W6": ("<", 0.001, 133.0),
               "W7": ("=", 0.004, 26.0), "W8": ("=", 0.470, 94.0), "W9": ("=", 0.003, 25.0),
               "W10": ("=", 0.042, 44.0), "W11": ("=", 0.820, 412.0), "W12": ("=", 0.095, 301.0),
               "W13": ("=", 0.001, 23.0), "W14": ("<", 0.001, 127.0), "W15": ("=", 0.003, 29.0),
               "W16": ("=", 0.045, 34.0)}

p1 = {}
for cid, (frame, ca, cb, scope, desc) in P1_SPEC.items():
    a, b = pair_vec(frame, ca, cb)
    d = b - a                                   # novel (or later presentation) minus habituated
    nz = d[d != 0]
    routes = rs.three_routes(a, b, paired=True, n_perm=N_PERM, random_state=SEED)
    # permute (novel - habituated) so the statistic, the bootstrap and the paper's
    # wording all point the same way
    perm = rs.permutation_test(b, a, paired=True, n_perm=N_PERM, random_state=SEED)
    boot = rs.bootstrap_ci(d, np.mean, n_boot=N_BOOT, random_state=SEED)
    sgn = rs.sign_test(d)
    pos = float(np.sum(stats.rankdata(np.abs(nz))[nz > 0].sum())) if len(nz) else np.nan
    neg = float(np.sum(stats.rankdata(np.abs(nz))[nz < 0].sum())) if len(nz) else np.nan
    p1[cid] = dict(scope=scope, desc=desc, n_pairs=int(len(d)), n_nonzero=int(len(nz)),
                   mean_diff=float(np.mean(d)), median_diff=float(np.median(d)),
                   rank_biserial=float((pos - neg) / (pos + neg)) if len(nz) else np.nan,
                   boot_lo=boot["ci_low"], boot_hi=boot["ci_high"],
                   routes=routes, perm=perm, sign_p=sgn["p_value"],
                   floor=rs.min_achievable_p(len(nz), paired=True),
                   sign_floor=rs.min_sign_test_p(len(nz)),
                   r_V=float(RW.loc[cid, "statistic"]), r_p=float(RW.loc[cid, "p_value"]),
                   r_exact=bool(RW.loc[cid, "exact"]),
                   groups={ca: a, cb: b})
# the four pooled comparisons repeated at the dog level
for cid, (ca, cb) in {"W1": ("habit1", "habit3"), "W2": ("habit4", "novel"),
                      "W3": ("habit3", "novel"), "W4": ("habit3", "habit4")}.items():
    a, b = pair_vec(dogmeans, ca, cb)
    d = b - a
    nz = d[d != 0]
    p1[cid]["dog"] = dict(n=int(len(d)), V=float(RW.loc[cid + "_dog", "statistic"]),
                          p=float(RW.loc[cid + "_dog", "p_value"]),
                          routes=rs.three_routes(a, b, paired=True, n_perm=N_PERM, random_state=SEED),
                          perm=rs.permutation_test(b, a, paired=True, n_perm=N_PERM, random_state=SEED),
                          boot=rs.bootstrap_ci(d, np.mean, n_boot=N_BOOT, random_state=SEED),
                          mean_diff=float(np.mean(d)), floor=rs.min_achievable_p(len(nz), paired=True),
                          groups={ca: a, cb: b})

# ===========================================================================
# Part 2 — three nested units live in every sheet
# ===========================================================================
SHEETS8 = ["stressbehav", "joybehav", "owner", "box", "stranger", "upright", "sit", "prone"]
OUTCOME = {"stressbehav": "stress behaviours", "joybehav": "joy behaviours",
           "owner": "proximity to owner", "box": "proximity to odour box",
           "stranger": "proximity to stranger", "upright": "upright posture",
           "sit": "sitting posture", "prone": "prone posture"}
COND = {"C": "blank", "J": "joy", "S": "stress"}
P2 = {}
for sh in SHEETS8:
    d = x2.parse(sh)
    rawd = d.loc[d["Subject"].notna() & d["Odor"].notna(), ["Subject", "Sex", "Odor"] + PTS].copy()
    rawd["n_points"] = rawd[PTS].notna().sum(axis=1)
    rawd["dog_total"] = rawd[PTS].sum(axis=1, skipna=True)
    rawd["dog_mean"] = rawd["dog_total"] / rawd["n_points"]
    proc = d.loc[d["Frequency"].notna(), ["Condition", "Frequency"]].copy()
    mf = (d.loc[d["Frequency_MF"].notna(), ["Condition_MF", "Frequency_MF", "Sex_MF"]].copy()
          if "Frequency_MF" in d.columns else None)
    # the processed column is literally the per-sample-point mean over dogs
    recomputed_pt = np.concatenate([rawd.loc[rawd.Odor == c, PTS].mean(axis=0, skipna=True).values
                                    for c in ["C", "J", "S"]])
    P2[sh] = dict(raw=rawd, proc=proc, mf=mf,
                  proc_matches_mean=bool(np.nanmax(np.abs(proc["Frequency"].values - recomputed_pt)) < 1e-12),
                  dogs_per_point={c: rawd.loc[rawd.Odor == c, PTS].notna().sum().tolist()
                                  for c in ["C", "J", "S"]},
                  raw_integer=bool(np.all(np.mod(rawd[PTS].stack().values, 1) == 0)),
                  raw_total=float(rawd[PTS].stack().sum()), n_cells=int(rawd[PTS].notna().sum().sum()),
                  n_dogs={c: int((rawd.Odor == c).sum()) for c in ["C", "J", "S"]})
POSTURE_DATA = {s: P2[s]["raw_total"] for s in ["upright", "sit", "prone"]}
POSTURE_PAPER = {"upright": 169.0, "sit": 20.0, "prone": 94.0}
POSTURE_CELLS = P2["upright"]["n_cells"]
N_DOGS_P2 = sum(P2["stressbehav"]["n_dogs"].values())
SEX_CELLS = {c: dict(P2["stressbehav"]["raw"].loc[P2["stressbehav"]["raw"].Odor == c, "Sex"]
                     .value_counts()) for c in ["C", "J", "S"]}

# --- the 11 reported Dunn comparisons, at the dog level --------------------
DUNN_SPEC = [("D1", "stressbehav", ("C", "J"), "<", 0.001), ("D2", "stressbehav", ("C", "S"), "<", 0.001),
             ("D3", "joybehav", ("C", "J"), "<", 0.001), ("D4", "joybehav", ("C", "S"), "<", 0.001),
             ("D5", "owner", ("C", "S"), "<", 0.001), ("D6", "owner", ("J", "S"), "<", 0.001),
             ("D7", "stranger", ("C", "S"), "=", 0.002), ("D8", "upright", ("C", "J"), "<", 0.001),
             ("D9", "upright", ("C", "S"), "=", 0.004), ("D10", "prone", ("C", "J"), "<", 0.001),
             ("D11", "prone", ("C", "S"), "=", 0.005)]
dunn = {}
for cid, sh, pair, rel, prep in DUNN_SPEC:
    rawd = P2[sh]["raw"]
    a = rawd.loc[rawd.Odor == pair[0], "dog_mean"].values
    b = rawd.loc[rawd.Odor == pair[1], "dog_mean"].values
    dr = RDUNN[(RDUNN.sheet == sh) & (RDUNN.comparison == f"{pair[0]} - {pair[1]}")].iloc[0]
    pwr = RPW[(RPW.sheet == sh) & (RPW.comparison == f"{pair[0]} - {pair[1]}")].iloc[0]
    dunn[cid] = dict(sheet=sh, pair=pair, p_relation=rel, p_reported=prep,
                     z=float(dr["z"]), p_half=float(dr["p_half"]), p_two=float(dr["p_two_sided"]),
                     p_half_holm=float(dr["p_half_holm"]), p_two_holm=float(dr["p_two_holm"]),
                     dog_W=float(pwr["W"]), dog_wilcox_p=float(pwr["wilcox_p"]),
                     dog_welch_p=float(pwr["welch_p"]),
                     routes=rs.three_routes(a, b, n_perm=N_PERM, random_state=SEED),
                     perm=rs.permutation_test(a, b, n_perm=N_PERM, random_state=SEED),
                     boot=rs.bootstrap_ci(b - np.mean(a) if len(a) == len(b) else b, np.mean,
                                          n_boot=N_BOOT, random_state=SEED),
                     mean_diff=float(np.mean(b) - np.mean(a)),
                     floor=rs.min_achievable_p(len(a), len(b)), n1=len(a), n2=len(b),
                     groups={COND[pair[0]]: a, COND[pair[1]]: b})

# --- the 8 omnibus tests ---------------------------------------------------
OMNI_SPEC = {"K1": ("stressbehav", "<", 0.001, 20.96), "K2": ("joybehav", "<", 0.001, 23.98),
             "K3": ("owner", "<", 0.001, 23.11), "K4": ("box", "=", 0.370, 1.035),
             "K5": ("stranger", "=", 0.016, 8.28), "K6": ("upright", "<", 0.001, 16.89),
             "K7": ("sit", "=", 0.079, 5.08), "K8": ("prone", "<", 0.001, 15.6326)}
omni = {}
for cid, (sh, rel, prep, stat_rep) in OMNI_SPEC.items():
    r = ROMNI.loc[sh]
    omni[cid] = dict(sheet=sh, p_relation=rel, p_reported=prep, stat_reported=stat_rep,
                     stat_name=r["statistic_name"], stat_recomputed=float(r["statistic"]),
                     method=r["method"], p_recomputed=float(r["p_value"]), n_proc=int(r["n"]),
                     shapiro_W=float(RSHAP.loc[sh, "W"]), shapiro_p=float(RSHAP.loc[sh, "p_value"]),
                     dog_kw_H=float(RDOG.loc[sh, "kw_H"]), dog_kw_p=float(RDOG.loc[sh, "kw_p"]),
                     dog_aov_F=float(RDOG.loc[sh, "aov_F"]), dog_aov_p=float(RDOG.loc[sh, "aov_p"]),
                     dog_qp_F=float(RDOGQ.loc[sh, "quasipois_F"]), dog_qp_p=float(RDOGQ.loc[sh, "quasipois_p"]),
                     dispersion=float(RDOGQ.loc[sh, "dispersion"]))

# --- the 4 reported sex x condition coefficients ---------------------------
GLM_SPEC = {"G1": ("stressbehav", "S", "Condition_MFS:Sex_MFM", -0.281, 0.100, -2.796, "=", 0.007, -39.927,
                   "stress behaviours: F > M under stress odour"),
            "G2": ("joybehav", "S", "Condition_MFS:Sex_MFM", 0.943, 0.181, 5.201, "<", 0.001, 45.290,
                   "joy behaviours: F < M under stress odour"),
            "G3": ("owner", "S", "Condition_MFS:Sex_MFM", -0.300, 0.099, -3.044, "=", 0.003, -42.484,
                   "proximity to owner: F > M under stress odour"),
            "G4": ("stranger", "J", "Condition_MFJ:Sex_MFM", 0.339, 0.122, 2.780, "=", 0.007, -11.885,
                   "proximity to stranger: M > F under joy odour")}
glm = {}
for cid, (sh, cond, term, b_rep, se_rep, t_rep, rel, p_rep, aic_rep, desc) in GLM_SPEC.items():
    g = RGLM[(RGLM.sheet == sh) & (RGLM.term == term)].iloc[0]
    s = RSIMP[(RSIMP.sheet == sh) & (RSIMP.condition == cond)].iloc[0]
    dl = RSEXDOG[(RSEXDOG.sheet == sh) & (RSEXDOG.condition == cond)].iloc[0]
    q = RGLMQ[(RGLMQ.sheet == sh) & (RGLMQ.term == f"Odor{cond}:SexMale")]
    rawd = P2[sh]["raw"]
    sub = rawd[rawd.Odor == cond]
    a = sub.loc[sub.Sex == "Female", "dog_mean"].values
    bvals = sub.loc[sub.Sex == "Male", "dog_mean"].values
    glm[cid] = dict(sheet=sh, condition=cond, term=term, desc=desc,
                    beta_reported=b_rep, se_reported=se_rep, t_reported=t_rep,
                    p_relation=rel, p_reported=p_rep, aic_reported=aic_rep,
                    beta=float(g["estimate"]), se=float(g["std_error"]), t=float(g["t_value"]),
                    p=float(g["p_value"]), aic=float(g["aic"]), resid_df=int(g["resid_df"]), n_mf=int(g["n"]),
                    simple_est=float(s["estimate"]), simple_t=float(s["t_value"]), simple_p=float(s["p_value"]),
                    dog_n_f=int(dl["n_female"]), dog_n_m=int(dl["n_male"]),
                    dog_mean_f=float(dl["mean_female"]), dog_mean_m=float(dl["mean_male"]),
                    dog_wilcox_p=float(dl["wilcox_p"]), dog_welch_p=float(dl["welch_p"]),
                    dog_qp_p=float(q["p_value"].iloc[0]) if len(q) else np.nan,
                    routes=rs.three_routes(a, bvals, n_perm=N_PERM, random_state=SEED),
                    perm=rs.permutation_test(a, bvals, n_perm=N_PERM, random_state=SEED),
                    floor=rs.min_achievable_p(len(a), len(bvals)),
                    groups={"female": a, "male": bvals})

# ===========================================================================
# STATCHECK — every reported statistic against its recomputation
# ===========================================================================
sc_rows = []
for cid, (rel, prep, vrep) in P1_REPORTED.items():
    e = p1[cid]
    win = 0.5 * 10 ** (-len(str(prep).split(".")[-1]))
    p_ok = (e["r_p"] < prep + win) if rel == "<" else (abs(e["r_p"] - prep) <= win)
    sc_rows.append(dict(claim_id=cid, test="Wilcoxon signed-rank (R wilcox.test, paired)",
                        statistic_name="V", statistic_reported=vrep, statistic_recomputed=e["r_V"],
                        statistic_match=bool(abs(vrep - e["r_V"]) < 1e-9),
                        p_relation=rel, p_reported=prep, p_recomputed=e["r_p"],
                        p_consistent=bool(p_ok), decision_error=bool((prep < ALPHA) != (e["r_p"] < ALPHA)),
                        n_used=e["n_pairs"], exact=e["r_exact"]))
for cid, e in omni.items():
    rel, prep = e["p_relation"], e["p_reported"]
    win = 0.5 * 10 ** (-len(str(prep).split(".")[-1]))
    p_ok = (e["p_recomputed"] < prep + win) if rel == "<" else (abs(e["p_recomputed"] - prep) <= win)
    sc_rows.append(dict(claim_id=cid, test=f"{e['method']} on 12 sample-point means per condition",
                        statistic_name=e["stat_name"], statistic_reported=e["stat_reported"],
                        statistic_recomputed=e["stat_recomputed"],
                        statistic_match=bool(abs(e["stat_reported"] - e["stat_recomputed"]) < 0.005),
                        p_relation=rel, p_reported=prep, p_recomputed=e["p_recomputed"],
                        p_consistent=bool(p_ok),
                        decision_error=bool((prep < ALPHA) != (e["p_recomputed"] < ALPHA)),
                        n_used=e["n_proc"], exact=False))
for cid, e in dunn.items():
    rel, prep = e["p_relation"], e["p_reported"]
    win = 0.5 * 10 ** (-len(str(prep).split(".")[-1]))
    p_ok = (e["p_half"] < prep + win) if rel == "<" else (abs(e["p_half"] - prep) <= win)
    sc_rows.append(dict(claim_id=cid, test="Dunn post-hoc (dunn.test default: halved, unadjusted)",
                        statistic_name="z", statistic_reported=np.nan, statistic_recomputed=e["z"],
                        statistic_match=pd.NA, p_relation=rel, p_reported=prep,
                        p_recomputed=e["p_half"], p_consistent=bool(p_ok),
                        decision_error=bool((prep < ALPHA) != (e["p_half"] < ALPHA)),
                        n_used=e["n1"] + e["n2"], exact=False))
for cid, e in glm.items():
    rel, prep = e["p_relation"], e["p_reported"]
    win = 0.5 * 10 ** (-len(str(prep).split(".")[-1]))
    p_ok = (e["p"] < prep + win) if rel == "<" else (abs(e["p"] - prep) <= win)
    sc_rows.append(dict(claim_id=cid, test=f"Gaussian GLM coefficient {e['term']} (72 sex x condition x point means)",
                        statistic_name="t", statistic_reported=e["t_reported"], statistic_recomputed=e["t"],
                        statistic_match=bool(abs(abs(e["t_reported"]) - abs(e["t"])) < 0.001),
                        p_relation=rel, p_reported=prep, p_recomputed=e["p"], p_consistent=bool(p_ok),
                        decision_error=bool((prep < ALPHA) != (e["p"] < ALPHA)),
                        n_used=e["n_mf"], exact=False))
statcheck = pd.DataFrame(sc_rows)
SC = dict(n=len(statcheck), stat_mismatch=int((statcheck["statistic_match"] == False).sum()),  # noqa: E712
          p_inconsistent=int((statcheck["p_consistent"] == False).sum()),                      # noqa: E712
          decision_errors=int((statcheck["decision_error"] == True).sum()))                    # noqa: E712

# ===========================================================================
# MULTIPLICITY — the full inventory the script makes visible
# ===========================================================================
INVENTORY = [("Wilcoxon signed-rank (Part 1)", SCRIPT_CALLS["wilcox.test"], 16,
              "3 pooled + 6 by pairing + 12 by presentation order + 1 format control"),
             ("Side-bias / season regression coefficients", SCRIPT_CALLS["lm"] * 3, 0,
              "6 lm models x (side + 2 season contrasts); reported only as p > 0.05"),
             ("Shapiro-Wilk route selection", SCRIPT_CALLS["shapiro.test"], 0,
              "one per Part 2 outcome; decides aov vs kruskal.test"),
             ("Part 2 omnibus", SCRIPT_CALLS["kruskal.test"] + SCRIPT_CALLS["aov"], 8,
              "7 kruskal.test + 1 aov"),
             ("Dunn pairwise comparisons", SCRIPT_CALLS["dunn.test"] * 3, 11,
              "6 dunn.test families x 3 pairwise comparisons"),
             ("GLM coefficients", SCRIPT_CALLS["glm summary"] * 5, 4,
              "5 models x (2 condition + 1 sex + 2 interaction) coefficients")]
M_RUN = sum(n for _, n, _, _ in INVENTORY)
M_PAPER = sum(n for _, _, n, _ in INVENTORY)
recomputed_p = ([p1[c]["r_p"] for c in P1_REPORTED] +
                [RW.loc[t, "p_value"] for t in ["W11b", "O5", "O6", "O7", "O8", "O11"]] +
                RLM[RLM.term != "(Intercept)"]["p_value"].tolist() + RSHAP["p_value"].tolist() +
                ROMNI["p_value"].tolist() + RDUNN["p_half"].tolist() +
                RGLM[RGLM.term != "(Intercept)"]["p_value"].tolist())
recomputed_p = [float(v) for v in recomputed_p if np.isfinite(v)]
recount = mp.discovery_recount(recomputed_p, alpha=ALPHA)
jb = mp.just_below_rate(recomputed_p)
reported_p = ([P1_REPORTED[c][1] for c in P1_REPORTED] + [omni[c]["p_reported"] for c in omni] +
              [dunn[c]["p_reported"] for c in dunn] + [glm[c]["p_reported"] for c in glm])
recount_reported = mp.discovery_recount(reported_p, alpha=ALPHA)
FAMILIES = {"Part 1 order-level (12 tests)": ["W7", "W8", "W9", "W10", "W13", "W15", "W16"],
            "Part 2 omnibus (8 tests)": list(omni), "Dunn (18 comparisons)": list(dunn),
            "GLM coefficients (25)": list(glm)}
bh_family = {}
for fam, ids in FAMILIES.items():
    if fam.startswith("Part 1 order"):
        ps = [RW.loc[t, "p_value"] for t in ["W7", "W8", "W9", "W10", "O5", "O6", "O7", "O8",
                                             "W13", "W15", "O11", "W16"]]
        lab = ["W7", "W8", "W9", "W10", "O5", "O6", "O7", "O8", "W13", "W15", "O11", "W16"]
    elif fam.startswith("Part 2 omnibus"):
        ps, lab = [omni[c]["p_recomputed"] for c in omni], list(omni)
    elif fam.startswith("Dunn"):
        ps = RDUNN["p_half"].tolist()
        lab = [f"{r.sheet}:{r.comparison}" for r in RDUNN.itertuples()]
    else:
        sub = RGLM[RGLM.term != "(Intercept)"]
        ps, lab = sub["p_value"].tolist(), [f"{r.sheet}:{r.term}" for r in sub.itertuples()]
    tab = mp.benjamini_hochberg(ps, alpha=ALPHA)
    tab["label"] = lab
    bh_family[fam] = tab
BH_LOOKUP = {}
for fam, tab in bh_family.items():
    for _, r in tab.iterrows():
        BH_LOOKUP[str(r["label"])] = bool(r["reject"])

# ===========================================================================
# THE CLAIM LEDGER
# ---------------------------------------------------------------------------
# One row per numerically reported claim, one verdict each. The cascade used:
#   not-reproducible        the test at the unit that was randomised, or the
#                           contrast that matches the claim's wording, falls on
#                           the other side of alpha from the reported p
#   assumption-dependent    parametric / rank / permutation disagree at alpha
#   supported-with-caveat   the conclusion stands, but the unit, the p-value
#                           convention, a mis-stated statistic or a missing
#                           correction qualifies it
#   supported               reproduces, all routes agree, nothing qualifies it
# ===========================================================================
HEADLINE = {
    "H1": ("dogs investigate a novel odour longer than a repeated one", ["W1", "W2", "W3", "W4"]),
    "H2": ("dogs distinguish joy from baseline odour", ["W5", "W6", "W7", "W8", "W9", "W10"]),
    "H3": ("dogs distinguish joy from stress, not stress from baseline",
           ["W11", "W12", "W13", "W14", "W15", "W16"]),
    "H4": ("stress odour moves dogs toward the owner and away from the stranger",
           ["K3", "K5", "D5", "D6", "D7"]),
    "H5": ("conspecific odours shift behaviour and posture relative to blank",
           ["K1", "K2", "K4", "K6", "K7", "K8", "D1", "D2", "D3", "D4", "D8", "D9", "D10", "D11"]),
    "H6": ("male and female dogs respond differently to the odour conditions", ["G1", "G2", "G3", "G4"])}
CLAIM_TO_HEADLINE = {c: h for h, (_, ids) in HEADLINE.items() for c in ids}
rows = []


def routes_line(routes):
    return "; ".join(f"{r.route} {r.test} p = {r.p_value:.4g}" for r in routes.itertuples())


for cid, (rel, prep, vrep) in P1_REPORTED.items():
    e = p1[cid]
    sc = statcheck[statcheck.claim_id == cid].iloc[0]
    unanimous = bool(e["routes"].attrs.get("unanimous", True))
    reasons = [f"reported V = {vrep:g}, p {rel} {prep:g}; R wilcox.test on the authors' data gives "
               f"V = {e['r_V']:g}, p = {e['r_p']:.3g} over {e['n_pairs']} pairs "
               f"({e['n_nonzero']} non-zero, {'exact' if e['r_exact'] else 'normal approximation'})"]
    if not sc["statistic_match"]:
        reasons.append(f"reported statistic differs: V = {vrep:g} against V = {e['r_V']:g}")
    if not sc["p_consistent"]:
        reasons.append(f"reported p = {prep:g} is outside the rounding window of the recomputed "
                       f"p = {e['r_p']:.4g}")
    reasons.append(f"three routes: {routes_line(e['routes'])}" +
                   ("" if unanimous else " - these disagree at alpha"))
    reasons.append(f"sign-flip floor at n = {e['n_nonzero']} is {e['floor']:.2g}; sign test p = "
                   f"{e['sign_p']:.3g} (floor {e['sign_floor']:.2g})")
    verdict = None
    if "dog" in e:
        dg = e["dog"]
        reasons.insert(1, f"the test pools {e['n_pairs']} trials from {N_DOGS_P1} dogs "
                          f"(ICC of the discrimination difference = {INDEP['discrim_diff']['icc']:.3f}, "
                          f"design effect {INDEP['discrim_diff']['design_effect']:.2f}, effective n "
                          f"{INDEP['discrim_diff']['effective_n']:.1f}); repeating it on the 43 dog means "
                          f"gives V = {dg['V']:g}, p = {dg['p']:.3g}")
        if (prep < ALPHA) != (dg["p"] < ALPHA):
            verdict = "not-reproducible"
            reasons.append(f"the conclusion changes at the unit that was randomised: reported p {rel} "
                           f"{prep:g} against dog-level p = {dg['p']:.3g}")
        else:
            verdict = "supported-with-caveat"
            reasons.append("the conclusion survives at the dog level, so the pooling inflates the "
                           "evidence but does not create the result")
    if verdict is None:
        fam_ok = BH_LOOKUP.get(cid)
        if not unanimous:
            verdict = "assumption-dependent"
        elif not sc["statistic_match"] or not sc["p_consistent"]:
            verdict = "supported-with-caveat"
        elif prep < ALPHA and fam_ok is False:
            verdict = "supported-with-caveat"
            reasons.append("does not survive Benjamini-Hochberg across the 12 order-level tests")
        else:
            verdict = "supported"
        if prep < ALPHA and fam_ok is False and verdict == "supported-with-caveat" and \
                "Benjamini" not in reasons[-1]:
            reasons.append("does not survive Benjamini-Hochberg across the 12 order-level tests")
    rows.append(dict(claim_id=cid, claim_short=e["desc"][:60], headline=CLAIM_TO_HEADLINE[cid],
                     verdict=verdict, unit_tested=e["scope"],
                     paper_statistic=f"Wilcoxon signed-rank V = {vrep:g}",
                     n_reported=e["n_pairs"], n_units=N_DOGS_P1 if "dog" in e else e["n_pairs"],
                     effect_reported=np.nan, effect_recomputed=e["mean_diff"],
                     ci_low=e["boot_lo"], ci_high=e["boot_hi"],
                     p_relation=rel, p_reported=prep, p_recomputed=e["r_p"],
                     p_permutation=e["perm"]["p_value"], p_floor_design=e["floor"],
                     floor_basis=f"two-sided sign-flip floor 2/2^n at n = {e['n_nonzero']} non-zero pairs",
                     survives_bh_in_family=BH_LOOKUP.get(cid, pd.NA),
                     reasons=" | ".join(reasons)))

for cid, e in omni.items():
    sc = statcheck[statcheck.claim_id == cid].iloc[0]
    sig_rep = e["p_reported"] < ALPHA
    reasons = [f"reported {e['stat_name']} = {e['stat_reported']:g}, p {e['p_relation']} "
               f"{e['p_reported']:g}; recomputed on the authors' data: {e['stat_name']} = "
               f"{e['stat_recomputed']:.4f}, p = {e['p_recomputed']:.3g}",
               f"the test ran on {e['n_proc']} values - the 12 per-condition sample-point means, each a "
               f"mean over {P2[e['sheet']]['n_dogs']['C']}/{P2[e['sheet']]['n_dogs']['J']}/"
               f"{P2[e['sheet']]['n_dogs']['S']} dogs, so the dogs are averaged away before testing",
               f"the same comparison on the {N_DOGS_P2} dogs' own means: Kruskal-Wallis H = "
               f"{e['dog_kw_H']:.3f}, p = {e['dog_kw_p']:.3g}; one-way ANOVA p = {e['dog_aov_p']:.3g}; "
               f"quasi-Poisson on per-dog counts p = {e['dog_qp_p']:.3g} "
               f"(dispersion {e['dispersion']:.2f})",
               f"Shapiro-Wilk on the 36 pooled means W = {e['shapiro_W']:.3f}, p = {e['shapiro_p']:.3g} - "
               f"this is the value that selected {e['method']}; it describes the shape of a three-group "
               f"mixture, not of any group or of the residuals, and it is not treated as a gate here"]
    if sig_rep != (e["dog_kw_p"] < ALPHA):
        verdict = "not-reproducible"
        reasons.append(f"the conclusion changes at the unit that was randomised: reported p "
                       f"{e['p_relation']} {e['p_reported']:g} against dog-level p = {e['dog_kw_p']:.3g}")
    elif not sig_rep:
        verdict = "supported"
        reasons.append("a null result that holds at both units and under all three dog-level routes")
    else:
        verdict = "supported-with-caveat"
        reasons.append("the conclusion survives at the dog level; the reported p-value nonetheless "
                       "describes sample points rather than dogs")
    rows.append(dict(claim_id=cid, claim_short=f"{OUTCOME[e['sheet']]} ~ odour condition"[:60],
                     headline=CLAIM_TO_HEADLINE[cid], verdict=verdict,
                     unit_tested="12 sample-point means per condition",
                     paper_statistic=f"{e['stat_name']} = {e['stat_reported']:g} ({e['method']})",
                     n_reported=e["n_proc"], n_units=N_DOGS_P2,
                     effect_reported=np.nan, effect_recomputed=np.nan, ci_low=np.nan, ci_high=np.nan,
                     p_relation=e["p_relation"], p_reported=e["p_reported"],
                     p_recomputed=e["p_recomputed"], p_permutation=np.nan,
                     p_floor_design=1.0 / math.comb(24, 9) / math.comb(15, 9),
                     floor_basis="1 / distinct 9/9/6 dog-level relabellings",
                     survives_bh_in_family=BH_LOOKUP.get(cid, pd.NA),
                     reasons=" | ".join(reasons)))

for cid, e in dunn.items():
    unanimous = bool(e["routes"].attrs.get("unanimous", True))
    sig_rep = e["p_reported"] < ALPHA
    reasons = [f"reported p {e['p_relation']} {e['p_reported']:g}; dunn.test's default output is the "
               f"halved (one-tailed) unadjusted p, and that is what reproduces: z = {e['z']:.3f}, "
               f"halved p = {e['p_half']:.3g}, two-sided p = {e['p_two']:.3g}, "
               f"two-sided Holm p = {e['p_two_holm']:.3g}",
               f"the comparison is between {e['n1']} and {e['n2']} dogs but was run on 12 vs 12 "
               f"sample-point means; on the dogs' own means Wilcoxon rank-sum p = "
               f"{e['dog_wilcox_p']:.3g}, Welch t p = {e['dog_welch_p']:.3g}",
               f"three routes at the dog level: {routes_line(e['routes'])}" +
               ("" if unanimous else " - these disagree at alpha"),
               f"dog-level permutation floor for {e['n1']} vs {e['n2']} is {e['floor']:.2g}"]
    if sig_rep != (e["dog_wilcox_p"] < ALPHA):
        verdict = "not-reproducible"
        reasons.append(f"the conclusion changes at the dog level: reported p {e['p_relation']} "
                       f"{e['p_reported']:g} against p = {e['dog_wilcox_p']:.3g} on "
                       f"{e['n1']} vs {e['n2']} dogs")
    elif not unanimous:
        verdict = "assumption-dependent"
    else:
        verdict = "supported-with-caveat"
        reasons.append("the comparison survives at the dog level; the reported value is still a halved, "
                       "unadjusted p from a 12-vs-12 test on sample-point means")
    rows.append(dict(claim_id=cid, claim_short=f"{OUTCOME[e['sheet']]}: {COND[e['pair'][0]]} vs "
                                               f"{COND[e['pair'][1]]}"[:60],
                     headline=CLAIM_TO_HEADLINE[cid], verdict=verdict,
                     unit_tested="12 vs 12 sample-point means",
                     paper_statistic="Dunn post-hoc (no statistic reported)",
                     n_reported=24, n_units=e["n1"] + e["n2"],
                     effect_reported=np.nan, effect_recomputed=e["mean_diff"],
                     ci_low=np.nan, ci_high=np.nan,
                     p_relation=e["p_relation"], p_reported=e["p_reported"], p_recomputed=e["p_half"],
                     p_permutation=e["perm"]["p_value"], p_floor_design=e["floor"],
                     floor_basis=f"two-sided label-permutation floor 2/C({e['n1'] + e['n2']},"
                                 f"{min(e['n1'], e['n2'])}) at the dog level",
                     survives_bh_in_family=BH_LOOKUP.get(f"{e['sheet']}:{e['pair'][0]} - {e['pair'][1]}", pd.NA),
                     reasons=" | ".join(reasons)))

for cid, e in glm.items():
    unanimous = bool(e["routes"].attrs.get("unanimous", True))
    cond_name = COND[e["condition"]]
    reasons = [f"reported beta = {e['beta_reported']:g}, SE = {e['se_reported']:g}, t = "
               f"{e['t_reported']:g}, p {e['p_relation']} {e['p_reported']:g}, AIC = {e['aic_reported']:g}; "
               f"refitting gives beta = {e['beta']:.4f}, SE = {e['se']:.4f}, t = {e['t']:.3f}, p = "
               f"{e['p']:.3g}, AIC = {e['aic']:.3f} on {e['n_mf']} rows "
               f"(residual df = {e['resid_df']}) - an exact reproduction",
               f"the coefficient is the {e['term']} interaction, a difference of differences against the "
               f"blank condition, while the sentence claims a sex difference within the {cond_name} "
               f"condition; that contrast is estimate {e['simple_est']:.3f}, t = {e['simple_t']:.3f}, "
               f"p = {e['simple_p']:.3g}",
               f"the {cond_name} condition contains {e['dog_n_f']} females and {e['dog_n_m']} males; on "
               f"their own means Wilcoxon rank-sum p = {e['dog_wilcox_p']:.3g}, Welch t p = "
               f"{e['dog_welch_p']:.3g}, quasi-Poisson interaction p = {e['dog_qp_p']:.3g}",
               f"three routes at the dog level: {routes_line(e['routes'])}" +
               ("" if unanimous else " - these disagree at alpha"),
               f"dog-level permutation floor for {e['dog_n_f']} vs {e['dog_n_m']} is {e['floor']:.4f}, so "
               f"the reported p {e['p_relation']} {e['p_reported']:g} is "
               f"{'below' if e['p_reported'] < e['floor'] else 'above'} the smallest value relabelling "
               f"those dogs could produce",
               "the outcome is a mean of per-dog counts over 3-5 dogs, not a count, so the missing "
               "family argument is not the issue; the unit is"]
    if (e["p_reported"] < ALPHA) != (e["dog_wilcox_p"] < ALPHA):
        verdict = "not-reproducible"
        reasons.append(f"the conclusion changes at the unit that was randomised: reported p "
                       f"{e['p_relation']} {e['p_reported']:g} against p = {e['dog_wilcox_p']:.3g} on "
                       f"{e['dog_n_f'] + e['dog_n_m']} dogs")
    elif not unanimous:
        verdict = "assumption-dependent"
    else:
        verdict = "supported-with-caveat"
    rows.append(dict(claim_id=cid, claim_short=e["desc"][:60], headline="H6", verdict=verdict,
                     unit_tested="72 sex x condition x sample-point means",
                     paper_statistic=f"Gaussian GLM beta = {e['beta_reported']:g}, t = {e['t_reported']:g}",
                     n_reported=e["n_mf"], n_units=e["dog_n_f"] + e["dog_n_m"],
                     effect_reported=e["beta_reported"], effect_recomputed=e["beta"],
                     ci_low=e["beta"] - 1.96 * e["se"], ci_high=e["beta"] + 1.96 * e["se"],
                     p_relation=e["p_relation"], p_reported=e["p_reported"], p_recomputed=e["p"],
                     p_permutation=e["perm"]["p_value"], p_floor_design=e["floor"],
                     floor_basis=f"two-sided label-permutation floor 2/C({e['dog_n_f'] + e['dog_n_m']},"
                                 f"{min(e['dog_n_f'], e['dog_n_m'])}) at the dog level",
                     survives_bh_in_family=BH_LOOKUP.get(f"{e['sheet']}:{e['term']}", pd.NA),
                     reasons=" | ".join(reasons)))

ledger = pd.DataFrame(rows)
LEDGER_COLS = ["claim_id", "claim_short", "headline", "verdict", "unit_tested", "paper_statistic",
               "n_reported", "n_units", "effect_reported", "effect_recomputed", "ci_low", "ci_high",
               "p_relation", "p_reported", "p_recomputed", "p_permutation", "p_floor_design",
               "floor_basis", "survives_bh_in_family", "reasons"]
ledger = ledger[LEDGER_COLS]
assert len(ledger) == 39 and ledger["verdict"].notna().all()
VERDICTS = ledger["verdict"].value_counts().to_dict()

ledger.to_csv(os.path.join(SLUG, "claim_ledger.csv"), index=False)
statcheck.to_csv(os.path.join(SLUG, "statcheck.csv"), index=False)

# ===========================================================================
# FIGURES — the fixed six-figure spine, now with the authors' data behind it
# ===========================================================================
_style()
FIGS = {}
VD = dict(zip(ledger.claim_id, ledger.verdict))

FIGS["F1_design_map"] = os.path.join(FIGDIR, "F1_design_map.png")
_f1 = (long.dropna(subset=["discrim_diff"])
       .rename(columns={"discrim_diff": "discrimination difference (s)", "pairing": "odour pairing"}))
rf.fig_design_map(_f1, "Subject", "odour pairing", "discrimination difference (s)",
                  independence=INDEP["discrim_diff"], path=FIGS["F1_design_map"],
                  unit_label="subject dog")

owner_raw = P2["owner"]["raw"]
owner_groups = {COND[c]: owner_raw.loc[owner_raw.Odor == c, "dog_mean"].values for c in ["C", "J", "S"]}
FIGS["F2_raw_behind_claim"] = os.path.join(FIGDIR, "F2_raw_behind_claim.png")
rf.fig_raw_behind_claim(owner_groups, "proximity to owner (per-dog mean)",
                        path=FIGS["F2_raw_behind_claim"],
                        test_name=f"24 dogs: H(2) = {omni['K3']['dog_kw_H']:.2f}, "
                                  f"p = {omni['K3']['dog_kw_p']:.2f}",
                        p_value=omni["K3"]["p_reported"],
                        reported_error="S.E.M. across 12 sample points, not across dogs")

w2 = p1["W2"]["dog"]
FIGS["F3_assumption_shape"] = os.path.join(FIGDIR, "F3_assumption_shape.png")
rf.fig_assumption_qq(owner_groups, path=FIGS["F3_assumption_shape"],
                     paired_differences=(w2["groups"]["novel"] - w2["groups"]["habit4"]),
                     value_label="proximity to owner (per-dog mean); inset: per-dog novel - habituated (s)")

FIGS["F4_claim_ledger"] = os.path.join(FIGDIR, "F4_claim_ledger.png")
f4 = rf.fig_claim_ledger(ledger, path=None,
                         effect_label="recomputed effect (seconds, score, or GLM beta — see unit column)")
f4.suptitle("F4 - Claim ledger", fontsize=12, x=0.01, ha="left", y=1.004)
f4.savefig(FIGS["F4_claim_ledger"], dpi=300, bbox_inches="tight")
plt.close(f4)

FIGS["F5_null_distribution"] = os.path.join(FIGDIR, "F5_null_distribution.png")
rf.fig_null_distribution(w2["perm"], boot=w2["boot"], path=FIGS["F5_null_distribution"],
                         statistic_label="mean per-dog difference, novel - habituated (s), 43 dogs")

FIGS["F6_multiplicity"] = os.path.join(FIGDIR, "F6_multiplicity.png")
f6 = rf.fig_multiplicity(recomputed_p, alpha=ALPHA, path=None)
f6.axes[1].set_ylim(0, max(recount["summary"]["n_significant"]) * 1.32)
f6.savefig(FIGS["F6_multiplicity"], dpi=300, bbox_inches="tight")
plt.close(f6)

PANELS = {
    "H1": dict(cid="W2", groups={"habituated, 4th presentation": w2["groups"]["habit4"],
                                 "novel odour": w2["groups"]["novel"]},
               label="investigation time (s), per-dog mean of up to 3 trials",
               perm=w2["perm"], boot=w2["boot"], routes=w2["routes"]),
    "H2": dict(cid="W5", groups=p1["W5"]["groups"], label="investigation time (s), joy-baseline pairing",
               perm=p1["W5"]["perm"], boot=None, routes=p1["W5"]["routes"]),
    "H3": dict(cid="W13", groups=p1["W13"]["groups"],
               label="investigation time (s), stress habituated / joy novel",
               perm=p1["W13"]["perm"], boot=None, routes=p1["W13"]["routes"]),
    "H4": dict(cid="D5", groups=dunn["D5"]["groups"], label="proximity to owner, per-dog mean",
               perm=dunn["D5"]["perm"], boot=None, routes=dunn["D5"]["routes"]),
    "H5": dict(cid="D3", groups=dunn["D3"]["groups"], label="joy behaviours, per-dog mean",
               perm=dunn["D3"]["perm"], boot=None, routes=dunn["D3"]["routes"]),
    "H6": dict(cid="G1", groups=glm["G1"]["groups"],
               label="stress behaviours under the stress odour, per-dog mean",
               perm=glm["G1"]["perm"], boot=None, routes=glm["G1"]["routes"]),
}
for hid, spec in PANELS.items():
    p = os.path.join(FIGDIR, f"panel_{hid}.png")
    rf.fig_claim_panel(f"{hid} / {spec['cid']} - {HEADLINE[hid][0]}", spec["groups"], spec["label"],
                       spec["perm"], boot=spec["boot"], routes=spec["routes"],
                       verdict=VD[spec["cid"]], path=p)
    FIGS[f"panel_{hid}"] = p
for k in list(FIGS):
    assert os.path.exists(FIGS[k]), k
plt.close("all")

# ===========================================================================
# A few descriptive checks the data now settle directly
# ===========================================================================
age = num(orders["Age_months"]).dropna()
AGE = dict(n=int(len(age)), mean_months=float(age.mean()), sd_months=float(age.std(ddof=1)),
           sem_months=float(age.std(ddof=1) / math.sqrt(len(age))),
           min_months=float(age.min()), max_months=float(age.max()),
           paper_pm=40.1, paper_mean_years=5.3)
SEASON = orders["Season"].value_counts().to_dict()

# ===========================================================================
# REPORT
# ===========================================================================
_ids = os.path.join(SLUG, "artifact_ids.json")
ART = json.load(open(_ids)) if os.path.exists(_ids) else {}


def img(key, caption):
    src = ("{{artifact:art_" + ART[key] + "}}") if key in ART else FIGS[key].replace("\\", "/")
    return f"![{caption}]({src})"


def md_table(df, cols=None, fmt="{:.4g}"):
    d = df[cols] if cols else df
    lines = ["| " + " | ".join(str(c) for c in d.columns) + " |",
             "|" + "|".join(["---"] * len(d.columns)) + "|"]
    for _, r in d.iterrows():
        cells = []
        for v in r:
            if isinstance(v, (bool, np.bool_)):
                cells.append("yes" if v else "no")
            elif v is None or (isinstance(v, float) and not np.isfinite(v)) or v is pd.NA:
                cells.append("—")
            elif isinstance(v, float):
                cells.append(fmt.format(v))
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


VORDER = ["supported", "supported-with-caveat", "assumption-dependent",
          "underpowered-resolution", "not-reproducible", "contradicted", "not-assessable"]
vcount = pd.DataFrame([dict(verdict=v, n=int(VERDICTS.get(v, 0))) for v in VORDER])
flip_p2 = [c for c in list(omni) + list(dunn) if ledger.set_index("claim_id").loc[c, "verdict"] == "not-reproducible"]
P = []
P.append(f"""# Statistical review — Wang & Horowitz (2026), *Dogs (Canis familiaris) distinguish conspecific emotional chemosignals*

*Sci Rep* **16**:11176 · doi:10.1038/s41598-026-41426-1 · review slug `{SLUG}`
Method: the six-stage `statistical-review` protocol; figures on `data-visualization`.
**Revision 2** — the authors' raw data arrived after the first pass. Every claim that was
`not-assessable` in revision 1 has been re-derived against their own files and re-issued.
Section 2.1 states which revision-1 inferences the data confirmed and which it overturned.

## 1. Header — what was supplied

| input | bytes | sha256 (first 16) | corresponds to |
|---|---|---|---|
| `paper_fulltext.txt` | {INPUTS['paper']['bytes']:,} | `{INPUTS['paper']['sha256'][:16]}` | the article text (PMC extraction; figure images still absent) |
| `Part_1.xlsx` | {INPUTS['Part_1.xlsx']['bytes']:,} | `{INPUTS['Part_1.xlsx']['sha256'][:16]}` | Figs 3–5: investigation times, 3 sheets |
| `Part_2.xlsx` | {INPUTS['Part_2.xlsx']['bytes']:,} | `{INPUTS['Part_2.xlsx']['sha256'][:16]}` | Fig 6 and the sex analyses: 8 behavioural sheets + 10 per-sex sheets |
| `statistical_analyses.R` | {INPUTS['statistical_analyses.R']['bytes']:,} | `{INPUTS['statistical_analyses.R']['sha256'][:16]}` | the analysis as run |
| `README.md` | {INPUTS['README.md']['bytes']:,} | `{INPUTS['README.md']['sha256'][:16]}` | variable dictionary |

**What the files contain.** `Part_1.xlsx` holds three sheets of the same trials at two different
units: `all` ({len(allsheet)} rows, one per odour-pairing trial, **no subject column**),
`emotionalpairings` ({len(pairings)} rows, one per dog) and `specificorders` ({len(orders)} rows, one
per dog). Dog identity was restored for the {len(allsheet)} pooled rows by rebuilding them from
`specificorders`; the reconstruction reproduces all four columns of `all` exactly as multisets
(`habit1`, `habit3`, `habit4`, `novel`), which is what licenses the dog-level re-analysis below.
Each sheet of `Part_2.xlsx` is three column blocks of different lengths side by side, not one table:
{N_DOGS_P2} rows of per-dog 5-s scores, {len(P2['stressbehav']['proc'])} rows of
`Condition`/`Frequency`, and {len(P2['stressbehav']['mf'])} rows of
`Condition_MF`/`Frequency_MF`/`Sex_MF`. The `Frequency` column is the mean over dogs for each
sample point — verified to machine precision for all eight outcomes — so the three blocks are three
different units of analysis, and the script picks a different one in each section.

**Coverage.** All {len(ledger)} numerically reported claims are in the ledger with a verdict, and
all were re-derived from the authors' data. Six comparisons the paper describes but never reports
numerically were recovered from the script and are reported here for the first time. Still out of
reach: the figure images (not supplied, so plotted values and significance brackets in Figs 3–6 were
not checked), and the design questions in section 9.

## 2. Summary of findings

{md_table(vcount)}

Every reported statistic reproduces. `statcheck.csv` checks {SC['n']} reported
(statistic, p) pairs against recomputation: {SC['decision_errors']} decision errors,
{SC['p_inconsistent']} p-values outside their rounding window, {SC['stat_mismatch']} statistic
mismatches. The three exceptions are all last-digit matters with no bearing on any conclusion:
W3 prints *V* = {P1_REPORTED['W3'][2]:g} against a recomputed {p1['W3']['r_V']:g}; W2 prints
*p* = {P1_REPORTED['W2'][1]:g} for {p1['W2']['r_p']:.5f}; W16 prints *p* = {P1_REPORTED['W16'][1]:g}
for {p1['W16']['r_p']:.5f}. What fails is not the arithmetic but the unit the arithmetic was done on.

**Finding 1 — Part 2's battery is a test on time bins, and most of it does not survive the dogs
(stages 2–3).** Every Kruskal–Wallis and Dunn comparison in Part 2 ran on the 12 per-condition
sample-point means, each itself a mean over {P2['owner']['n_dogs']['C']}/{P2['owner']['n_dogs']['J']}/{P2['owner']['n_dogs']['S']}
dogs. Re-running the same comparisons on the {N_DOGS_P2} dogs' own means:
{sum(1 for c in omni if omni[c]['p_reported'] < ALPHA and omni[c]['dog_kw_p'] >= ALPHA)} of
{sum(1 for c in omni if omni[c]['p_reported'] < ALPHA)} significant omnibus results and
{sum(1 for c in dunn if dunn[c]['dog_wilcox_p'] >= ALPHA)} of {len(dunn)} Dunn comparisons no longer
cross alpha, including the paper's central attachment result — proximity to owner, reported
*p* < 0.001, dog-level Kruskal–Wallis *p* = {omni['K3']['dog_kw_p']:.3f} (ANOVA
{omni['K3']['dog_aov_p']:.3f}, quasi-Poisson on per-dog counts {omni['K3']['dog_qp_p']:.3f}). Joy
behaviours is the one outcome that survives ({omni['K2']['dog_kw_p']:.3f},
{omni['K2']['dog_aov_p']:.3f}, {omni['K2']['dog_qp_p']:.3f}). Resolution: analyse the per-dog scores
that are already in the deposit, or model the 5-s bins with dog as a random effect.

**Finding 2 — the four sex claims cite an interaction but are worded as a simple effect
(stage 3).** Each reported coefficient is the `Condition × Sex` term, a difference of differences
against the blank condition, while each sentence claims a sex difference *within* one odour
condition. For two of the four the claimed contrast is not significant:
proximity to owner under stress, {glm['G3']['simple_est']:+.3f}, *p* = {glm['G3']['simple_p']:.3f}
(reported *p* = {glm['G3']['p_reported']:g}), and proximity to stranger under joy,
{glm['G4']['simple_est']:+.3f}, *p* = {glm['G4']['simple_p']:.3f} (reported
{glm['G4']['p_reported']:g}). At the dog level none of the four survives: Wilcoxon rank-sum
*p* = {', '.join(f"{glm[c]['dog_wilcox_p']:.2f}" for c in ['G1','G2','G3','G4'])} on
{glm['G1']['dog_n_f']}–{glm['G1']['dog_n_m']} dogs per cell. Resolution: report the simple effects
the sentences describe, with the dog as the unit.

**Finding 3 — Part 1 holds, but its own negative control flips (stages 2 and 5).** The pooled
Part 1 tests mix units with the rest of the section ({len(allsheet)} trials against
{N_DOGS_P1} dogs), yet the intraclass correlation of the discrimination difference is
{INDEP['discrim_diff']['icc']:.3f} (effective n {INDEP['discrim_diff']['effective_n']:.0f} of
{len(allsheet)}) — so for the discrimination comparison pooling barely inflates precision, though
for the dishabituation difference it does (ICC {INDEP['dishab_diff']['icc']:.2f}, effective n
{INDEP['dishab_diff']['effective_n']:.0f}). Habituation, discrimination and dishabituation all
survive on the 43 dog means
(*p* = {p1['W1']['dog']['p']:.2g}, {p1['W2']['dog']['p']:.3f}, {p1['W3']['dog']['p']:.2g}). The
exception is the format control: "no dishabituation due to the change in presentation format" is
*p* = {p1['W4']['r_p']:.3f} on trials but *p* = {p1['W4']['dog']['p']:.3f} on dogs, with
investigation rising from the third to the fourth presentation. That control is what separates
dishabituation to the novel odour from a response to seeing two boxes. Resolution: report the
dog-level control, or rest the dishabituation claim on the novel-vs-simultaneous comparison (W2)
alone.

## 2.1 Revision 1 inferences, checked against the data

| revision-1 inference (text only) | outcome |
|---|---|
| Pooled Part 1 tests ran on ~119–122 trial-level pairs, not 43 dogs (from *V* ≤ n(n+1)/2) | **confirmed** — sheet `all` is exactly {len(allsheet)} trial rows with no subject column |
| Part 2 tests ran on 12 sample-point means per condition (from the H ceiling of 23.00 for 9/9/6 dogs) | **confirmed** — `Frequency` is the per-sample-point mean over dogs, 12 per condition |
| GLM residual df ≈ 46–79 for G1, inconsistent with 24 dogs or 36 means | **confirmed and pinned** — residual df = {glm['G1']['resid_df']} on {glm['G1']['n_mf']} sex × condition × sample-point means |
| The four sex p-values lie below the 0.0159 dog-level permutation floor | **confirmed** — dog-level tests give *p* ≥ {min(glm[c]['dog_wilcox_p'] for c in glm):.2f} |
| "Unreported comparisons showed no statistical significance" | **confirmed** — the 6 recovered tests give *p* = {', '.join(f"{RW.loc[t,'p_value']:.2f}" for t in ['W11b','O5','O6','O7','O8','O11'])} |
| Posture counts fall 5 short of 24 × 12 = 288 | **refined** — {POSTURE_CELLS} cells were codable ({288 - POSTURE_CELLS} missing 60-s cells), and the data give upright = {POSTURE_DATA['upright']:.0f} against {POSTURE_PAPER['upright']:.0f} printed, so the gap is {288 - POSTURE_CELLS} uncoded cells plus a one-count difference in the paper |
| The age dispersion "± 40.1 months" must be an SD, not an SEM | **confirmed** — `Age_months` gives mean {AGE['mean_months']:.1f}, SD {AGE['sd_months']:.1f}, SEM {AGE['sem_months']:.1f} over n = {AGE['n']} |
| Pooling would inflate the evidence for the Part 1 claims | **overturned** — ICC ≈ {INDEP['discrim_diff']['icc']:.2f}; effective n is {INDEP['discrim_diff']['effective_n']:.0f} of {len(allsheet)}, and all three claims survive at the dog level. The unit problem that matters is in Part 2, where averaging over dogs *created* the significance |
| The missing `family` argument makes these count models Gaussian-on-counts | **overturned** — `Frequency_MF` is a mean of per-dog counts over 3–5 dogs, not a count, so the family is not the issue; the unit is |
| F2, F3, F5 and all six per-claim panels are not assessable | **overturned** — all are built from the authors' data in this revision |

## 3. Claim ledger — F4

{img('F4_claim_ledger', 'F4 claim ledger: reported vs recomputed p-values and recomputed effects')}

{md_table(ledger[['claim_id', 'claim_short', 'verdict', 'unit_tested', 'p_relation', 'p_reported', 'p_recomputed', 'p_permutation', 'p_floor_design']], fmt='{:.3g}')}

`p_recomputed` is the authors' own test re-run on their data; `p_permutation` is the label- or
sign-flip permutation at the unit that was randomised. Full reasons, effect sizes, bootstrap
intervals and the floor basis per row are in [`claim_ledger.csv`](claim_ledger.csv); the
statistic-by-statistic check is in [`statcheck.csv`](statcheck.csv).
""")

P.append(f"""## 4. Design and sampling — F1, F2

{img('F1_design_map', 'F1 design map: observations per dog, units per pairing, effective n')}

**Population.** {N_DOGS_P1} pet dogs recruited from one laboratory's owner mailing list in New York
City, screened on six owner-reported criteria including being comfortable with new people and
treat-motivated; mean age {AGE['mean_months']:.1f} months (SD {AGE['sd_months']:.1f}, range
{AGE['min_months']:.0f}–{AGE['max_months']:.0f}), collected in three periods
({', '.join(f'{k} n={v}' for k, v in SEASON.items())}). Odour stimuli came from **one** donor dog, so
every discrimination claim is a claim about discriminating that donor's three sample types. The
claims are supported for socially confident, food-motivated urban pet dogs whose owners volunteer
for cognition studies; generalising beyond that needs a sample not filtered on sociability.

**What is n, and of what.** F1 shows the Part 1 structure: {len(allsheet)} trials from
{N_DOGS_P1} dogs ({', '.join(f'{v} dogs x {k}' for k, v in sorted(TRIALS_PER_DOG.items(), reverse=True))}
trials). The ICC of the discrimination difference is {INDEP['discrim_diff']['icc']:.3f} and of the
dishabituation difference {INDEP['dishab_diff']['icc']:.3f}; of the novel-odour investigation time
itself {INDEP['novel']['icc']:.3f}. Clustering is therefore mild for the differences the claims rest
on, which is why the pooled and dog-level answers agree for W1–W3. The reporting inconsistency
stands on its own: the three pooled tests use trial rows while every by-pairing and by-order test in
the same section uses dog rows, and the paper does not say so.

Part 2 nests three units in every sheet — {N_DOGS_P2} dogs (blank {P2['owner']['n_dogs']['C']},
joy {P2['owner']['n_dogs']['J']}, stress {P2['owner']['n_dogs']['S']}), {len(P2['owner']['proc'])}
condition-by-sample-point means, and {len(P2['owner']['mf'])} condition-by-sex-by-sample-point means
— and the script tests the second for the omnibus comparisons and the third for the sex models. The
dogs never appear as the unit. A further wrinkle: the number of dogs behind each sample-point mean
is not constant. For the 60-s bin it is
{', '.join(f"{COND[c]} {P2['stressbehav']['dogs_per_point'][c][-1]}" for c in ['C', 'J', 'S'])}
instead of {P2['stressbehav']['n_dogs']['C']}/{P2['stressbehav']['n_dogs']['J']}/{P2['stressbehav']['n_dogs']['S']},
so the 12 values entering each test are means over different denominators and are not equally
precise.

{img('F2_raw_behind_claim', 'F2 raw data behind the owner-proximity claim: every dog, SD against SEM')}

**What the error bars describe.** F2 is the figure the paper does not contain. Left: every dog's
own score for the claim in H4. Right: mean ± SD beside mean ± SEM. Figs 3–6 of the paper plot
S.E.M., and Fig 6 states its unit explicitly as "means of the 12 sample points", so the published
interval describes the precision of a time-bin average, not variability between dogs — which is the
quantity the claims are about. The subject-age dispersion ("± 40.1 months") names neither quantity;
the deposit settles it as the SD ({AGE['sd_months']:.1f} months; the SEM is
{AGE['sem_months']:.1f}).

**Internal arithmetic.** One posture is scored per coded sample point. The deposit contains
{POSTURE_CELLS} coded cells ({288 - POSTURE_CELLS} of the {N_DOGS_P2} × 12 possible are missing at
60 s) and the posture counts total
upright {POSTURE_DATA['upright']:.0f} + sit {POSTURE_DATA['sit']:.0f} + prone
{POSTURE_DATA['prone']:.0f} = {sum(POSTURE_DATA.values()):.0f}, exactly the number of coded cells.
The Discussion prints upright {POSTURE_PAPER['upright']:.0f}, one fewer. The subject accounting
closes: {N_DOGS_P2} analysed in Part 2 + 2 excluded for not investigating + 3 discontinued + 14 not
offered Part 2 = {N_DOGS_P1}.

**Confounders the paper tests.** Two rooms, an owner-distraction task added after the first 14
owners, and three collection periods, entered into six side-bias regressions reported only as
*p* > 0.05. Re-running them: {int((RLM[RLM.term != '(Intercept)']['p_value'] < ALPHA).sum())} of
{int((RLM.term != '(Intercept)').sum())} coefficients across the six models reach *p* < {ALPHA},
and the largest model R² is {RLM['r_squared'].max():.3f}. A non-significant test is not a
demonstration of no effect, and with {'/'.join(str(v) for v in SEASON.values())} dogs per period
these models have little resolution; the coefficients and intervals belong in the paper.

## 5. Method and assumptions — F3

{img('F3_assumption_shape', 'F3 distribution shape: Q-Q per condition and on the paired differences')}

Shape is described here, never graded, and no test in this paper is called invalid on the strength
of a normality diagnostic.

**The route was chosen by a normality test, and that choice is invisible in the paper.**
`statistical_analyses.R` runs `shapiro.test` on all eight Part 2 outcomes and branches on the
result: `box` (*W* = {RSHAP.loc['box', 'W']:.3f}, *p* = {RSHAP.loc['box', 'p_value']:.3f}) goes to
`aov`, the other seven to `kruskal.test`. Two of those seven sit just below the threshold — owner
*p* = {RSHAP.loc['owner', 'p_value']:.3f} and stranger *p* = {RSHAP.loc['stranger', 'p_value']:.3f}
— so which route those outcomes took was decided by a diagnostic within 0.01 of its cut-off. Two
observations about the diagnostic itself, both descriptive: it is computed on the 36 values pooled
across the three conditions, so it describes a three-group mixture rather than any group or the
model residuals, and a mixture of three groups with different locations departs from normality
*because* the groups differ; and at n = 36 the test has modest power either way. The consequence
worth reporting is procedural: a reader of the paper cannot tell which comparisons were parametric.

**The Dunn p-values are halved and unadjusted.** `dunn.test` prints one-tailed values by default,
and those are the numbers in the paper: for proximity to stranger the paper reports
*p* = {dunn['D7']['p_reported']:g} against a halved value of {dunn['D7']['p_half']:.4f} and a
two-sided value of {dunn['D7']['p_two']:.4f}; for prone posture, {dunn['D11']['p_reported']:g}
against {dunn['D11']['p_half']:.4f} and {dunn['D11']['p_two']:.4f}. All eleven reported Dunn values
match the halved convention and none matches the two-sided one. No decision changes — the largest
two-sided value is {max(d['p_two'] for d in dunn.values()):.3f} — and no adjustment was applied
either, though `dunn.test` can apply one. Resolution: state the convention, or report two-sided
adjusted values.

**Interaction against simple effect.** The sex models are
`glm(Frequency_MF ~ Condition_MF * Sex_MF)` with no family argument, i.e. Gaussian least squares.
The outcome is a mean of per-dog counts over 3–5 dogs, so it is not a count and the missing family
argument is not the defect; the unit is, and the contrast is. Each reported coefficient is an
interaction term; the matching within-condition contrasts are
{'; '.join(f"{c} {glm[c]['simple_est']:+.3f} (p = {glm[c]['simple_p']:.3g})" for c in ['G1', 'G2', 'G3', 'G4'])}.
G1 and G2 agree with their sentences; G3 and G4 do not.

**The three routes.** Every claim now carries parametric, rank-based and permutation routes at the
unit that was randomised; they appear per row in the ledger and on each per-claim panel. They agree
for most claims. Where they disagree the verdict is `assumption-dependent` and the disagreement is
named:
{'; '.join(f"{r.claim_id} ({r.claim_short})" for r in ledger[ledger.verdict == 'assumption-dependent'].itertuples()) or 'none'}.
The pattern in those rows is the one worth stating: the signed-rank and t routes reach alpha while
the **sign test does not**, because the effect is carried by a large difference in a minority of
dogs rather than by most dogs moving the same way. For W5, {int(np.sum(np.sign(p1['W5']['groups']['JBnovel'] - p1['W5']['groups']['JBhabit4']) > 0))}
of {p1['W5']['n_pairs']} dogs investigated the novel odour longer (sign test
*p* = {p1['W5']['sign_p']:.3f}) against Wilcoxon *p* = {p1['W5']['r_p']:.4f}. "Dogs investigate the
novel odour longer" and "most dogs investigate the novel odour longer" are different claims, and the
data support the first more strongly than the second.

## 6. Resolution and evidence — F5

{img('F5_null_distribution', 'F5 permutation null and bootstrap interval for the pooled discrimination claim')}

{md_table(pd.DataFrame([dict(claim=c, design=ledger.set_index('claim_id').loc[c, 'floor_basis'].split(' at ')[0],
                             n=int(ledger.set_index('claim_id').loc[c, 'n_units']),
                             smallest_possible_p=float(ledger.set_index('claim_id').loc[c, 'p_floor_design']),
                             reported_p=float(ledger.set_index('claim_id').loc[c, 'p_reported']))
                        for c in ['W2', 'W13', 'D5', 'D7', 'G1']]), fmt='{:.3g}')}

No design floor in this paper exceeds α = {ALPHA}, so no claim is `underpowered-resolution`: each
design could in principle have produced a significant result. Two resolution facts still shape how
the numbers read. The four sex-by-condition p-values are smaller than the smallest value a
relabelling of the dogs in that condition could produce ({glm['G1']['floor']:.4f} for
{glm['G1']['dog_n_f']} vs {glm['G1']['dog_n_m']} dogs), which is the arithmetic signature of a
p-value that comes from the model and the time bins rather than from the animals. And every Dunn
comparison involving the blank condition rests on {P2['owner']['n_dogs']['C']} dogs, whose
dog-level floor is {dunn['D5']['floor']:.2g}.

**Effect sizes and intervals, which the paper does not report.** The pooled discrimination effect is
{p1['W2']['dog']['mean_diff']:+.2f} s per dog (95% BCa CI
[{p1['W2']['dog']['boot']['ci_low']:.2f}, {p1['W2']['dog']['boot']['ci_high']:.2f}], 43 dogs,
{N_BOOT:,} resamples), matched-pairs rank-biserial {p1['W2']['rank_biserial']:+.2f} at trial level.
Dishabituation is {p1['W3']['dog']['mean_diff']:+.2f} s
[{p1['W3']['dog']['boot']['ci_low']:.2f}, {p1['W3']['dog']['boot']['ci_high']:.2f}]. The
owner-proximity difference behind H4 is {dunn['D5']['mean_diff']:+.2f} score points between
{dunn['D5']['n1']} blank and {dunn['D5']['n2']} stress dogs. F5 shows the permutation null and the
bootstrap distribution for the discrimination claim side by side; the permutation resolution is
1/{N_PERM + 1:,} = {1 / (N_PERM + 1):.1e}, so no p-value from it should be read below that.

## 7. Multiplicity — F6

{img('F6_multiplicity', 'F6 multiplicity: the full p-value inventory and hits under each correction')}

{md_table(pd.DataFrame([dict(analysis=a, tests_run=n, reported_in_paper=rep, basis=b)
                        for a, n, rep, b in INVENTORY]))}

The script produces **{M_RUN} p-values**; the paper reports {M_PAPER}. At α = {ALPHA} over
{M_RUN} tests about {ALPHA * M_RUN:.1f} significant results are expected by chance alone. Recounting
the {len(recomputed_p)} recomputed values that are defined:

{md_table(recount['summary'])}

*{recount['note']}*

Within families, under Benjamini–Hochberg:

{md_table(pd.DataFrame([dict(family=f, tests=len(t), survive_bh=int(t['reject'].sum())) for f, t in bh_family.items()]))}

The paper states no correction anywhere. Which claims need one is a judgement about the claim: the
three pooled Part 1 comparisons are pre-specified by the paradigm and read as primary, while the
Part 2 battery — two behaviour categories, three proximity targets, three postures, then a sex
interaction on five of them — is a screen, and "dogs exhibited distinct behavioural responses
depending on the emotional content of the odor" is a claim about that family.

**Distribution of the p-values.** {jb['n_in_suspect_window']} of {len(recomputed_p)} fall in
[0.045, 0.05) against {jb['n_in_reference_window']} in [0.04, 0.045) (one-sided binomial
*p* = {jb['binomial_p_excess']:.3g}; the function reports this as
{'interpretable' if jb['interpretable'] else 'not interpretable'} at this scale). Three limits
apply and apply here: the diagnostic is a property of a *collection* of p-values and says nothing
about any individual result; it needs a reasonable number of them; and it cannot distinguish
selective reporting from legitimate exploratory analysis. **No intent is attributed here, and none
can be.** The inventory above — the full list of tests run, now visible because the script was
deposited — is what a reader needs, and depositing it is what made this recount possible.

## 8. Per-claim supplement

One panel per headline claim: the raw observations at the unit that was randomised, the permutation
null, and the three routes with the design floor shaded.

{chr(10).join(img('panel_' + h, h + ' ' + HEADLINE[h][0]) + chr(10) for h in HEADLINE)}
""")

P.append(f"""## 9. What was not assessed

{ac.uncheckable_assumptions_note()}

Specific to this review:

- **The figure images.** Figs 1–6 were not supplied, so no plotted value, error bar or significance
  bracket was checked against the text or the data.
- **Emotion induction.** Whether the three collection events produced the intended emotional states
  in the donor dog, and whether the samples differ chemically, are not statistical questions. The
  paper's validation is the owner's behavioural judgement of one donor.
- **Blinding and randomisation as executed.** Both are described in the methods (coder blind to
  condition, randomised pairing order and side); nothing in the deposit records the randomisation
  sequence, so the description was taken as given.
- **Video coding reliability.** No second coder, and no inter-rater agreement statistic, is reported
  or deposited; the per-dog scores were taken as given.
- **Why the 60-s sample point is missing for {288 - POSTURE_CELLS} dogs**, and whether missingness
  relates to the outcome. The README attributes "null" cells to early trial termination.
- **Nothing was subsampled or approximated.** Every claim was recomputed on the full deposited data.
  Where a reported p-value is a bound ("p < 0.001"), the recomputed value is given exactly and the
  bound is kept in `p_relation`.

## 10. Reproduction

| item | value |
|---|---|
| step 1 | `Rscript handoff/recompute_authors_tests.R` — the authors' own routines on their data, plus the dog-level and count-family alternatives |
| step 2 | `python review_{SLUG}.py` — routes, floors, intervals, figures, ledger, this report |
| seeds | `SEED = {SEED}`, `n_perm = {N_PERM:,}`, `n_boot = {N_BOOT:,}` (BCa) |
| permutation scheme | sign-flip within dog for paired Part 1 claims; label permutation between conditions for Part 2 claims; statistic = difference in means, two-sided |
| Part_1.xlsx sha256 | `{INPUTS['Part_1.xlsx']['sha256']}` |
| Part_2.xlsx sha256 | `{INPUTS['Part_2.xlsx']['sha256']}` |
| statistical_analyses.R sha256 | `{INPUTS['statistical_analyses.R']['sha256']}` |
| paper text sha256 | `{INPUTS['paper']['sha256']}` |
| R | {open(os.path.join(RRES, 'r_version.txt')).read().strip()} |
| python | {sys.version.split()[0]} |
| numpy / pandas / scipy / matplotlib | {np.__version__} / {pd.__version__} / {__import__('scipy').__version__} / {matplotlib.__version__} |

One defect in the deposited script, verified against the columns actually present:
`glm(Frequenc_MFy ~ Condition_MF * Sex_MF, data = box)` misspells the outcome
({'present' if TYPO else 'not present'} in the file; `box` contains `Frequency_MF`, not
`Frequenc_MFy`). The line immediately after it summarises the correctly spelled model, so no
reported number depends on it, but `source()`-ing the script stops there. Reproduction also needs
`install.packages` lines removed and the `dunn.test` dependency available; this review recomputes
Dunn's z directly instead, which is why both p-value conventions could be compared.

Derivations used above, each reproducible by hand: *V* ≤ n(n+1)/2 for a signed-rank statistic over
n non-zero pairs; the exact two-sided sign-flip floor 2/2ⁿ; the two-sided label-permutation floor
2/C(n₁+n₂, n₁); design effect 1 + (k−1)·ICC and effective n = rows / design effect.
""")

report = "\n".join(P)
open(os.path.join(SLUG, "REVIEW_REPORT.md"), "w", encoding="utf-8").write(report)

summary = dict(verdicts=VERDICTS, statcheck=SC, m_run=M_RUN, m_paper=M_PAPER,
               n_pvalues_recomputed=len(recomputed_p),
               recount={k: v for k, v in zip(recount["summary"]["correction"],
                                             recount["summary"]["n_significant"])},
               icc={k: round(v["icc"], 4) for k, v in INDEP.items()},
               effective_n={k: round(v["effective_n"], 1) for k, v in INDEP.items()},
               part1_dog={c: dict(trial_p=p1[c]["r_p"], dog_p=p1[c]["dog"]["p"]) for c in ["W1", "W2", "W3", "W4"]},
               part2_dog={c: dict(reported=omni[c]["p_reported"], recomputed=omni[c]["p_recomputed"],
                                  dog_kw=omni[c]["dog_kw_p"], dog_aov=omni[c]["dog_aov_p"],
                                  dog_qp=omni[c]["dog_qp_p"]) for c in omni},
               dunn={c: dict(reported=dunn[c]["p_reported"], half=dunn[c]["p_half"],
                             two=dunn[c]["p_two"], dog=dunn[c]["dog_wilcox_p"]) for c in dunn},
               glm={c: dict(reported=glm[c]["p_reported"], refit=glm[c]["p"], simple=glm[c]["simple_p"],
                            dog=glm[c]["dog_wilcox_p"], resid_df=glm[c]["resid_df"]) for c in glm},
               posture=dict(data=POSTURE_DATA, paper=POSTURE_PAPER, cells=POSTURE_CELLS),
               age=AGE, figures=FIGS, report_chars=len(report), jb=jb)
os.makedirs("handoff", exist_ok=True)
json.dump(summary, open("handoff/review_summary_v2.json", "w"), indent=1, default=str)
print(json.dumps({k: summary[k] for k in ["verdicts", "statcheck", "m_run", "m_paper",
                                          "n_pvalues_recomputed", "recount", "icc", "effective_n",
                                          "part1_dog", "report_chars"]}, indent=1, default=str))





