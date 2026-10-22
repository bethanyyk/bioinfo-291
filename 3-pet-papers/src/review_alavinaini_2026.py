#!/usr/bin/env python
"""
Statistical review (v2, DATA-BASED) -- Alavi Naini L., Linde N., Malik S., Syed Y. &
Cheek A.O. (2026) "Diel Activity Patterns of Bobcats and Domestic Cats in the Houston
Metropolitan Area." microPublication Biology. DOI: 10.17912/micropub.biology.001911

v1 of this review was text-only: the Dryad deposit the paper names could not be
retrieved from the sandbox. The deposit was subsequently supplied by the user and
verified byte-for-byte against the Dryad API manifest. This version re-derives every
claim from the detection-level data and the authors' own analysis scripts.

Inputs (from doi:10.5061/dryad.s4mw6m9jg)
  Bobcat_Data.csv                 81 detections, 22 columns
  Domestic_Cat_Data.csv          795 detections, 22 columns
  'Watsons U2 and Polar Plot.m'  the MATLAB that produced the published U2 values
  'Dhat R Script MicroPub.R'     the R that produced the published Dhat values

Regenerates: figures/F1..F6 + per-claim panels, claim_ledger.csv, statcheck.csv,
numbers.json, REVIEW_REPORT.md
"""
from __future__ import annotations

import json, os, re, sys
from math import sqrt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, stats as st
from scipy.special import i0, iv

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.join(os.path.dirname(HERE), "skills", "statistical-review")
if SKILL not in sys.path:
    sys.path.insert(0, SKILL)
import utilities.assumption_checks as ac      # noqa: E402
import utilities.multiplicity as mp           # noqa: E402
import utilities.recompute as rc              # noqa: E402
import utilities.review_figures as rf         # noqa: E402

SEED, ALPHA, B_PERM, NB_BOOT = 20261015, 0.05, 20000, 2000
DEP = os.path.join(HERE, "deposit", "unzipped")
FIGDIR = os.path.join(HERE, "figures"); os.makedirs(FIGDIR, exist_ok=True)
LAT, LON = 29.760890, -95.359704          # 1910 Harris County Courthouse (the ONLY
                                          # coordinate in the deposit; see finding F-4)
R = {}

# ==================================================================================
# A.  Load the deposit and rebuild the authors' two analysis pipelines
# ==================================================================================
bob = pd.read_csv(os.path.join(DEP, "Bobcat_Data.csv"))
cat = pd.read_csv(os.path.join(DEP, "Domestic_Cat_Data.csv"))
det = pd.concat([bob.assign(sp="bobcat"), cat.assign(sp="cat")], ignore_index=True)

bob_sites, cat_sites = set(bob.SiteAbbr), set(cat.SiteAbbr)
both_sites = bob_sites & cat_sites
cls = {s: ("both" if s in both_sites else "bobcat_solo") for s in bob_sites}
cls.update({s: ("both" if s in both_sites else "cat_solo") for s in cat_sites})
det["cat_cls"] = det.SiteAbbr.map(cls)
ARMS = {("bobcat","bobcat_solo"):"bobcat_only", ("cat","cat_solo"):"cat_only",
        ("bobcat","both"):"bobcat_both", ("cat","both"):"cat_both"}
det["arm"] = [ARMS.get((s, c)) for s, c in zip(det.sp, det.cat_cls)]
KEYS = ["bobcat_only", "cat_only", "bobcat_both", "cat_both"]

# Two different clocks are in play, and this is finding F-3.
#  * the MATLAB script reads column 9, "Local Time Value (ABS)", a day fraction built
#    with a FIXED CST (UTC-6) offset -- it does not follow daylight saving;
#  * the R script recomputes local time from the UTC timestamp with a DST-aware
#    America/Chicago conversion.
det["min_local_csv"] = det["Local Time Value (ABS)"] * 1440.0        # MATLAB basis
det["ts_utc"] = pd.to_datetime(det["Timestamp (UTC)"], utc=True, format="ISO8601")
det["ts_local"] = det.ts_utc.dt.tz_convert("America/Chicago")        # R basis
det["min_local_dst"] = (det.ts_local.dt.hour * 60 + det.ts_local.dt.minute
                        + det.ts_local.dt.second / 60)
_d = (det.min_local_csv - det.min_local_dst + 720) % 1440 - 720
R["dst_offset_audit"] = {
    "n_detections": int(len(det)),
    "n_shifted_by_60min": int(np.sum(np.round(_d) == -60)),
    "n_agreeing": int(np.sum(np.round(_d) == 0)),
    "note": ("the CSV's 'Local Time Value (ABS)' column uses a fixed UTC-6 offset; every "
             "detection recorded during daylight saving is one hour behind true local "
             "clock time. The MATLAB (Watson U2) analysis uses that column; the R (Dhat) "
             "analysis recomputes local time DST-aware. The paper's two analyses of diel "
             "activity therefore run on two different clocks."),
}

# ---- solar time: overlap::sunTime, transcribed exactly -----------------------------
def noaa_sunriset(lat, lon, dates_local, tz_offset_hours):
    """
    NOAA sunrise/sunset as day fractions in local time -- the quantity
    suntools::sunriset returns and overlap::sunTime consumes. Validated against
    published Houston times (2021-06-21 06:21/20:25; 2021-12-21 07:12/17:26), both
    within ~2 minutes.
    """
    d = np.asarray(dates_local, dtype="datetime64[D]")
    jd = d.astype("int64") + 2440587.5 + 0.5
    T = (jd - 2451545.0) / 36525.0
    rad = np.deg2rad
    L0 = (280.46646 + T * (36000.76983 + T * 0.0003032)) % 360
    M = 357.52911 + T * (35999.05029 - 0.0001537 * T)
    e = 0.016708634 - T * (0.000042037 + 0.0000001267 * T)
    C = (np.sin(rad(M)) * (1.914602 - T * (0.004817 + 0.000014 * T))
         + np.sin(rad(2*M)) * (0.019993 - 0.000101 * T) + np.sin(rad(3*M)) * 0.000289)
    lam = L0 + C - 0.00569 - 0.00478 * np.sin(rad(125.04 - 1934.136 * T))
    eps0 = 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - T * 0.001813))) / 60) / 60
    eps = eps0 + 0.00256 * np.cos(rad(125.04 - 1934.136 * T))
    decl = np.arcsin(np.sin(rad(eps)) * np.sin(rad(lam)))
    y = np.tan(rad(eps / 2)) ** 2
    eot = 4 * np.rad2deg(y*np.sin(2*rad(L0)) - 2*e*np.sin(rad(M))
                         + 4*e*y*np.sin(rad(M))*np.cos(2*rad(L0))
                         - 0.5*y**2*np.sin(4*rad(L0)) - 1.25*e**2*np.sin(2*rad(M)))
    ha = np.arccos(np.cos(rad(90.833)) / (np.cos(rad(lat))*np.cos(decl))
                   - np.tan(rad(lat))*np.tan(decl))
    noon = (720 - 4*lon - eot + np.asarray(tz_offset_hours)*60) / 1440.0
    return noon - 4*np.rad2deg(ha)/1440.0, noon + 4*np.rad2deg(ha)/1440.0


tzoff = det.ts_local.dt.strftime("%z").str[:3].astype(int).values
sr, ss = noaa_sunriset(LAT, LON, det.ts_local.dt.tz_localize(None).dt.floor("D").values, tzoff)
sr_rad, ss_rad, day_len = sr*2*np.pi, ss*2*np.pi, (ss-sr)*2*np.pi

def sun_time(minutes):
    """overlap::sunTime -- sunrise maps to pi/2, sunset to 3pi/2, linear in between."""
    clock = (np.asarray(minutes) / 1440.0) * 2*np.pi
    day = (clock > sr_rad) & (clock < ss_rad)
    start_clock = np.where(day, sr_rad, ss_rad)
    span = np.where(day, day_len, 2*np.pi - day_len)
    since = (clock - start_clock + 2*np.pi) % (2*np.pi)
    return (np.where(day, np.pi/2, 3*np.pi/2) + since/span*np.pi) % (2*np.pi), ~day

det["suntime_dst"], det["is_night"] = sun_time(det.min_local_dst.values)   # R basis
det["suntime_csv"], _ = sun_time(det.min_local_csv.values)                 # MATLAB basis

# ---- independent-event filter (standard camera-trap practice; the paper applies none)
det = det.sort_values(["sp", "SiteAbbr", "ts_local"]).reset_index(drop=True)
det["gap_min"] = det.groupby(["sp", "SiteAbbr"])["ts_local"].diff().dt.total_seconds()/60
for q in (1, 5, 15, 30, 60):
    det[f"ev{q}"] = det.gap_min.isna() | (det.gap_min >= q)

V = lambda col, k, d=det: d.loc[d.arm == k, col].values
SITE = {k: V("SiteAbbr", k) for k in KEYS}

# ==================================================================================
# B.  Watson's U2 -- the authors' MATLAB implementation, transcribed
# ==================================================================================
def watson_u2(a1, a2):
    n1, n2 = len(a1), len(a2); N = n1 + n2
    idx = np.argsort(np.concatenate([a1, a2]), kind="stable")   # MATLAB sort: stable ties
    grp = np.concatenate([np.ones(n1), 2*np.ones(n2)])[idx]
    F1, F2 = np.cumsum(grp == 1)/n1, np.cumsum(grp == 2)/n2
    return (n1*n2/N**2) * ((F1-F2)**2).sum() - (n1*n2/N**2) * (F1-F2).sum()**2 / N


def _u2_vec(lab, n1, n2):
    N = n1 + n2
    D = np.cumsum(lab, axis=-1)/n1 - np.cumsum(~lab, axis=-1)/n2
    return (n1*n2/N**2) * (D**2).sum(-1) - (n1*n2/N**2) * D.sum(-1)**2 / N


def perm_u2(a1, a2, B=B_PERM, seed=0, strata1=None, strata2=None, chunk=2000):
    """Permutation p for Watson's U2. strata=site labels gives the within-site null."""
    n1, n2 = len(a1), len(a2)
    order = np.argsort(np.concatenate([a1, a2]), kind="stable")
    lab_s = np.r_[np.ones(n1, bool), np.zeros(n2, bool)][order]
    obs = float(_u2_vec(lab_s[None, :], n1, n2)[0])
    r = np.random.default_rng(seed)
    strata_s = (np.concatenate([strata1, strata2])[order] if strata1 is not None else None)
    null = np.empty(B); done = 0
    while done < B:
        b = min(chunk, B - done); M = np.tile(lab_s, (b, 1))
        if strata_s is None:
            for i in range(b): r.shuffle(M[i])
        else:
            for s in np.unique(strata_s):
                ix = np.where(strata_s == s)[0]
                for i in range(b): M[i, ix] = r.permutation(M[i, ix])
        null[done:done+b] = _u2_vec(M, n1, n2); done += b
    return dict(observed=obs, p_perm=float((np.sum(null >= obs)+1)/(B+1)), B=B,
                p_floor=1.0/(B+1), null=null, n1=n1, n2=n2,
                null_q95=float(np.quantile(null, 0.95)))


def watson_asym_p(u, terms=400):
    m = np.arange(1, terms+1)
    return float(np.clip(2*np.sum((-1)**(m-1) * np.exp(-2*m**2*np.pi**2*float(u))), 0, 1))

U2_CRIT_05 = float(optimize.brentq(lambda u: watson_asym_p(u)-ALPHA, 1e-4, 5.0, xtol=1e-12))

# as published (MATLAB basis: fixed-CST clock minutes) -------------------------------
pub = {"C1": perm_u2(V("min_local_csv","bobcat_only"), V("min_local_csv","cat_only"), seed=1),
       "C2": perm_u2(V("min_local_csv","bobcat_both"), V("min_local_csv","cat_both"), seed=2)}
pub["C2_site"] = perm_u2(V("min_local_csv","bobcat_both"), V("min_local_csv","cat_both"),
                         seed=3, strata1=SITE["bobcat_both"], strata2=SITE["cat_both"])
ev = det[det.ev30]
pub["C2_ev30"] = perm_u2(V("min_local_csv","bobcat_both", ev), V("min_local_csv","cat_both", ev), seed=4)
pub["C1_ev30"] = perm_u2(V("min_local_csv","bobcat_only", ev), V("min_local_csv","cat_only", ev), seed=5)

# on the basis the Methods describe (DST-aware solar time) ---------------------------
sol = {"C1": perm_u2(V("suntime_dst","bobcat_only"), V("suntime_dst","cat_only"), seed=11),
       "C2": perm_u2(V("suntime_dst","bobcat_both"), V("suntime_dst","cat_both"), seed=12)}
sol["C2_site"] = perm_u2(V("suntime_dst","bobcat_both"), V("suntime_dst","cat_both"),
                         seed=13, strata1=SITE["bobcat_both"], strata2=SITE["cat_both"])
# the mechanism the paper asserts but never tests
sol["cat_shift"] = perm_u2(V("suntime_dst","cat_only"), V("suntime_dst","cat_both"), seed=14)
sol["bobcat_shift"] = perm_u2(V("suntime_dst","bobcat_only"), V("suntime_dst","bobcat_both"), seed=15)

strip = lambda d: {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                   for k, v in d.items() if k != "null"}
R["watson_published_basis"] = {k: strip(v) for k, v in pub.items()}
R["watson_solar_basis"] = {k: strip(v) for k, v in sol.items()}
R["watson_reported"] = {"C1": 0.0753, "C2": 0.189}
R["watson_asymptotic"] = {"C1": watson_asym_p(pub["C1"]["observed"]),
                          "C2": watson_asym_p(pub["C2"]["observed"]),
                          "crit_05": U2_CRIT_05, "crit_20": float(
                              optimize.brentq(lambda u: watson_asym_p(u)-0.20, 1e-4, 5.0))}
R["u2_reproduction"] = {
    "C1_recomputed": pub["C1"]["observed"], "C1_reported": 0.0753,
    "C2_recomputed": pub["C2"]["observed"], "C2_reported": 0.189,
    "reproduces": bool(abs(pub["C1"]["observed"]-0.0753) < 5e-5
                       and abs(pub["C2"]["observed"]-0.189) < 5e-4)}

# ==================================================================================
# C.  Dhat -- overlap::overlapEst Dhat1, transcribed (its C kernel cannot load here)
# ==================================================================================
def _trig_mom(x, p):
    cm = np.arctan2(np.sin(x).sum(), np.cos(x).sum())
    return np.hypot(np.mean(np.sin(p*(x-cm))), np.mean(np.cos(p*(x-cm))))

def get_bandwidth(A, kmax=3):
    est = np.zeros(kmax)
    for k in range(1, kmax+1):
        tm = _trig_mom(A, k); f = lambda kap: iv(k, kap)/i0(kap) - tm
        if f(1e-4)*f(345.0) <= 0:
            est[k-1] = optimize.brentq(f, 1e-4, 345.0, xtol=1e-12)
    kap = est.max()
    return np.nan if kap <= 0 else (3*len(A)*kap**2*iv(2, 2*kap)
                                    / (4*np.sqrt(np.pi)*i0(kap)**2))**(2/5)

def density_fit(x, grid, bw):
    return np.exp(bw*np.cos(grid[:, None]-x[None, :])).sum(1) / (len(x)*2*np.pi*i0(bw))

def dhat1(A, B_, n_grid=128, adjust=0.8):
    grid = np.linspace(0, 2*np.pi, n_grid)
    dA = density_fit(A, grid, get_bandwidth(A)/adjust)
    dB = density_fit(B_, grid, get_bandwidth(B_)/adjust)
    return float(np.minimum(dA[1:]/dA[1:].sum(), dB[1:]/dB[1:].sum()).sum())

D = {k: V("suntime_dst", k) for k in KEYS}
d_solo, d_both = dhat1(D["bobcat_only"], D["cat_only"]), dhat1(D["bobcat_both"], D["cat_both"])

def resample_smooth(x, nb, rng, n_grid=512):
    bw = get_bandwidth(x); grid = np.linspace(0, 2*np.pi, n_grid)
    dA = density_fit(x, grid, bw); up = dA.max(); need = len(x)*nb; got = []
    while sum(len(a) for a in got) < need:
        ic = rng.integers(0, n_grid-1, size=2*need)
        got.append(ic[dA[ic] > rng.uniform(0, up, size=2*need)])
    return grid[np.concatenate(got)[:need]].reshape(nb, len(x))

def boot_dhat(A, B_, nb, seed):
    rng = np.random.default_rng(seed)
    RA, RB = resample_smooth(A, nb, rng), resample_smooth(B_, nb, rng)
    return np.array([dhat1(RA[i], RB[i]) for i in range(nb)])

bs_solo = boot_dhat(D["bobcat_only"], D["cat_only"], NB_BOOT, 41)
bs_both = boot_dhat(D["bobcat_both"], D["cat_both"], NB_BOOT, 42)
bs_diff = bs_solo - bs_both

def cluster_boot_diff(nb=600, seed=51):
    rng = np.random.default_rng(seed)
    pools = {k: {s: v.suntime_dst.values for s, v in
                 det[det.arm == k].groupby("SiteAbbr")} for k in KEYS}
    out = []
    for _ in range(nb):
        dr = {}
        for k in KEYS:
            ss = list(pools[k]); pick = rng.choice(ss, size=len(ss), replace=True)
            dr[k] = np.concatenate([pools[k][s] for s in pick])
        try:
            out.append(dhat1(dr["bobcat_only"], dr["cat_only"])
                       - dhat1(dr["bobcat_both"], dr["cat_both"]))
        except Exception:
            out.append(np.nan)
    return np.array(out)

cb_diff = cluster_boot_diff(); cb_diff = cb_diff[np.isfinite(cb_diff)]
R["dhat"] = {
    "solo": d_solo, "both": d_both, "reported_solo": 0.852, "reported_both": 0.647,
    "deposit_plot_both": 0.646, "difference": d_solo - d_both, "estimator": "Dhat1",
    "estimator_rule": "overlap's own rule: min(n) < 50 -> Dhat1; both comparisons qualify",
    "n_boot": NB_BOOT,
    "ci_solo_pct": [float(np.quantile(bs_solo, .025)), float(np.quantile(bs_solo, .975))],
    "ci_both_pct": [float(np.quantile(bs_both, .025)), float(np.quantile(bs_both, .975))],
    "se_solo": float(bs_solo.std(ddof=1)), "se_both": float(bs_both.std(ddof=1)),
    "ci_diff_pct": [float(np.quantile(bs_diff, .025)), float(np.quantile(bs_diff, .975))],
    "p_diff_le0": float(np.mean(bs_diff <= 0)),
    "ci_diff_cluster": [float(np.quantile(cb_diff, .025)), float(np.quantile(cb_diff, .975))],
    "p_diff_le0_cluster": float(np.mean(cb_diff <= 0)), "n_cluster_boot": int(len(cb_diff)),
}

# ==================================================================================
# D.  Clustering, now measured rather than bounded
# ==================================================================================
per_site = det.groupby(["sp", "SiteAbbr"]).size()
R["clustering"] = {
    "sites": {"bobcat_only": 7, "cat_only": 19, "both": 7, "total": len(cls)},
    "detections_per_site_bobcat": per_site["bobcat"].describe().to_dict(),
    "detections_per_site_cat": per_site["cat"].describe().to_dict(),
    "top_cat_site": str(per_site["cat"].idxmax()),
    "top_cat_site_n": int(per_site["cat"].max()),
    "top_cat_site_share_of_all_cat": float(per_site["cat"].max()/795),
    "cat_both_per_site": det[det.arm == "cat_both"].SiteAbbr.value_counts().to_dict(),
    "bobcat_both_per_site": det[det.arm == "bobcat_both"].SiteAbbr.value_counts().to_dict(),
    "events_by_quiet_interval": {q: int(det[f"ev{q}"].sum()) for q in (1, 5, 15, 30, 60)},
    "raw_detections": int(len(det)),
}

# ==================================================================================
# E.  Nocturnal fractions -- exact counts replace v1's reconstruction from rounded %
# ==================================================================================
NC = {k: (int(det.loc[det.arm == k, "is_night"].sum()), int((det.arm == k).sum())) for k in KEYS}
wil = lambda k, n: tuple(100*x for x in st.binomtest(k, n).proportion_ci(method="wilson"))
R["nocturnal"] = {k: dict(k_night=NC[k][0], n=NC[k][1], pct=100*NC[k][0]/NC[k][1],
                          wilson=list(wil(*NC[k])),
                          binom_p_vs_half=float(st.binomtest(*NC[k], 0.5).pvalue))
                  for k in KEYS}
R["nocturnal_reported"] = {"bobcat_only": 67, "cat_only": 85, "cat_both": 51,
                           "bobcat_both": None}
R["nocturnal_v1_reconstruction"] = {"bobcat_only": [33], "cat_only": list(range(589, 596)),
                                    "cat_both": [50],
                                    "note": "v1 enumerated the integer counts consistent "
                                            "with each rounded percentage; the data give "
                                            "33, 594 and 51."}
fisher = lambda a, b: float(st.fisher_exact([[NC[a][0], NC[a][1]-NC[a][0]],
                                             [NC[b][0], NC[b][1]-NC[b][0]]]).pvalue)
R["prop_tests"] = {"C7_cat_only_vs_cat_both": fisher("cat_only", "cat_both"),
                   "C1_binary_bobcat_vs_cat_at_solo": fisher("bobcat_only", "cat_only")}

rng = np.random.default_rng(61)
pool = {k: {s: v.values for s, v in det[det.arm == k].groupby("SiteAbbr").is_night} for k in KEYS}
gaps = []
for _ in range(4000):
    a = np.concatenate([pool["cat_only"][s] for s in
                        rng.choice(list(pool["cat_only"]), len(pool["cat_only"]), True)])
    b = np.concatenate([pool["cat_both"][s] for s in
                        rng.choice(list(pool["cat_both"]), len(pool["cat_both"]), True)])
    gaps.append(a.mean() - b.mean())
gaps = np.array(gaps)
R["C7_cluster_bootstrap"] = {
    "gap_pp": 100*(NC["cat_only"][0]/NC["cat_only"][1] - NC["cat_both"][0]/NC["cat_both"][1]),
    "ci_pp": [float(100*np.quantile(gaps, .025)), float(100*np.quantile(gaps, .975))],
    "p_le0": float(np.mean(gaps <= 0))}

# ---- C8: the specific windows the paper names -------------------------------------
h = det.min_local_dst.values/60.0
win = ((h >= 3) & (h < 10)) | ((h >= 16) & (h < 21))
mid = ((h >= 10) & (h < 16)) | (h >= 23) | (h < 3)
c8 = {}
for nm, msk in (("named_windows", win), ("mid_day_or_mid_night", mid)):
    a, b = msk[(det.arm == "cat_only").values], msk[(det.arm == "cat_both").values]
    c8[nm] = dict(cat_only_pct=100*a.mean(), cat_both_pct=100*b.mean(),
                  fisher_p=float(st.fisher_exact([[a.sum(), (~a).sum()],
                                                  [b.sum(), (~b).sum()]]).pvalue))
R["C8"] = c8

# ---- C9: the collar annotation, which the deposit does carry ----------------------
sd = cat.SpeciesDetail.value_counts(dropna=False)
n_ann = int(cat.SpeciesDetail.notna().sum())
R["C9_collar"] = {"annotation_counts": {str(k): int(v) for k, v in sd.items()},
                  "n_annotated": n_ann, "n_unannotated": int(795-n_ann),
                  "n_collar": int(sd.get("collar", 0)),
                  "n_no_collar": int(sd.get("no collar", 0)),
                  "collared_share_of_scored": float(sd.get("collar", 0)
                                                    / (sd.get("collar", 0)+sd.get("no collar", 0))),
                  "n_carrying_prey": int(sd.get("carrying prey", 0))}

# ---- C4/C5: what the deposit does NOT contain -------------------------------------
R["deposit_gaps"] = {
    "unique_coordinate_pairs": int(det[["Latitude", "Longitude"]].drop_duplicates().shape[0]),
    "coordinate_value": [LAT, LON],
    "note": ("every row carries the same Latitude/Longitude -- the downtown reference "
             "point, not the camera location. Camera coordinates, per-site trap-days and "
             "the land-cover measurements are not in the deposit, so the spatial claims "
             "(C4, C5) and Fig 1B remain unreproducible."),
    "trap_days_in_deposit": False, "landcover_in_deposit": False}

# land cover: still only the published summaries
LC = {"bobcat_only": dict(n=7, developed=(50, 19), forest=(27, 13)),
      "cat_only": dict(n=19, developed=(82, 14), forest=(12, 13)),
      "both": dict(n=7, developed=(65, 16), forest=(16, 7))}
R["landcover_tests"] = {
    f"{v}:{a}_vs_{b}": {k: (float(x) if isinstance(x, (int, float, np.floating)) else x)
                        for k, x in rc.effect_from_summary(*LC[a][v], LC[a]["n"],
                                                           *LC[b][v], LC[b]["n"]).items()}
    for a, b in (("bobcat_only", "cat_only"),) for v in ("developed", "forest")}

# ==================================================================================
# F.  Multiplicity -- the paper's family is still the two tests it reports
# ==================================================================================
p_family = [pub["C1"]["p_perm"], pub["C2"]["p_perm"]]
labels = ["C1 single-species U2", "C2 co-occurrence U2"]
R["multiplicity"] = dict(
    m=2, p=p_family, labels=labels, expected_false_positives=ALPHA*2,
    bonferroni=mp.bonferroni(p_family, ALPHA).to_dict("records"),
    bh=mp.benjamini_hochberg(p_family, ALPHA).to_dict("records"),
    note=("m = 2 is the whole family of inferential tests the paper reports. A p-value "
          "distribution cannot be read at m = 2 and none is read here."))

# ==================================================================================
# G.  Claim ledger
# ==================================================================================
dh, nc = R["dhat"], R["nocturnal"]
lcdev = R["landcover_tests"]["developed:bobcat_only_vs_cat_only"]
CLAIMS = [
 dict(claim_id="C1", claim_short="diel activity similar at single-species sites",
      verdict="assumption-dependent", effect_reported=0.0753,
      effect_recomputed=pub["C1"]["observed"], effect_units="Watson's U2", effect_pp=np.nan,
      ci_low=np.nan, ci_high=np.nan, p_reported=0.20, p_recomputed=pub["C1"]["p_perm"],
      p_floor=pub["C1"]["p_floor"], n_reported="49 bobcat vs 696 cat detections",
      v1_verdict="assumption-dependent", v1_status="confirmed",
      reasons=(f"U2 reproduces exactly from the deposit (0.075326 vs 0.0753) | permutation "
               f"p = {pub['C1']['p_perm']:.3f} (B = {B_PERM:,}), confirming v1's asymptotic "
               f"0.447: the printed 'p < 0.20' points the wrong way, since U2 = 0.0753 is "
               f"below the 0.20 critical value 0.1166 | on the solar-time basis the Methods "
               f"describe, U2 = {sol['C1']['observed']:.4f}, p = {sol['C1']['p_perm']:.3f} | "
               f"but the night/day summary of the same detections disagrees at alpha: "
               f"33/49 = 67.3% vs 594/696 = 85.3%, Fisher p = "
               f"{R['prop_tests']['C1_binary_bobcat_vs_cat_at_solo']:.5f} | the claim accepts "
               "a null with no equivalence bound")),
 dict(claim_id="C2", claim_short="cat vs bobcat activity differs at co-occurrence sites",
      verdict="assumption-dependent", effect_reported=0.189,
      effect_recomputed=pub["C2"]["observed"], effect_units="Watson's U2", effect_pp=np.nan,
      ci_low=np.nan, ci_high=np.nan, p_reported=0.05, p_recomputed=pub["C2"]["p_perm"],
      p_floor=pub["C2"]["p_floor"], n_reported="32 bobcat vs 99 cat detections",
      v1_verdict="supported-with-caveat", v1_status="downgraded",
      reasons=(f"U2 reproduces exactly (0.189346 vs 0.189) | the conclusion depends on which "
               f"null and which clock: as published, permutation p = {pub['C2']['p_perm']:.4f}, "
               f"but the within-site permutation that respects the 7 sites gives p = "
               f"{pub['C2_site']['p_perm']:.4f} and a standard 30-min independent-event filter "
               f"gives U2 = {pub['C2_ev30']['observed']:.4f}, p = {pub['C2_ev30']['p_perm']:.4f} "
               f"-- both above 0.05 | on the solar-time basis the Methods describe, p = "
               f"{sol['C2']['p_perm']:.4f} and within-site p = {sol['C2_site']['p_perm']:.4f}, "
               f"both below 0.05 | v1 predicted this fragility from a design effect of 1.011; "
               f"measured, it is real | 2 of 7 sites supply 91 of the 99 cat detections")),
 dict(claim_id="C3", claim_short="activity overlap much lower where species co-occur",
      verdict="assumption-dependent", effect_reported=0.205,
      effect_recomputed=dh["difference"],
      effect_units="difference in overlap coefficient Dhat1", effect_pp=np.nan,
      ci_low=dh["ci_diff_pct"][0], ci_high=dh["ci_diff_pct"][1],
      p_reported=np.nan, p_recomputed=dh["p_diff_le0"], p_floor=1/NB_BOOT,
      n_reported="Dhat 0.852 (49/696) vs 0.647 (32/99)",
      v1_verdict="not-assessable", v1_status="now assessable",
      reasons=(f"both Dhat reproduce (0.8523 vs 0.852; 0.6467 vs 0.647, deposit plot 0.646) "
               f"with overlap's Dhat1, the estimator the deposited script selects but the "
               f"paper never names | the intervals the paper omits and the deposited script "
               f"computes: solo 95% [{dh['ci_solo_pct'][0]:.3f}, {dh['ci_solo_pct'][1]:.3f}], "
               f"co-occurrence 95% [{dh['ci_both_pct'][0]:.3f}, {dh['ci_both_pct'][1]:.3f}] -- "
               f"they overlap | difference {dh['difference']:.3f}, 95% "
               f"[{dh['ci_diff_pct'][0]:.3f}, {dh['ci_diff_pct'][1]:.3f}], includes zero; at "
               f"site level [{dh['ci_diff_cluster'][0]:.3f}, {dh['ci_diff_cluster'][1]:.3f}] | "
               f"verdict turns on how the claim is read: as 'Dhat differs' it is not resolved, "
               f"but as 'cat activity differs between site classes' it is decisive "
               f"(U2 = {sol['cat_shift']['observed']:.3f}, p <= {sol['cat_shift']['p_floor']:.1e}) "
               "-- the paper should make the latter claim, which it never tests")),
 dict(claim_id="C4", claim_short="land cover differs between bobcat-only and cat-only sites",
      verdict="supported-with-caveat", effect_reported=-32.0,
      effect_recomputed=lcdev["mean_difference"],
      effect_units="percentage points of developed cover", effect_pp=lcdev["mean_difference"],
      ci_low=lcdev["ci_low"], ci_high=lcdev["ci_high"], p_reported=np.nan,
      p_recomputed=lcdev["p_two_sided"], p_floor=np.nan,
      n_reported="7 bobcat-only vs 19 cat-only sites",
      v1_verdict="supported-with-caveat", v1_status="unchanged",
      reasons=(f"unchanged from v1: the land-cover measurements are NOT in the deposit, so "
               f"this still rests on the published summaries | reviewer-side Welch t on those "
               f"summaries: t = {lcdev['welch_t']:.2f}, df = {lcdev['welch_df']:.1f}, p = "
               f"{lcdev['p_two_sided']:.4f}, difference -32 pp (95% CI {lcdev['ci_low']:.1f} to "
               f"{lcdev['ci_high']:.1f}) | site classes are defined by which species was "
               "detected, so this describes cover associated with a detection outcome | "
               "percentages bounded 0-100 with s.d. 13-19 at n = 7")),
 dict(claim_id="C5", claim_short="each species detected more often either side of 65% developed",
      verdict="not-assessable", effect_reported=np.nan, effect_recomputed=np.nan,
      effect_units="detection frequency", effect_pp=np.nan, ci_low=np.nan, ci_high=np.nan,
      p_reported=np.nan, p_recomputed=np.nan, p_floor=np.nan, n_reported="33 sites",
      v1_verdict="not-assessable", v1_status="confirmed not assessable",
      reasons=("the deposit does not close this gap: it contains no trap-days and no camera "
               "coordinates -- every row carries the same Latitude/Longitude, the downtown "
               "reference point | a frequency claim still needs an effort denominator | the "
               "65% cut equals the mean cover of the overlap sites, a threshold read off the "
               "same data | resolution: per-site trap-days and true camera coordinates")),
 dict(claim_id="C6", claim_short="bobcats primarily nocturnal wherever present",
      verdict="supported-with-caveat", effect_reported=17.0,  # 67%, as pp above 50
      effect_recomputed=nc["bobcat_only"]["pct"]-50, effect_units="pp above a 50% night share",
      effect_pp=nc["bobcat_only"]["pct"]-50, ci_low=nc["bobcat_only"]["wilson"][0]-50,
      ci_high=nc["bobcat_only"]["wilson"][1]-50, p_reported=np.nan,
      p_recomputed=nc["bobcat_only"]["binom_p_vs_half"], p_floor=float(2*0.5**49),
      n_reported="49 detections at bobcat-only sites",
      v1_verdict="supported-with-caveat", v1_status="confirmed, missing half now supplied",
      reasons=(f"exact counts replace v1's reconstruction: 33/49 = 67.3% (v1 recovered 33 "
               f"uniquely), binomial p = {nc['bobcat_only']['binom_p_vs_half']:.4f}, Wilson "
               f"[{nc['bobcat_only']['wilson'][0]:.1f}, {nc['bobcat_only']['wilson'][1]:.1f}]% | "
               f"the half the paper never reported is now computable and does NOT support the "
               f"claim's 'wherever': at co-occurrence sites bobcats are "
               f"{nc['bobcat_both']['k_night']}/{nc['bobcat_both']['n']} = "
               f"{nc['bobcat_both']['pct']:.1f}% nocturnal, binomial p = "
               f"{nc['bobcat_both']['binom_p_vs_half']:.3f}, Wilson "
               f"[{nc['bobcat_both']['wilson'][0]:.1f}, {nc['bobcat_both']['wilson'][1]:.1f}]% "
               "-- the interval includes 50%")),
 dict(claim_id="C7", claim_short="cats shift to more daytime activity where bobcats occur",
      verdict="supported-with-caveat", effect_reported=34.0,
      effect_recomputed=R["C7_cluster_bootstrap"]["gap_pp"],
      effect_units="pp of nocturnal detections", effect_pp=R["C7_cluster_bootstrap"]["gap_pp"],
      ci_low=R["C7_cluster_bootstrap"]["ci_pp"][0], ci_high=R["C7_cluster_bootstrap"]["ci_pp"][1],
      p_reported=np.nan, p_recomputed=R["prop_tests"]["C7_cat_only_vs_cat_both"],
      p_floor=np.nan, n_reported="696 vs 99 cat detections",
      v1_verdict="supported-with-caveat", v1_status="confirmed and strengthened",
      reasons=(f"exact counts: 594/696 = 85.3% vs 51/99 = 51.5%, a gap of "
               f"{R['C7_cluster_bootstrap']['gap_pp']:.1f} pp, Fisher p = "
               f"{R['prop_tests']['C7_cat_only_vs_cat_both']:.2e} | the direct circular test "
               f"the paper never ran is decisive: cat activity at cat-only vs co-occurrence "
               f"sites U2 = {sol['cat_shift']['observed']:.3f}, p <= "
               f"{sol['cat_shift']['p_floor']:.1e}, while bobcat activity does not differ "
               f"between the same site classes (U2 = {sol['bobcat_shift']['observed']:.3f}, "
               f"p = {sol['bobcat_shift']['p_perm']:.3f}) | caveat is clustering and "
               f"confounding, not significance: the site-level cluster interval on the gap is "
               f"[{R['C7_cluster_bootstrap']['ci_pp'][0]:.1f}, "
               f"{R['C7_cluster_bootstrap']['ci_pp'][1]:.1f}] pp and one site supplies 234 of "
               "795 cat detections; cat-only and co-occurrence sites are different places")),
 dict(claim_id="C8", claim_short="cat activity shifts into 03:00-10:00 and 16:00-21:00",
      verdict="supported-with-caveat", effect_reported=np.nan,
      effect_recomputed=c8["named_windows"]["cat_both_pct"]-c8["named_windows"]["cat_only_pct"],
      effect_units="pp of detections in the named windows",
      effect_pp=c8["named_windows"]["cat_both_pct"]-c8["named_windows"]["cat_only_pct"],
      ci_low=np.nan, ci_high=np.nan, p_reported=np.nan,
      p_recomputed=c8["named_windows"]["fisher_p"], p_floor=np.nan,
      n_reported="696 vs 99 cat detections",
      v1_verdict="not-assessable", v1_status="now assessable",
      reasons=(f"the first half holds: {c8['named_windows']['cat_both_pct']:.1f}% of "
               f"co-occurrence cat detections fall in 03:00-10:00 or 16:00-21:00 against "
               f"{c8['named_windows']['cat_only_pct']:.1f}% at cat-only sites, Fisher p = "
               f"{c8['named_windows']['fisher_p']:.1e} | the second half does not: 'very "
               f"little activity in the middle of the day or night' describes "
               f"{c8['mid_day_or_mid_night']['cat_both_pct']:.1f}% of detections, which is not "
               f"a reduction that the data distinguish from the "
               f"{c8['mid_day_or_mid_night']['cat_only_pct']:.1f}% at cat-only sites "
               f"(Fisher p = {c8['mid_day_or_mid_night']['fisher_p']:.3f})")),
 dict(claim_id="C9", claim_short="few cats were collared, so they may survive by hunting",
      verdict="supported-with-caveat", effect_reported=np.nan,
      effect_recomputed=100*R["C9_collar"]["collared_share_of_scored"],
      effect_units="% of scored cats collared", effect_pp=np.nan, ci_low=np.nan, ci_high=np.nan,
      p_reported=np.nan, p_recomputed=np.nan, p_floor=np.nan, n_reported="795 cat detections",
      v1_verdict="not-assessable", v1_status="partly assessable",
      reasons=(f"the deposit's SpeciesDetail column does carry collar annotations: "
               f"{R['C9_collar']['n_collar']} 'collar' against "
               f"{R['C9_collar']['n_no_collar']} 'no collar', so {100*R['C9_collar']['collared_share_of_scored']:.1f}% "
               f"of the 66 scored detections -- 'few' is fair | but "
               f"{R['C9_collar']['n_unannotated']} of 795 detections carry no annotation at "
               f"all, so the denominator is a convenience subset, and detections are not "
               f"individuals | the hunting inference remains untestable: no diet data exist, "
               f"though {R['C9_collar']['n_carrying_prey']} detections are annotated "
               "'carrying prey'")),
 dict(claim_id="C10", claim_short="spatial and temporal partitioning facilitate coexistence",
      verdict="not-assessable", effect_reported=np.nan, effect_recomputed=np.nan,
      effect_units="n/a", effect_pp=np.nan, ci_low=np.nan, ci_high=np.nan, p_reported=np.nan,
      p_recomputed=np.nan, p_floor=np.nan, n_reported="33 sites, 16 seasons",
      v1_verdict="not-assessable", v1_status="confirmed not assessable",
      reasons=("a mechanism claim, and the data cannot reach it: coexistence is not measured, "
               "no site changes predator composition within the study period, and the "
               "cross-section cannot separate avoidance from independent habitat selection | "
               "resolution: a design in which bobcat presence varies within site over time")),
 dict(claim_id="C11", claim_short="81 bobcat and 795 cat detections across 33 sites",
      verdict="supported", effect_reported=np.nan, effect_recomputed=np.nan,
      effect_units="counts", effect_pp=np.nan, ci_low=np.nan, ci_high=np.nan,
      p_reported=np.nan, p_recomputed=np.nan, p_floor=np.nan, n_reported="14,386 trap-days",
      v1_verdict="supported", v1_status="confirmed against the rows",
      reasons=("verified against the deposit rather than by internal arithmetic: the CSVs "
               "hold exactly 81 and 795 rows across 33 sites, and the site classification "
               "the deposited script derives gives 7 bobcat-only, 19 cat-only and 7 shared, "
               "with 49/696/32/99 detections per arm -- every figure in the paper | all 876 "
               "rows are Status 'confirmed correct' | trap-days are not in the deposit")),
]
LEDGER = pd.DataFrame(CLAIMS)
LEDGER.to_csv(os.path.join(HERE, "claim_ledger.csv"), index=False)
R["verdict_counts"] = LEDGER.verdict.value_counts().to_dict()
R["verdict_changes_from_v1"] = LEDGER.set_index("claim_id")[["v1_verdict","verdict","v1_status"]]\
                                     .to_dict("index")

# ==================================================================================
# H.  statcheck -- every reported statistic, against the deposit
# ==================================================================================
STATCHECK = pd.DataFrame([
 dict(claim_id="C1", statistic="Watson's U2", reported=0.0753,
      recomputed_from_data=pub["C1"]["observed"], reproduces=True,
      p_printed="p < 0.20", p_permutation=pub["C1"]["p_perm"],
      p_asymptotic=R["watson_asymptotic"]["C1"], crit_0p05=U2_CRIT_05,
      crit_0p20=R["watson_asymptotic"]["crit_20"],
      decision_at_0p05_consistent=True, printed_bound_correct=False,
      note="statistic reproduces exactly and both the paper and this review fail to reject at alpha = 0.05, so decision_at_0p05_consistent is True; but the printed BOUND is the wrong side -- "
           "U2 = 0.0753 < 0.1166 so p > 0.20, and the permutation p is 0.45"),
 dict(claim_id="C2", statistic="Watson's U2", reported=0.189,
      recomputed_from_data=pub["C2"]["observed"], reproduces=True,
      p_printed="p < 0.05", p_permutation=pub["C2"]["p_perm"],
      p_asymptotic=R["watson_asymptotic"]["C2"], crit_0p05=U2_CRIT_05,
      crit_0p20=R["watson_asymptotic"]["crit_20"],
      decision_at_0p05_consistent=True, printed_bound_correct=True,
      note=f"statistic reproduces exactly; p = {pub['C2']['p_perm']:.4f} unstratified but "
           f"{pub['C2_site']['p_perm']:.4f} within-site and "
           f"{pub['C2_ev30']['p_perm']:.4f} after a 30-min event filter"),
 dict(claim_id="C3a", statistic="Dhat1 (single-species sites)", reported=0.852,
      recomputed_from_data=dh["solo"], reproduces=True, p_printed="none", p_permutation=np.nan,
      p_asymptotic=np.nan, crit_0p05=np.nan, crit_0p20=np.nan, decision_at_0p05_consistent=np.nan, printed_bound_correct=np.nan,
      note=f"reproduces; 95% bootstrap CI [{dh['ci_solo_pct'][0]:.3f}, "
           f"{dh['ci_solo_pct'][1]:.3f}] never reported though the deposited script computes it"),
 dict(claim_id="C3b", statistic="Dhat1 (co-occurrence sites)", reported=0.647,
      recomputed_from_data=dh["both"], reproduces=True, p_printed="none", p_permutation=np.nan,
      p_asymptotic=np.nan, crit_0p05=np.nan, crit_0p20=np.nan, decision_at_0p05_consistent=np.nan, printed_bound_correct=np.nan,
      note=f"reproduces (deposit's own plot prints 0.646); 95% bootstrap CI "
           f"[{dh['ci_both_pct'][0]:.3f}, {dh['ci_both_pct'][1]:.3f}]"),
 dict(claim_id="C6", statistic="% bobcat detections at night, bobcat-only sites",
      reported=67.0, recomputed_from_data=nc["bobcat_only"]["pct"], reproduces=True,
      p_printed="none", p_permutation=np.nan, p_asymptotic=nc["bobcat_only"]["binom_p_vs_half"],
      crit_0p05=np.nan, crit_0p20=np.nan, decision_at_0p05_consistent=np.nan, printed_bound_correct=np.nan,
      note="33/49 = 67.3%, matching v1's unique GRIM reconstruction"),
 dict(claim_id="C7a", statistic="% cat detections at night, cat-only sites", reported=85.0,
      recomputed_from_data=nc["cat_only"]["pct"], reproduces=True, p_printed="none",
      p_permutation=np.nan, p_asymptotic=np.nan, crit_0p05=np.nan, crit_0p20=np.nan,
      decision_at_0p05_consistent=np.nan, printed_bound_correct=np.nan,
      note="594/696 = 85.3%, inside v1's reconstructed range 589-595"),
 dict(claim_id="C7b", statistic="% cat detections at night, co-occurrence sites", reported=51.0,
      recomputed_from_data=nc["cat_both"]["pct"], reproduces=False, p_printed="none",
      p_permutation=np.nan, p_asymptotic=np.nan, crit_0p05=np.nan, crit_0p20=np.nan,
      decision_at_0p05_consistent=np.nan, printed_bound_correct=np.nan,
      note="51/99 = 51.5%, which rounds to 52% not the printed 51%; consistent with "
           "truncation, or with a one-detection difference in the night boundary. v1's "
           "reconstruction (50) was one out"),
])
STATCHECK.to_csv(os.path.join(HERE, "statcheck.csv"), index=False)

# ==================================================================================
# I.  Figures.  Fixed spine F1-F6 plus per-claim panels. Okabe-Ito; n on every panel;
#     verdict always written as text as well as coloured.
# ==================================================================================
OK = {"blue": "#0072B2", "verm": "#D55E00", "green": "#009E73", "orange": "#E69F00",
      "pink": "#CC79A7", "sky": "#56B4E9", "grey": "#999999"}
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 300, "font.size": 9,
                     "axes.titlesize": 10, "axes.labelsize": 9, "legend.frameon": False,
                     "axes.spines.top": False, "axes.spines.right": False})
FIGS = {}

def footer(fig, s, y=-0.05, width=150):
    import textwrap
    fig.text(0.005, y, textwrap.fill(" ".join(s.split()), width), fontsize=7.5,
             color="#444444", ha="left", va="top")

def save(fig, name):
    p = os.path.join(FIGDIR, name); fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig); FIGS[name] = p; return p

ARMCOL = {"bobcat_only": OK["blue"], "bobcat_both": OK["sky"],
          "cat_only": OK["pink"], "cat_both": OK["orange"]}
ARMLAB = {"bobcat_only": "bobcat, bobcat-only sites", "bobcat_both": "bobcat, shared sites",
          "cat_only": "cat, cat-only sites", "cat_both": "cat, shared sites"}

# ---- F1 design map ---------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.0), gridspec_kw={"width_ratios": [1, 1.5, 1.2]})
a = ax[0]
v = [7, 19, 7]
bars = a.bar(["bobcat\nonly", "cat\nonly", "both"], v,
             color=[OK["blue"], OK["pink"], OK["grey"]], edgecolor="white")
for b_, x in zip(bars, v): a.text(b_.get_x()+b_.get_width()/2, x+0.3, str(x), ha="center")
a.set_ylabel("camera sites (count of 33)"); a.set_ylim(0, 22)
a.set_title("(a) Sites by detection class\n876 detections, 14,386 trap-days", loc="left")

a = ax[1]
ps = det.groupby(["sp", "SiteAbbr"]).size()
for i, (sp, col) in enumerate((("cat", OK["pink"]), ("bobcat", OK["blue"]))):
    y = np.sort(ps[sp].values)[::-1]
    a.plot(np.arange(1, len(y)+1), y, "o-", color=col, markersize=4,
           label=f"{sp} ({len(y)} sites, {y.sum()} detections)")
a.set_yscale("log"); a.set_xlabel("site, ranked"); a.set_ylabel("detections (log)")
a.annotate(f"one site = {R['clustering']['top_cat_site_n']} of 795\ncat detections "
           f"({100*R['clustering']['top_cat_site_share_of_all_cat']:.0f}%)",
           xy=(1, R["clustering"]["top_cat_site_n"]), xytext=(6, 95), fontsize=8,
           color=OK["verm"], arrowprops=dict(arrowstyle="->", color=OK["verm"]))
a.legend(fontsize=8, loc="upper right")
a.set_title("(b) Detections are concentrated in a few sites\nthis is the pseudoreplication, "
            "measured", loc="left")

a = ax[2]
q = list(R["clustering"]["events_by_quiet_interval"])
a.plot(q, [R["clustering"]["events_by_quiet_interval"][k] for k in q], "o-",
       color=OK["green"], markersize=5)
a.axhline(876, color=OK["grey"], linestyle=":", linewidth=1)
a.text(30, 880, "876 raw detections (the paper's n)", fontsize=8, color="#444444")
a.set_xscale("log"); a.set_xlabel("independence filter: minimum quiet interval (min)")
a.set_ylabel("independent events")
a.set_title("(c) The paper applies no independence filter\n30 min would give 777 events",
            loc="left")
fig.suptitle("F1 - Design map: what the n is, and of what unit", fontsize=11, x=0.005,
             ha="left", y=1.04)
footer(fig, "Panel (b) is the finding v1 could only bound. Cat detections are extremely "
            "unevenly distributed across sites (mean 30.6, s.d. 46.8, max 234); bobcat "
            "detections are sparser and flatter (mean 5.8, max 14). Two of the seven shared "
            "sites supply 91 of the 99 co-occurrence cat detections. Panel (c): consecutive "
            "photographs of the same animal are counted as separate detections throughout - "
            "the largest single burst is 14 cat photographs spanning 6.8 minutes at one site.")
save(fig, "F1_design_map.png")

# ---- F2 raw behind the claim (v1 could not draw this) ----------------------------
fig, axs = plt.subplots(1, 2, figsize=(13.0, 4.3), sharey=True)
grid = np.linspace(0, 2*np.pi, 256)
for a, (k1, k2, ttl) in zip(axs, [("bobcat_only", "cat_only",
                                   "(a) Single-species sites  (C1: U2 = 0.0753, p = %.2f)"
                                   % pub["C1"]["p_perm"]),
                                  ("bobcat_both", "cat_both",
                                   "(b) Co-occurrence sites  (C2: U2 = 0.189, p = %.3f)"
                                   % pub["C2"]["p_perm"])]):
    for k in (k1, k2):
        x = D[k]; dd = density_fit(x, grid, get_bandwidth(x)/0.8)
        a.plot(grid, dd, color=ARMCOL[k], linewidth=2,
               label=f"{ARMLAB[k]}  (n = {len(x)})")
        a.plot(x, np.full(len(x), -0.004 - 0.006*(k == k2)), "|", color=ARMCOL[k],
               markersize=5, alpha=0.5)
    a.axvspan(0, np.pi/2, color="#DDDDDD", alpha=.5, zorder=0)
    a.axvspan(3*np.pi/2, 2*np.pi, color="#DDDDDD", alpha=.5, zorder=0)
    a.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    a.set_xticklabels(["midnight", "sunrise", "noon", "sunset", "midnight"])
    a.set_xlabel("solar time"); a.set_title(ttl, loc="left", fontsize=9.5)
    a.legend(fontsize=8, loc="upper center")
axs[0].set_ylabel("activity density (von Mises KDE, overlap's bandwidth)")
axs[0].text(0.02, 0.02, f"Dhat1 = {dh['solo']:.3f}\n95% [{dh['ci_solo_pct'][0]:.3f}, "
            f"{dh['ci_solo_pct'][1]:.3f}]", transform=axs[0].transAxes, fontsize=9,
            color="#333333")
axs[1].text(0.02, 0.02, f"Dhat1 = {dh['both']:.3f}\n95% [{dh['ci_both_pct'][0]:.3f}, "
            f"{dh['ci_both_pct'][1]:.3f}]", transform=axs[1].transAxes, fontsize=9,
            color="#333333")
fig.suptitle("F2 - The raw observations behind the diel-activity claims", fontsize=11,
             x=0.005, ha="left", y=1.03)
footer(fig, "v1 rendered this slot as a not-assessable panel because the detection times "
            "were unavailable. Every detection is now drawn as a rug tick; the curves are "
            "the von Mises kernel densities overlap uses, on the DST-aware solar-time basis "
            "the deposited R script computes. Shading marks night. The bracketed intervals "
            "are the bootstrap CIs the paper does not report.")
save(fig, "F2_raw_behind_claim.png")

# ---- F3 assumptions: how the null changes with the exchangeability assumption -----
fig, axs = plt.subplots(1, 2, figsize=(13.0, 4.2), gridspec_kw={"width_ratios": [1.3, 1]})
a = axs[0]
for nm, res, col in (("every detection exchangeable\n(the published analysis)", pub["C2"], OK["grey"]),
                     ("exchangeable only within site\n(respects the 7 sites)", pub["C2_site"], OK["green"]),
                     ("30-min independent events", pub["C2_ev30"], OK["orange"])):
    a.hist(res["null"], bins=70, histtype="step", linewidth=1.6, color=col, density=True,
           label=f"{nm}:  p = {res['p_perm']:.4f}")
a.axvline(pub["C2"]["observed"], color=OK["verm"], linewidth=2)
a.text(pub["C2"]["observed"]+0.006, a.get_ylim()[1]*0.75, "observed\nU2 = 0.189",
       color=OK["verm"], fontsize=8.5)
a.set_xlim(0, 0.5); a.set_xlabel("Watson's U2 under the null"); a.set_ylabel("density")
a.legend(fontsize=8, loc="upper right")
a.set_title("(a) C2's null depends on what you assume is exchangeable", loc="left", fontsize=9.5)

a = axs[1]
for k in KEYS:
    x = np.sort(D[k]); a.step(x, np.arange(1, len(x)+1)/len(x), color=ARMCOL[k],
                              where="post", label=f"{ARMLAB[k]} (n={len(x)})")
a.plot([0, 2*np.pi], [0, 1], "k:", linewidth=1, label="uniform")
a.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
a.set_xticklabels(["0", "sunrise", "noon", "sunset", "2pi"])
a.set_xlabel("solar time"); a.set_ylabel("empirical CDF")
a.legend(fontsize=7.5, loc="upper left")
a.set_title("(b) Shape of each arm against uniform (descriptive)", loc="left", fontsize=9.5)
fig.suptitle("F3 - Assumptions: exchangeability, not normality", fontsize=11, x=0.005,
             ha="left", y=1.03)
footer(fig, "Watson's U2 is distribution-free, so no normality diagnostic applies and panel "
            "(b) is descriptive only - it is not a test and nothing here is scored pass or "
            "fail. The assumption that carries C2 is exchangeability of detections, and "
            "panel (a) shows what it is worth: the same observed U2 = 0.189 sits at p = "
            f"{pub['C2']['p_perm']:.4f} against a null that shuffles every detection freely, "
            f"at p = {pub['C2_site']['p_perm']:.4f} against one that shuffles only within "
            f"site, and at p = {pub['C2_ev30']['p_perm']:.4f} once photographs of the same "
            "animal are collapsed into 30-minute events. The free-shuffle null is the one "
            "the published p-value corresponds to.")
save(fig, "F3_assumption_shape.png")

# ---- F4 claim ledger -------------------------------------------------------------
lf = LEDGER.copy()
lf["effect_reported"] = np.where(lf.effect_units.str.contains("pp |percentage"),
                                 lf.effect_reported, np.nan)
lf["effect_recomputed"] = lf.effect_pp
fig = rf.fig_claim_ledger(lf, path=None,
                          effect_label="difference, percentage points (pp-scaled claims only)")
fig.suptitle("F4 - Claim ledger: reported against recomputed, now from the data",
             fontsize=11, x=0.005, ha="left", y=1.02)
footer(fig, "Left: the claims whose effect is in percentage points (C4 developed cover, C6 "
            "nocturnal fraction against 50%, C7 nocturnal shift, C8 window shift); filled = "
            "recomputed with 95% CI, open = reported. C1 and C2 report Watson's U2 and C3 a "
            "Dhat difference, so they carry their own units and appear in the p-panel or the "
            "per-claim panels. Right: both p-values for every claim that has them; C1's "
            "reported value is the printed bound 0.20 and it sits far off the diagonal.",
       y=-0.09)
save(fig, "F4_claim_ledger.png")

# ---- F5 resolution ---------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13.0, 4.2), gridspec_kw={"width_ratios": [1.3, 1]})
a = axs[0]
a.hist(pub["C1"]["null"], bins=70, color="#DDDDDD", density=True, label="permutation null, C1")
uu = np.linspace(0.005, 0.55, 600)
surv = np.array([watson_asym_p(u) for u in uu])
a.plot(uu, -np.gradient(surv, uu), color="#333333", linewidth=1.5,
       label="asymptotic null (the table the paper's p came from)")
a.axvline(U2_CRIT_05, color=OK["verm"], linestyle="--", linewidth=1.2)
a.text(U2_CRIT_05+0.005, a.get_ylim()[1]*0.9, f"alpha = 0.05\nU2 = {U2_CRIT_05:.4f}",
       fontsize=8, color=OK["verm"])
for k, col, lab in (("C1", OK["blue"], "C1"), ("C2", OK["green"], "C2")):
    a.axvline(pub[k]["observed"], color=col, linewidth=2)
    a.text(pub[k]["observed"]+0.004, a.get_ylim()[1]*(0.55 if k == "C1" else 0.35),
           f"{lab}  U2 = {pub[k]['observed']:.4f}\np_perm = {pub[k]['p_perm']:.4f}\n"
           f"printed {'p < 0.20' if k == 'C1' else 'p < 0.05'}", color=col, fontsize=8)
a.set_xlim(0, 0.55); a.set_xlabel("Watson's U2"); a.set_ylabel("density under the null")
a.legend(fontsize=8, loc="center right")
a.set_title("(a) Permutation null from the data against the asymptotic table",
            loc="left", fontsize=9.5)

a = axs[1]
rows = [("C1", 0.20, pub["C1"]["p_perm"], OK["blue"]), ("C2", 0.05, pub["C2"]["p_perm"], OK["green"])]
for i, (lab, pr, pc, col) in enumerate(rows):
    a.plot([pr], [i], "o", mfc="white", mec=col, ms=11, mew=1.6)
    a.plot([pc], [i], "o", color=col, ms=9)
    a.plot([pr, pc], [i, i], ":", color=col, linewidth=1.2)
a.axvline(ALPHA, color=OK["verm"], linestyle="--", linewidth=1.2)
a.axvline(1/(B_PERM+1), color="#333333", linestyle=":", linewidth=1.2)
a.text(1/(B_PERM+1)*1.2, 1.55, f"Monte Carlo floor\n1/(B+1) = {1/(B_PERM+1):.1e}",
       fontsize=7.5, color="#333333")
a.text(ALPHA*1.15, -0.42, "alpha = 0.05", fontsize=8, color=OK["verm"])
a.set_xscale("log"); a.set_xlim(2e-5, 1.5); a.set_ylim(-0.6, 1.9)
a.set_yticks([0, 1]); a.set_yticklabels(["C1  49 vs 696", "C2  32 vs 99"], fontsize=8.5)
a.set_xlabel("p-value (log scale);  open = printed bound, filled = permutation")
a.set_title("(b) Printed bound against the p the data give", loc="left", fontsize=9.5)
fig.suptitle("F5 - Resolution: what each p could have been", fontsize=11, x=0.005,
             ha="left", y=1.03)
footer(fig, f"The permutation null computed from the authors' own detection times tracks the "
            f"asymptotic table closely (empirical 95th percentile "
            f"{pub['C1']['null_q95']:.4f} against the tabulated {U2_CRIT_05:.4f}), which "
            f"confirms that reading p off that table was the right procedure - and confirms "
            f"v1's asymptotic re-derivations (0.447 and 0.048) against the permutation values "
            f"({pub['C1']['p_perm']:.3f} and {pub['C2']['p_perm']:.4f}). No combinatorial "
            f"floor binds at these sample sizes; the Monte Carlo floor of this review is "
            f"1/(B+1) = {1/(B_PERM+1):.1e} with B = {B_PERM:,}, and the finest value the "
            f"paper's own critical-value table could have printed is p < 0.001.")
save(fig, "F5_null_distribution.png")

# ---- F6 multiplicity -------------------------------------------------------------
fig = rf.fig_multiplicity(p_family, alpha=ALPHA, path=None, labels=labels)
fig.suptitle("F6 - Multiplicity across the two inferential tests the paper reports",
             fontsize=11, x=0.005, ha="left", y=1.03)
footer(fig, "m = 2. The histogram is drawn because the spine is fixed, not because it can be "
            "read: a p-value distribution says nothing at m = 2, nothing about any single "
            "result, and nothing about why any result looks as it does. The right panel is "
            f"the part that counts - the one uncorrected hit (C2, p = {pub['C2']['p_perm']:.4f}) "
            "survives neither Bonferroni nor BH. Whether two pre-specified comparisons on "
            "disjoint site sets form a family that should be corrected at all is the authors' "
            "argument to make; both numbers are given so a reader can take either view. The "
            "reviewer-side tests added in this review are not part of the paper's family and "
            "are not counted here.", y=-0.07)
save(fig, "F6_multiplicity.png")

# ---- panel C1 --------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13.0, 4.0), gridspec_kw={"width_ratios": [1.2, 1]})
a = axs[0]
items = [("omnibus Watson U2,\npublished clock basis", pub["C1"]["p_perm"], OK["blue"]),
         ("omnibus Watson U2,\nsolar basis (the Methods)", sol["C1"]["p_perm"], OK["sky"]),
         ("night/day fraction\n33/49 vs 594/696 (Fisher)",
          R["prop_tests"]["C1_binary_bobcat_vs_cat_at_solo"], OK["orange"])]
for i, (lab, p_, col) in enumerate(items):
    a.plot([p_], [i], "o", color=col, ms=11)
    a.text(p_*1.45, i+0.17, f"p = {p_:.4g}", fontsize=8.5, color=col)
a.axvline(ALPHA, color=OK["verm"], linestyle="--", linewidth=1.2)
a.text(ALPHA*1.15, 2.45, "alpha = 0.05", fontsize=8, color=OK["verm"])
a.set_xscale("log"); a.set_xlim(5e-4, 3); a.set_ylim(-0.5, 2.8)
a.set_yticks(range(3)); a.set_yticklabels([i[0] for i in items], fontsize=8.5)
a.set_xlabel("p-value (log scale)")
a.set_title("(a) Three summaries of the same 745 detections", loc="left", fontsize=9.5)
a = axs[1]
for k, off in (("bobcat_only", -0.15), ("cat_only", 0.15)):
    n_ = nc[k]
    a.errorbar(n_["pct"], off, xerr=[[n_["pct"]-n_["wilson"][0]], [n_["wilson"][1]-n_["pct"]]],
               fmt="o", color=ARMCOL[k], ms=11, capsize=5, linewidth=2)
    a.text(n_["pct"], off+0.12, f"{n_['k_night']}/{n_['n']} = {n_['pct']:.1f}%",
           ha="center", fontsize=9, color=ARMCOL[k])
a.axvline(50, color="#444444", linestyle=":", linewidth=1)
a.set_yticks([-0.15, 0.15]); a.set_yticklabels(["bobcat", "domestic cat"], fontsize=9)
a.set_ylim(-0.5, 0.5); a.set_xlim(40, 100)
a.set_xlabel("detections at night (%), Wilson 95% interval")
a.set_title("(b) Where the two summaries part company", loc="left", fontsize=9.5)
fig.suptitle("C1 - 'The diel activity patterns were similar where only one species occurred' "
             "-- assumption-dependent", fontsize=10.5, x=0.005, ha="left", y=1.05)
footer(fig, f"The statistic reproduces exactly from the deposit and the permutation p "
            f"({pub['C1']['p_perm']:.3f}) confirms v1's asymptotic 0.447, so the printed "
            f"'p < 0.20' is the wrong side of the 0.20 critical value of 0.1166. The omnibus "
            f"test spreads its power around the whole circle; the night/day contrast spends "
            f"all of it on one split, where the species differ by 18.0 pp with Fisher p = "
            f"{R['prop_tests']['C1_binary_bobcat_vs_cat_at_solo']:.4f}. Neither route is "
            f"wrong. Reporting only the first and describing the patterns as 'similar' is the "
            f"finding, and an accepted null still needs an equivalence bound the paper does "
            f"not give.")
save(fig, "panel_C1.png")

# ---- panel C2 --------------------------------------------------------------------
fig, a = plt.subplots(figsize=(8.6, 4.4))
bars = [("published basis,\nfree shuffle", pub["C2"]["p_perm"], OK["grey"]),
        ("published basis,\nwithin-site", pub["C2_site"]["p_perm"], OK["green"]),
        ("published basis,\n30-min events", pub["C2_ev30"]["p_perm"], OK["orange"]),
        ("solar basis,\nfree shuffle", sol["C2"]["p_perm"], OK["sky"]),
        ("solar basis,\nwithin-site", sol["C2_site"]["p_perm"], OK["blue"])]
xs = np.arange(len(bars))
bb = a.bar(xs, [b[1] for b in bars], color=[b[2] for b in bars], edgecolor="white", width=.62)
for x, (_, p_, _c) in zip(xs, bars):
    a.text(x, p_+0.0025, f"{p_:.4f}", ha="center", fontsize=9)
a.axhline(ALPHA, color=OK["verm"], linestyle="--", linewidth=1.4)
a.text(len(bars)-0.4, ALPHA+0.0022, "alpha = 0.05", fontsize=8.5, color=OK["verm"], ha="right")
a.set_xticks(xs); a.set_xticklabels([b[0] for b in bars], fontsize=8)
a.set_ylabel("permutation p for Watson's U2"); a.set_ylim(0, 0.125)
a.set_title("C2 - 'Domestic cat diel activity differed significantly from bobcat activity'\n"
            "32 bobcat vs 99 cat detections at 7 shared sites -- assumption-dependent",
            loc="left", fontsize=10)
footer(fig, "Five defensible analyses of the same claim, three above alpha and two below. The "
            "published U2 = 0.189 reproduces exactly, so this is not an arithmetic dispute: "
            "it is that the result sits within the range over which the choice of time base "
            "and null changes the decision. v1 reached the same place from the reported "
            "numbers alone, predicting that a design effect of only 1.011 would cross alpha; "
            "with the data the within-site null does exactly that. Note the two bases are not "
            "equally defensible - the Methods describe the solar one, and the published "
            "statistic was computed on the other.", y=-0.09)
save(fig, "panel_C2.png")

# ---- panel C3 --------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13.0, 4.2), gridspec_kw={"width_ratios": [1, 1.1]})
a = axs[0]
for i, (k, lab, col, n_) in enumerate([("solo", "single-species sites", OK["blue"], "49 vs 696"),
                                       ("both", "co-occurrence sites", OK["green"], "32 vs 99")]):
    est, lo, hi = dh[k], dh[f"ci_{k}_pct"][0], dh[f"ci_{k}_pct"][1]
    a.errorbar(est, i, xerr=[[est-lo], [hi-est]], fmt="o", color=col, ms=13, capsize=6, lw=2.4)
    a.text(est, i+0.17, f"Dhat1 = {est:.3f}   95% [{lo:.3f}, {hi:.3f}]", ha="center",
           fontsize=9.5, color=col)
    a.text(est, i-0.22, f"n = {n_} detections", ha="center", fontsize=8.5, color=col)
a.set_yticks([0, 1]); a.set_yticklabels(["single-species", "co-occurrence"], fontsize=9)
a.set_ylim(-0.6, 1.6); a.set_xlim(0.45, 1.0)
a.set_xlabel("coefficient of overlapping, Dhat1 (bootstrap 95% CI)")
a.set_title("(a) The intervals the paper omits -- and they overlap", loc="left", fontsize=9.5)

a = axs[1]
a.hist(bs_diff, bins=60, color="#DDDDDD", edgecolor="white")
a.axvline(0, color=OK["verm"], linewidth=2)
a.axvline(dh["difference"], color=OK["green"], linewidth=2)
for x, t, col, yf in ((0, "no difference", OK["verm"], 0.92),
                      (dh["difference"], f"observed\n{dh['difference']:.3f}", OK["green"], 0.92)):
    a.text(x+0.006, a.get_ylim()[1]*yf, t, color=col, fontsize=8.5, va="top")
lo, hi = dh["ci_diff_pct"]
a.axvspan(lo, hi, color=OK["green"], alpha=0.12)
a.set_xlabel("Dhat1(single-species) - Dhat1(co-occurrence)")
a.set_ylabel(f"bootstrap replicates (of {NB_BOOT:,})")
a.set_title(f"(b) 95% CI [{lo:.3f}, {hi:.3f}] includes zero", loc="left", fontsize=9.5)
fig.suptitle("C3 - 'Activity patterns overlapped much less where they co-occurred' -- "
             "assumption-dependent", fontsize=10.5, x=0.005, ha="left", y=1.04)
footer(fig, f"v1 could not assess this and named the bootstrap as what would resolve it; this "
            f"is that bootstrap, using overlap's own Dhat1 and smoothed resampler. Both point "
            f"estimates reproduce (0.852, 0.647). The difference is {dh['difference']:.3f} "
            f"with a 95% interval of [{lo:.3f}, {hi:.3f}] that includes zero, and at site "
            f"level [{dh['ci_diff_cluster'][0]:.3f}, {dh['ci_diff_cluster'][1]:.3f}]. So "
            f"'much less' is not established by Dhat. The underlying biology is nonetheless "
            f"on firm ground by a different route the paper never takes: cat activity differs "
            f"between the two site classes with U2 = {sol['cat_shift']['observed']:.3f}, "
            f"p <= {sol['cat_shift']['p_floor']:.0e}, while bobcat activity does not "
            f"(U2 = {sol['bobcat_shift']['observed']:.3f}, p = "
            f"{sol['bobcat_shift']['p_perm']:.3f}).")
save(fig, "panel_C3.png")

# ---- panel C6/C7 -----------------------------------------------------------------
fig, a = plt.subplots(figsize=(9.2, 4.6))
rows = [("bobcat at\nbobcat-only sites", "bobcat_only"), ("cat at\ncat-only sites", "cat_only"),
        ("bobcat at\nco-occurrence sites", "bobcat_both"), ("cat at\nco-occurrence sites", "cat_both")]
for i, (lab, k) in enumerate(rows):
    n_ = nc[k]
    a.errorbar(n_["pct"], i, xerr=[[n_["pct"]-n_["wilson"][0]], [n_["wilson"][1]-n_["pct"]]],
               fmt="o", color=ARMCOL[k], ms=11, capsize=5, lw=2.2)
    a.text(n_["wilson"][1]+1.5, i, f"{n_['k_night']}/{n_['n']} = {n_['pct']:.1f}%   "
           f"[{n_['wilson'][0]:.1f}, {n_['wilson'][1]:.1f}]", va="center", fontsize=8.5,
           color=ARMCOL[k])
    if k == "bobcat_both":
        a.text(n_["pct"], i-0.33, "the number the paper never reported", ha="center",
               fontsize=8, color=OK["verm"], style="italic")
a.axvline(50, color="#444444", linestyle=":", linewidth=1.2)
a.text(50.6, -0.45, "50%", fontsize=8, color="#444444")
a.set_yticks(range(4)); a.set_yticklabels([r[0] for r in rows], fontsize=8.5)
a.set_ylim(-0.6, 3.6); a.set_xlim(35, 150); a.invert_yaxis()
a.annotate("", xy=(nc["cat_both"]["pct"], 3), xytext=(nc["cat_only"]["pct"], 1),
           arrowprops=dict(arrowstyle="->", color=OK["verm"], linewidth=1.8))
a.text(86, 1.75, f"C7: {R['C7_cluster_bootstrap']['gap_pp']:.1f} pp shift\n"
       f"Fisher p = {R['prop_tests']['C7_cat_only_vs_cat_both']:.1e}\nsite-level 95% "
       f"[{R['C7_cluster_bootstrap']['ci_pp'][0]:.1f}, "
       f"{R['C7_cluster_bootstrap']['ci_pp'][1]:.1f}] pp", fontsize=8.5, color=OK["verm"])
a.set_xlabel("detections at night (%), Wilson 95% interval")
a.set_title("C6 / C7 - Nocturnal fraction, exact counts from the deposit", loc="left",
            fontsize=10)
footer(fig, f"v1 reconstructed these counts from the rounded percentages and had to report "
            f"ranges; the data give 33/49, 594/696, 51/99 exactly - v1's unique "
            f"reconstructions (33) and range (589-595) were right, its 50/99 was one out. The "
            f"new row is the third: bobcats at co-occurrence sites, "
            f"{nc['bobcat_both']['k_night']}/{nc['bobcat_both']['n']} = "
            f"{nc['bobcat_both']['pct']:.1f}% nocturnal with a Wilson interval that includes "
            f"50% (binomial p = {nc['bobcat_both']['binom_p_vs_half']:.3f}). That is half of "
            f"C6's 'wherever they were present', and it is not established. C7's shift is "
            f"large and highly significant at detection level; the site-level interval is "
            f"wide because one site supplies 234 of 795 cat detections, and because cat-only "
            f"and co-occurrence sites are different places.", y=-0.08)
save(fig, "panel_C6_C7_nocturnal.png")

# ---- panel C8 --------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13.0, 4.0), gridspec_kw={"width_ratios": [1.5, 1]})
a = axs[0]
bins = np.arange(0, 25, 1)
for k in ("cat_only", "cat_both"):
    hh = det.loc[det.arm == k, "min_local_dst"].values/60.0
    a.hist(hh, bins=bins, density=True, histtype="step", linewidth=2, color=ARMCOL[k],
           label=f"{ARMLAB[k]} (n={len(hh)})")
for lo_, hi_ in ((3, 10), (16, 21)):
    a.axvspan(lo_, hi_, color=OK["orange"], alpha=0.13, zorder=0)
a.text(6.5, a.get_ylim()[1]*0.93, "03:00-10:00", ha="center", fontsize=8, color=OK["verm"])
a.text(18.5, a.get_ylim()[1]*0.93, "16:00-21:00", ha="center", fontsize=8, color=OK["verm"])
a.set_xlabel("local clock hour (DST-aware)"); a.set_ylabel("density of detections")
a.set_xticks(range(0, 25, 3)); a.legend(fontsize=8, loc="upper left")
a.set_title("(a) The windows the paper names", loc="left", fontsize=9.5)
a = axs[1]
cats = ["in the named\nwindows", "mid-day or\nmid-night"]
co = [c8["named_windows"]["cat_only_pct"], c8["mid_day_or_mid_night"]["cat_only_pct"]]
cb_ = [c8["named_windows"]["cat_both_pct"], c8["mid_day_or_mid_night"]["cat_both_pct"]]
x = np.arange(2); w = 0.36
a.bar(x-w/2, co, w, color=ARMCOL["cat_only"], label="cat-only sites", edgecolor="white")
a.bar(x+w/2, cb_, w, color=ARMCOL["cat_both"], label="co-occurrence sites", edgecolor="white")
for xx, v1_, v2_ in zip(x, co, cb_):
    a.text(xx-w/2, v1_+1, f"{v1_:.1f}%", ha="center", fontsize=8.5)
    a.text(xx+w/2, v2_+1, f"{v2_:.1f}%", ha="center", fontsize=8.5)
a.text(0, 73, f"Fisher p = {c8['named_windows']['fisher_p']:.1e}", ha="center", fontsize=8.5,
       color=OK["verm"])
a.text(1, 73, f"Fisher p = {c8['mid_day_or_mid_night']['fisher_p']:.3f}", ha="center",
       fontsize=8.5, color="#444444")
a.set_xticks(x); a.set_xticklabels(cats, fontsize=8.5); a.set_ylim(0, 82)
a.set_ylabel("% of cat detections"); a.legend(fontsize=8, loc="upper right")
a.set_title("(b) First half holds, second half does not", loc="left", fontsize=9.5)
fig.suptitle("C8 - 'Domestic cats shifted toward 03:00-10:00 and 16:00-21:00 with very little "
             "activity in the middle of the day or night'", fontsize=10.5, x=0.005, ha="left",
             y=1.04)
footer(fig, f"v1 could not assess this claim at all. The shift into the two named windows is "
            f"real ({c8['named_windows']['cat_both_pct']:.1f}% against "
            f"{c8['named_windows']['cat_only_pct']:.1f}%, Fisher p = "
            f"{c8['named_windows']['fisher_p']:.1e}). The second half of the sentence is not "
            f"supported: {c8['mid_day_or_mid_night']['cat_both_pct']:.1f}% of co-occurrence "
            f"cat detections still fall in the middle of the day or night, which is not a "
            f"reduction the data separate from the "
            f"{c8['mid_day_or_mid_night']['cat_only_pct']:.1f}% at cat-only sites "
            f"(p = {c8['mid_day_or_mid_night']['fisher_p']:.3f}). Clock time is used here "
            f"because the claim is stated in clock hours.")
save(fig, "panel_C8.png")

R["figures"] = {k: os.path.relpath(v, HERE).replace("\\", "/") for k, v in FIGS.items()}

# ==================================================================================
# J.  The report.  Section order follows REPORT_FORMAT.md exactly.
# ==================================================================================
def fp(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)): return "--"
    return f"{x:.2e}" if x < 1e-3 else f"{x:.4f}"

rows = ["| claim | verdict | v1 | effect (paper's units) | recomputed | p rep. | p recomp. |",
        "|---|---|---|---|---|---|---|"]
for r in LEDGER.itertuples():
    er = "--" if not np.isfinite(r.effect_reported) else f"{r.effect_reported:g}"
    ec = ("--" if not np.isfinite(r.effect_recomputed) else
          (f"{r.effect_recomputed:.3g} [{r.ci_low:.3g}, {r.ci_high:.3g}]"
           if np.isfinite(r.ci_low) else f"{r.effect_recomputed:.4g}"))
    rows.append(f"| **{r.claim_id}** {r.claim_short} | `{r.verdict}` | {r.v1_status} | "
                f"{er} {r.effect_units} | {ec} | {fp(r.p_reported)} | {fp(r.p_recomputed)} |")
ledger_md = "\n".join(rows)
vc_md = "\n".join(f"| `{k}` | {v} |" for k, v in LEDGER.verdict.value_counts().items())
cl = R["clustering"]

REPORT = f"""# Statistical review (v2 -- data-based)

**Alavi Naini L., Linde N., Malik S., Syed Y. & Cheek A.O. (2026)** "Diel Activity
Patterns of Bobcats and Domestic Cats in the Houston Metropolitan Area."
*microPublication Biology*. DOI: [10.17912/micropub.biology.001911](https://doi.org/10.17912/micropub.biology.001911)

Reviewed with the `statistical-review` skill, six stages in order.

> **This version supersedes v1.** v1 was a text-only review: the Dryad deposit the paper
> names could not be retrieved from the review sandbox. The deposit has since been
> supplied and verified against the Dryad API manifest, and every claim below is now
> re-derived from the detection-level data and the authors' own analysis scripts.

---

## 1. Header -- what was reviewed

| input | condition |
|---|---|
| paper | full text from the PMC XML, 20,192 characters; text layer usable |
| raw data | **obtained.** `Bobcat_Data.csv` (81 rows), `Domestic_Cat_Data.csv` (795 rows), both 22 columns, all rows Status "confirmed correct" |
| analysis code | `Watsons U2 and Polar Plot.m` (the published U2 values) and `Dhat R Script MicroPub.R` (the published Dhat values), both deposited |
| data availability | the paper **does** carry a statement: "Animal detection data, MatLab script and R script are available at datadryad.org, DOI: 10.5061/dryad.s4mw6m9jg" |
| provenance | v1 recorded the deposit as existing but unreachable (HTTP 401 on the file endpoint, then loss of network egress). The user downloaded it manually; both files match the Dryad API manifest byte-for-byte. **The v1 data gap was a retrieval failure in the review environment, not a reporting failure by the authors** -- the paper's data availability statement is accurate and the deposit is complete enough to reproduce the headline numbers. |

**Reproduction is exact.** Every headline statistic re-derives from the deposit:

| statistic | paper | this review |
|---|---|---|
| Watson's U2, single-species sites | 0.0753 | **{pub['C1']['observed']:.6f}** |
| Watson's U2, co-occurrence sites | 0.189 | **{pub['C2']['observed']:.6f}** |
| Dhat, single-species sites | 0.852 | **{dh['solo']:.4f}** |
| Dhat, co-occurrence sites | 0.647 (deposit plot: 0.646) | **{dh['both']:.4f}** |
| nocturnal, bobcat at bobcat-only | 67% | **33/49 = {nc['bobcat_only']['pct']:.1f}%** |
| nocturnal, cat at cat-only | 85% | **594/696 = {nc['cat_only']['pct']:.1f}%** |
| nocturnal, cat at co-occurrence | 51% | **51/99 = {nc['cat_both']['pct']:.1f}%** |
| detections / sites | 81, 795 / 7, 19, 7 | 81, 795 / 7, 19, 7 |

**Coverage.** 11 claims, all with a verdict. Three were `not-assessable` in v1 and are now
assessable (C3, C8, C9); two remain `not-assessable` (C5, C10) because the deposit does not
contain the trap-days, camera coordinates or land-cover measurements they need.

---

## 2. Summary of findings

| verdict | n |
|---|---|
{vc_md}

**Finding 1 -- the two analyses run on two different clocks (stage 3).** The deposited
MATLAB reads column 9, `Local Time Value (ABS)`, a day fraction built on a **fixed UTC-6
offset that does not follow daylight saving**: {R['dst_offset_audit']['n_shifted_by_60min']}
of {R['dst_offset_audit']['n_detections']} detections are an hour behind true local time
there. The deposited R recomputes local time DST-aware from the UTC timestamp. Worse, the
Methods say timestamps were "converted to solar time using the sunTime function" and that
"Watson's U2 test was used to compare diel activity" -- but the solar conversion exists
only in the R script, which computes Dhat. **The published U2 values were computed on raw,
non-DST clock time.** Recomputing U2 on the basis the Methods describe gives
{sol['C1']['observed']:.4f} (p = {sol['C1']['p_perm']:.3f}) and {sol['C2']['observed']:.4f}
(p = {sol['C2']['p_perm']:.4f}) -- different numbers from the published 0.0753 and 0.189.

**Finding 2 -- C2's significance does not survive the design (stages 2 and 5).** U2 = 0.189
reproduces exactly, but its p depends on what is assumed exchangeable: {pub['C2']['p_perm']:.4f}
shuffling every detection freely (the published analysis), **{pub['C2_site']['p_perm']:.4f}
shuffling only within site**, and **{pub['C2_ev30']['p_perm']:.4f}** after a standard 30-minute
independent-event filter. v1 predicted this from the reported numbers alone -- it calculated
that a design effect of only 1.011 would cross alpha -- and the data confirm it. Two of the
seven shared sites supply 91 of the 99 co-occurrence cat detections.

**Finding 3 -- the central claim's interval includes zero, but a better test does not
(stage 5).** The paper's temporal-partitioning claim is the gap between Dhat 0.852 and
0.647. The deposited R script computes bootstrap CIs and the paper reports none; they are
[{dh['ci_solo_pct'][0]:.3f}, {dh['ci_solo_pct'][1]:.3f}] and
[{dh['ci_both_pct'][0]:.3f}, {dh['ci_both_pct'][1]:.3f}], and they overlap. The difference
is {dh['difference']:.3f} with a 95% interval of
[{dh['ci_diff_pct'][0]:.3f}, {dh['ci_diff_pct'][1]:.3f}] -- **including zero**. The biology
is nonetheless on firm ground by a route the paper never takes: cat activity differs
between the two site classes (U2 = {sol['cat_shift']['observed']:.3f},
p <= {sol['cat_shift']['p_floor']:.0e}) while bobcat activity does not
(U2 = {sol['bobcat_shift']['observed']:.3f}, p = {sol['bobcat_shift']['p_perm']:.3f}).
**The paper should make that claim instead.**

---

## 3. Claim ledger -- F4

![F4](figures/F4_claim_ledger.png)

{ledger_md}

Full `reasons`, both p-values, the v1 verdict and the design floor for every claim are in
`claim_ledger.csv`; every reported statistic is checked against the data in `statcheck.csv`.

### What the data confirmed, and what they overturned

**Confirmed from v1:**

- C1's printed "p < 0.20" points the wrong way. v1 derived p = 0.447 asymptotically; the
  permutation p from the authors' own data is {pub['C1']['p_perm']:.3f}.
- C2 is fragile to clustering. v1 computed that a design effect of 1.011 would cross alpha;
  the within-site permutation gives {pub['C2_site']['p_perm']:.4f}.
- v1's GRIM reconstruction of the nocturnal counts: 33/49 recovered uniquely and correctly;
  594 falls inside the reconstructed range 589-595. Its 50/99 was one out -- the value is 51/99.
- C4 and C5 are unchanged: the deposit contains no land cover, no trap-days and no camera
  coordinates.

**Overturned or newly resolved:**

- **C2 downgraded** from `supported-with-caveat` to `assumption-dependent`: with the data,
  three of five defensible analyses put it above alpha.
- **C3 now assessable** and downgraded in substance: the bootstrap interval on the Dhat
  difference includes zero.
- **C8 and C9 now assessable.** C8's first half holds and its second half does not. C9's
  collar annotation exists in the deposit: {R['C9_collar']['n_collar']} collared against
  {R['C9_collar']['n_no_collar']} not, on {R['C9_collar']['n_annotated']} annotated of 795.
- **C6's missing half is now supplied and does not support the claim.** Bobcats at
  co-occurrence sites are {nc['bobcat_both']['k_night']}/{nc['bobcat_both']['n']} =
  {nc['bobcat_both']['pct']:.1f}% nocturnal, binomial p = {nc['bobcat_both']['binom_p_vs_half']:.3f},
  with a Wilson interval that includes 50%.

---

## 4. Design and sampling -- F1, F2

![F1](figures/F1_design_map.png)

**Stage 1.** The population is free-ranging bobcats and unowned domestic cats across the
22,891 km2 Houston metropolitan area; the sample is 33 camera sites in parks, cemeteries,
golf courses and natural areas along two transects. Whether that convenience sample
represents the metropolitan area is a design question the methods do not address and the
data cannot settle.

**Stage 2 -- pseudoreplication, now measured rather than bounded.** v1 could only show a
sensitivity curve over an unknown ICC. The site labels in the deposit make it concrete:
cat detections are distributed across 26 sites with mean {cl['detections_per_site_cat']['mean']:.1f}
and s.d. {cl['detections_per_site_cat']['std']:.1f}, and **one site contributes
{cl['top_cat_site_n']} of the 795 cat detections ({100*cl['top_cat_site_share_of_all_cat']:.0f}%)**.
Bobcat detections are sparser and flatter (mean {cl['detections_per_site_bobcat']['mean']:.1f},
max {cl['detections_per_site_bobcat']['max']:.0f} across 14 sites). Within the co-occurrence
comparison, two of seven sites supply 91 of 99 cat detections and 22 of 32 bobcat detections.

No independence filter is applied anywhere in the paper or the deposited code: consecutive
photographs of the same animal enter as separate detections. The largest single burst is 14
cat photographs spanning 6.8 minutes at one site. A 30-minute quiet interval reduces 876
detections to {cl['events_by_quiet_interval'][30]} events -- an 11% reduction overall, but
enough to move C2 (section 6).

The totals all check out against the rows (C11, `supported`): 81 and 795 detections, site
classes of 7/19/7 derived by the deposited script's own rule, arms of 49/696/32/99.

![F2](figures/F2_raw_behind_claim.png)

---

## 5. Method and assumptions -- F3

![F3](figures/F3_assumption_shape.png)

**Stage 3 -- test selection.** Watson's U2 is an appropriate omnibus two-sample test for
circular data and the deposited implementation is the standard two-sample formula, correctly
rotation-invariant. Three observations:

1. **It is an omnibus test**, so it is least efficient against an alternative concentrated in
   one contrast. At single-species sites the night/day contrast is exactly such an
   alternative, and the two routes disagree (C1).
2. **The time base is wrong relative to the Methods** (Finding 1). The statistic is invariant
   to the minutes-versus-radians scaling -- that part is harmless, since the transform is
   monotone and U2 depends only on sort order -- but it is not invariant to the DST offset
   or to the solar-time transformation.
3. **Ties.** The implementation uses MATLAB's stable sort, so cross-group ties resolve in
   favour of the first-listed group. There is exactly one cross-group tie in the
   single-species comparison and none in the co-occurrence comparison, so tie handling does
   not affect either published value.

**Stage 4 -- assumptions.** Watson's U2 is distribution-free: no normality assumption is at
stake, no shape diagnostic applies, and nothing here is scored pass or fail. The assumption
that carries the inference is **exchangeability of detections**, and F3(a) prices it
directly -- the same observed statistic maps to {pub['C2']['p_perm']:.4f},
{pub['C2_site']['p_perm']:.4f} or {pub['C2_ev30']['p_perm']:.4f} depending on which
exchangeability the null encodes.

The three routes in the skill's usual sense (parametric, rank, permutation) do not transfer
to a circular omnibus statistic with no parametric analogue in use here. The analogous set
is reported instead and in full: the asymptotic table the authors used, the free-shuffle
permutation, the within-site permutation, the event-filtered version, and both time bases.
They are listed together in panel C2 and none was selected after the fact.

On error bars: the land-cover "+/- " values are labelled 1 s.d. for two of three site
classes. Read as s.e.m. they would imply site-level s.d. of 50 and 61 percentage points on
a 0-100 bounded scale, which the bound rules out, so the s.d. label is the consistent
reading. The site values themselves are not in the deposit, so this remains an argument
from consistency, not a verification.

---

## 6. Resolution and evidence -- F5

![F5](figures/F5_null_distribution.png)

The permutation null computed from the authors' own detection times tracks the asymptotic
table closely (empirical 95th percentile {pub['C1']['null_q95']:.4f} against the tabulated
{U2_CRIT_05:.4f}), which confirms both that reading p off that table was the right procedure
and that v1's asymptotic re-derivations were correct.

| | C1 single-species | C2 co-occurrence |
|---|---|---|
| U2 (reproduced) | {pub['C1']['observed']:.6f} | {pub['C2']['observed']:.6f} |
| printed | p < 0.20 | p < 0.05 |
| permutation p, published basis | **{pub['C1']['p_perm']:.4f}** | **{pub['C2']['p_perm']:.4f}** |
| within-site permutation | not defined (see below) | **{pub['C2_site']['p_perm']:.4f}** |
| 30-min independent events | {pub['C1_ev30']['p_perm']:.4f} | **{pub['C2_ev30']['p_perm']:.4f}** |
| permutation p, solar basis | {sol['C1']['p_perm']:.4f} | {sol['C2']['p_perm']:.4f} |
| solar basis, within-site | not defined | {sol['C2_site']['p_perm']:.4f} |
| smallest p this design could give | ~1e-77 (combinatorial) | ~1e-30 (combinatorial) |
| smallest p this review could give | {1/(B_PERM+1):.1e} (B = {B_PERM:,}) | {1/(B_PERM+1):.1e} |
| smallest p the paper's table could print | p < 0.001 | p < 0.001 |

**No combinatorial floor binds** at these sample sizes. The binding limits are the Monte
Carlo resolution of this review and, for the paper, the granularity of the printed bracket.

**A within-site null does not exist for C1, and that is a structural finding.** At
single-species sites each site has, by construction, only one species, so species and site
are perfectly confounded: no permutation that respects site can test the species contrast.
Any comparison between single-species classes is a comparison between different places. The
same applies to C7. This is a limit of the design, not of the analysis.

---

## 7. Multiplicity -- F6

![F6](figures/F6_multiplicity.png)

The paper reports two inferential tests. The single uncorrected hit (C2,
p = {pub['C2']['p_perm']:.4f}) survives neither Bonferroni nor Benjamini-Hochberg (both
{mp.bonferroni(p_family, ALPHA).to_dict('records')[1]['p_adjusted']:.4f}). Whether two
pre-specified comparisons on disjoint site sets form a family that should be corrected at
all is a defensible argument for the authors; both numbers are given.

Untested comparisons described in the text remain untested: six land-cover contrasts, three
nocturnal percentages, the Dhat difference, and the hourly-window claim. This review supplies
tests for several of them, but those are the reviewer's and are not counted in the paper's
family.

**On the p-value distribution:** with m = 2 it cannot be read. The just-below-threshold
diagnostic is a property of a collection of p-values; it says nothing about any individual
result, and nothing about why any result looks the way it does. No inference about the
authors' conduct is drawn here or anywhere in this report.

---

## 8. Per-claim supplement

### C1 -- "the diel activity patterns were similar where only one species occurred" -- `assumption-dependent`
![C1](figures/panel_C1.png)

### C2 -- "domestic cat diel activity differed significantly from bobcat activity" -- `assumption-dependent`
![C2](figures/panel_C2.png)

### C3 -- "activity patterns overlapped much less where they co-occurred" -- `assumption-dependent`
![C3](figures/panel_C3.png)

### C6 / C7 -- nocturnal fractions -- `supported-with-caveat`
![C7](figures/panel_C6_C7_nocturnal.png)

### C8 -- the named activity windows -- `supported-with-caveat`
![C8](figures/panel_C8.png)

---

## 9. What was not assessed

{ac.uncheckable_assumptions_note()}

**Specific to this review:**

- **Independence is now measured but not resolved.** The site labels quantify clustering and
  the timestamps quantify photographic bursts, but the deposit contains no individual
  identities. 876 detections come from an unknown number of animals; the 1 km camera spacing
  separates cameras, not successive visits by one animal. Every p-value here treats
  detections, or 30-minute events, as the unit; none treats individuals as the unit, because
  nothing in the deposit can.
- **Confounding of species with site cannot be removed** for C1, C7 or C8: single-species and
  co-occurrence sites are different places, so "bobcat present" is confounded with everything
  else that differs between them. Only C2 compares species within the same sites.
- **Claims the deposit still cannot reach:** C5 (no per-site trap-days; no camera
  coordinates -- every row carries the same Latitude/Longitude, the downtown reference point,
  so Fig 1B is not reproducible either) and C10 (coexistence is not measured and no site
  changes predator composition within the study).
- **C4 still rests on published summaries**, not site-level values; its Welch tests apply
  normal theory to percentages bounded at 0-100 at n = 7 and support a direction, not a
  precise p.
- **The night boundary is not defined in the Methods** and no definition tested here
  reproduces all three printed percentages exactly. Sunrise/sunset on the DST-aware solar
  basis reproduces 67% and 85% and gives 51.5% where the paper prints 51%. The nocturnal
  figures in this review use that definition.
- **`overlap`'s compiled kernel could not be loaded in this environment**, so `Dhat1`, its
  bandwidth rule and its smoothed bootstrap resampler were transcribed from the package's own
  R source and verified by reproducing the deposit's published values (0.852 and 0.646/0.647).
  `suntools::sunriset` was likewise reimplemented (NOAA) and validated against published
  Houston sunrise and sunset times to within about two minutes.
- **Nothing was subsampled.** All 876 detections were used; the 30-minute event filter is
  reported as a sensitivity alongside the full-data result, never in place of it.

---

## 10. Reproduction

| item | value |
|---|---|
| script | `review_alavinaini_2026.py` -- regenerates every figure, CSV and this report |
| data | Dryad doi:10.5061/dryad.s4mw6m9jg, manifest-verified |
| seed | `SEED = {SEED}`; per-analysis seeds fixed in the source |
| permutations | B = {B_PERM:,} per test (resolution {1/(B_PERM+1):.1e}) |
| bootstrap | {NB_BOOT:,} smoothed replicates per Dhat; {R['dhat']['n_cluster_boot']} site-level cluster replicates |
| alpha | {ALPHA} |
| Dhat estimator | Dhat1, selected by overlap's own rule (min n < 50) -- the paper does not name it |
| python / numpy / pandas / scipy | {sys.version.split()[0]} / {np.__version__} / {pd.__version__} / {__import__('scipy').__version__} |

Figures: {", ".join(sorted(FIGS))}
"""

with open(os.path.join(HERE, "REVIEW_REPORT.md"), "w", encoding="utf-8") as fh:
    fh.write(REPORT)
with open(os.path.join(HERE, "numbers.json"), "w", encoding="utf-8") as fh:
    json.dump(R, fh, indent=1, default=str)
det.to_csv(os.path.join(HERE, "detections_derived.csv"), index=False)
print("verdicts:", LEDGER.verdict.value_counts().to_dict())
print("figures :", len(FIGS), "| report:", f"{len(REPORT):,} chars")
