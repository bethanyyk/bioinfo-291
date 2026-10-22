"""
Part 2 -- analysis of the partner's agonist-antibody screen.

Input : data/raw/part2-partner-antibody/antibody_agonism.csv  (received, not modified)
Output: results/*.csv, figures/fig_antibody_agonism.png

Model: Tweedie GLM with log link (var_power selected by profile log-likelihood),
which handles the 113 exact zeros, the multiplicative site/kinetics structure and
the variance-tracks-mean behaviour that defeat an ordinary linear model.

Run: python workflows/analyze_antibody_agonism.py
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import statsmodels.api as sm, statsmodels.formula.api as smf
from scipy import stats
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "part2-partner-antibody" / "antibody_agonism.csv"
RES, FIG = ROOT / "results", ROOT / "figures"
SEED = 11

d = pd.read_csv(DATA)
for f in (d,):
    f["koff_c"] = f.log10_koff + 3.0      # centred near the optimum
    f["kon_c"]  = f.log10_kon - 5.5
    f["dens_c"] = f.receptor_density - 1.0

# ---- hold out 25%, stratified by site, BEFORE any model choice ------------
rng = np.random.default_rng(SEED)
test_idx = np.concatenate([rng.permutation(g.index.values)[:25] for _, g in d.groupby("site")])
te, tr = d.loc[np.sort(test_idx)].copy(), d.drop(index=test_idx).copy()

# ---- Tweedie variance power by profile log-likelihood on TRAIN -----------
SPEC = "agonism_pct ~ C(site) + kon_c + koff_c + I(koff_c**2) + dens_c"
prof = [(p, smf.glm(SPEC, tr, family=sm.families.Tweedie(
            link=sm.families.links.Log(), var_power=float(p))).fit().llf)
        for p in np.round(np.arange(1.10, 1.96, 0.05), 2)]
P_HAT = float(max(prof, key=lambda t: t[1])[0])
fam = sm.families.Tweedie(link=sm.families.links.Log(), var_power=P_HAT)

FINAL = SPEC + " + C(plate)"
WITH_INT = "agonism_pct ~ C(site)*(kon_c + koff_c + I(koff_c**2)) + dens_c + C(plate)"

# ---- held-out comparison --------------------------------------------------
def score(m):
    p = np.asarray(m.predict(te), float); y = np.asarray(te.agonism_pct, float)
    return dict(rmse=np.sqrt(np.mean((y-p)**2)), mae=np.mean(np.abs(y-p)),
                r2=1-np.sum((y-p)**2)/np.sum((y-y.mean())**2),
                tweedie_dev=float(fam.deviance(y, np.clip(p, 1e-9, None))),
                spearman=float(stats.spearmanr(y, p).statistic),
                top20_hits=int(len(set(np.argsort(-p)[:20]) & set(np.argsort(-y)[:20]))))

models = {"SLR: agonism ~ log10_KD (OLS)": smf.ols("agonism_pct ~ log10_KD", tr).fit(),
          "OLS: site+kin+interactions":    smf.ols("agonism_pct ~ C(site)*(kon_c + koff_c + I(koff_c**2)) + dens_c", tr).fit(),
          "Tweedie GLM, no plate":         smf.glm(SPEC, tr, family=fam).fit(),
          "Tweedie GLM + plate":           smf.glm(FINAL, tr, family=fam).fit(),
          "Tweedie GLM + plate + site:kin":smf.glm(WITH_INT, tr, family=fam).fit()}
pd.DataFrame([dict(model=k, **score(v)) for k, v in models.items()]).round(3) \
  .to_csv(RES/"model_comparison_heldout.csv", index=False)

# ---- is the site x kinetics interaction real, or the zero floor? ----------
# Restricted to sites 2-4 (where the floor barely binds) the kon slope is homogeneous.
mid = tr[tr.site.isin([2, 3, 4])]
mh = smf.glm("agonism_pct ~ C(site) + kon_c + koff_c + I(koff_c**2) + dens_c + C(plate)", mid, family=fam).fit()
mt = smf.glm("agonism_pct ~ C(site) + C(site):kon_c + koff_c + I(koff_c**2) + dens_c + C(plate)", mid, family=fam).fit()
F = ((mh.deviance - mt.deviance) / (mt.df_model - mh.df_model)) / (mt.pearson_chi2 / mt.df_resid)
print(f"sites 2-4, kon slope homogeneity: F={F:.2f}, "
      f"p={stats.f.sf(F, mt.df_model-mh.df_model, mt.df_resid):.3f}")

# ---- final fit on all 500 for reported coefficients -----------------------
mf = smf.glm(FINAL, d, family=fam).fit()
ci = mf.conf_int()
pd.DataFrame({"term": mf.params.index, "beta_log": mf.params.values,
              "lo": ci[0].values, "hi": ci[1].values,
              "mult": np.exp(mf.params.values), "p": mf.pvalues.values}).round(4) \
  .to_csv(RES/"final_model_coefficients.csv", index=False)

b1, b2 = mf.params["koff_c"], mf.params["I(koff_c ** 2)"]
V = mf.cov_params().loc[["koff_c","I(koff_c ** 2)"],["koff_c","I(koff_c ** 2)"]].values
grad = np.array([-1/(2*b2), b1/(2*b2**2)]); se = float(np.sqrt(grad @ V @ grad))
opt = -b1/(2*b2) - 3.0
print(f"Tweedie p={P_HAT}; koff optimum log10={opt:.2f} +/- {1.96*se:.2f}; "
      f"kon x{np.exp(mf.params['kon_c']):.2f} per log10")

# ---- plate-corrected shortlist -------------------------------------------
pl_eff = {1: 0.0, **{p: float(mf.params[f"C(plate)[T.{p}]"]) for p in range(2, 7)}}
cnt = d.plate.value_counts().sort_index()
offset = float(np.average([pl_eff[p] for p in cnt.index], weights=cnt.values))
d["agonism_plate_adj"] = d.agonism_pct / d.plate.map({p: np.exp(pl_eff[p]-offset) for p in range(1,7)})
# model-predicted agonism for the same antibody run on a typical plate
d["pred_typical_plate"] = np.exp(np.log(mf.predict(d)) - (d.plate.map(pl_eff) - offset))
d.sort_values("agonism_plate_adj", ascending=False).head(20)[
    ["ab_id","site","epitope_dist_nm","log10_kon","log10_koff","log10_KD","plate",
     "agonism_pct","agonism_plate_adj","pred_typical_plate"]].round(3).to_csv(
        RES/"top20_candidates.csv", index=False)
print("wrote results/ tables; see reports/antibody_agonism_report.md for the write-up")
