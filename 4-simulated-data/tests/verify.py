import pandas as pd, numpy as np, statsmodels.formula.api as smf
from statsmodels.stats.diagnostic import het_breuschpagan
from pathlib import Path
DATA = (Path(__file__).resolve().parents[1]
        / "data" / "processed" / "part1-mine-axolotl" / "axolotl_regeneration.csv")
d = pd.read_csv(DATA)
d["fgf8_c"]=d.fgf8_log2cpm-d.fgf8_log2cpm.mean(); d["shh_c"]=d.shh_log2cpm-d.shh_log2cpm.mean()
print("clipped at 100:", (d.regen_pct_day30>=100).sum(), " at 0:", (d.regen_pct_day30<=0).sum())
print("\nMarginal correlations with outcome:")
print(d[["age_months","water_temp_C","fgf8_log2cpm","shh_log2cpm","regen_pct_day30"]].corr()["regen_pct_day30"].round(3))
def show(f):
    m=smf.ols(f,d).fit(); print(f"\n{f}\n  R2={m.rsquared:.3f}"); print(m.params.round(3).to_string()); return m
show("regen_pct_day30 ~ age_months")
show("regen_pct_day30 ~ age_months + C(colony)")
naive=show("regen_pct_day30 ~ age_months + water_temp_C + fgf8_c + shh_c + C(colony)")
full=show("regen_pct_day30 ~ age_months + water_temp_C + I((water_temp_C-20)**2) + fgf8_c*shh_c + C(colony)")
print("\nBreusch-Pagan p (full model):", het_breuschpagan(full.resid, full.model.exog)[1])
r=full.resid; print("resid SD by colony:", r.groupby(d.colony).std().round(2).to_dict())
print("SE of age coef: alone", smf.ols("regen_pct_day30 ~ age_months",d).fit().bse.age_months.round(3),
      " with colony", full.bse.age_months.round(3))
# bootstrap instability of age coefficient
rng=np.random.default_rng(1); a=[]
for _ in range(300):
    s=d.sample(len(d),replace=True,random_state=int(rng.integers(1e9)))
    a.append(smf.ols("regen_pct_day30 ~ age_months + C(colony)",s).fit().params.age_months)
print("bootstrap age coef (with colony): 2.5/97.5%", np.percentile(a,[2.5,97.5]).round(2))
print("age range by colony:", d.groupby("colony").age_months.agg(["min","max"]).round(1).to_dict())
