"""
The Great Axolotl Regrowth Race -- simulated dataset (seeded, reproducible).
500 axolotls, 5 predictors, 1 outcome. Run: python simulate_axolotl.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

OUT_CSV = (Path(__file__).resolve().parents[1]
           / "data" / "processed" / "part1-mine-axolotl" / "axolotl_regeneration.csv")

rng = np.random.default_rng(2026)

# --- Colonies (recorded groups) -------------------------------------------
# Older, genetically diverse wild-rescue colony regenerates best;
# young inbred GFP line regenerates worst. Colony mean ages barely overlap.
colonies = {
    #  name               n    mean_age  sd_age  baseline_regen
    "GFP_inbred":       (150,   8.0,     2.5,    48.0),
    "Classic_lab":      (200,  16.0,     3.0,    62.0),
    "Xochimilco_rescue":(150,  30.0,     3.5,    76.0),
}

rows = []
for name, (n, mu_age, sd_age, base) in colonies.items():
    age = np.clip(rng.normal(mu_age, sd_age, n), 3, None)
    rows.append(pd.DataFrame({"colony": name, "age_months": age,
                              "_age_dev": age - mu_age, "_base": base}))
df = pd.concat(rows, ignore_index=True)
N = len(df)

# --- Tank temperature: symmetric around the optimum (20 C) ----------------
df["water_temp_C"] = rng.uniform(12, 28, N)

# --- Blastema gene expression (log2 CPM), independent of each other -------
fgf8_z = rng.normal(0, 1, N)
shh_z = rng.normal(0, 1, N)
df["fgf8_log2cpm"] = 6.0 + 1.0 * fgf8_z
df["shh_log2cpm"] = 4.0 + 1.2 * shh_z

# --- True mean model -------------------------------------------------------
mu = (df["_base"]
      - 1.0 * df["_age_dev"]                    # within-colony: older = slower
      - 0.30 * (df["water_temp_C"] - 20) ** 2   # inverted-U temperature optimum
      + 7.0 * fgf8_z * shh_z)                   # FGF8-SHH coordination

# --- Noise grows with age (heteroscedastic) --------------------------------
sd = 2.0 + 0.30 * df["age_months"]
df["regen_pct_day30"] = (mu + rng.normal(0, sd)).clip(0, 100)

df.insert(0, "axolotl_id", [f"AX{i+1:03d}" for i in range(N)])
df = df.sample(frac=1, random_state=7).reset_index(drop=True)  # shuffle rows
out = df[["axolotl_id", "colony", "age_months", "water_temp_C",
          "fgf8_log2cpm", "shh_log2cpm", "regen_pct_day30"]].round(2)
out.to_csv(OUT_CSV, index=False, lineterminator="\n")  # LF, so reruns are byte-identical on Windows too
print(out.head(), "\n", out.describe())
