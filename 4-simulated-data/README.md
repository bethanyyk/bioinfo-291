# The Great Axolotl Regrowth Race

## The story

The Ambystoma Regeneration Consortium amputated one forelimb from each of **500 axolotls**. Thirty days later, they measured how much of the limb had grown back. The animals came from three colonies:

- **GFP_inbred**: a young, inbred lab line that glows green.
- **Classic_lab**: the standard lab stock.
- **Xochimilco_rescue**: older, genetically diverse animals rescued from the canals of Mexico City.

Each animal lived in a tank set somewhere between 12 and 28 °C. On day 7, a tiny biopsy of the *blastema* (the bud of stem-like cells that rebuilds the limb) was sequenced for two key signalling genes: **Fgf8** and **Shh**.

A grant reviewer wrote: *"Your data clearly show that older axolotls regenerate better, and temperature and gene expression don't matter. Please confirm."*

**Your job:** build a model of regeneration and decide whether the reviewer is right.

## Data dictionary (`axolotl_regeneration.csv`, 500 rows)

| Column | Meaning |
| --- | --- |
| `axolotl_id` | Animal ID |
| `colony` | Source colony (GFP_inbred, Classic_lab, Xochimilco_rescue) |
| `age_months` | Age at amputation (months) |
| `water_temp_C` | Tank water temperature (°C) |
| `fgf8_log2cpm` | Fgf8 expression in day-7 blastema (log2 CPM) |
| `shh_log2cpm` | Shh expression in day-7 blastema (log2 CPM) |
| **`regen_pct_day30`** | **Outcome**: % of original limb length regrown by day 30 (capped at 100) |

Regenerate the data with `python simulate_axolotl.py` (seed 2026). `verify.py` reproduces the checks in the answer key.

---

## Answer key (spoilers, for the instructor)

The true generating model:

```
regen = colony_baseline − 1.0·(age − colony_mean_age)
        − 0.30·(temp − 20)²
        + 7.0·z_fgf8·z_shh
        + Normal(0, sd = 2 + 0.3·age)
```

| Trap | How it's built in | What a naive analysis finds | What the right model finds |
| --- | --- | --- | --- |
| **Different relationships between groups** (Simpson's paradox) | The older colony has a much higher baseline, but within every colony older animals regenerate *slower*. | `regen ~ age`: slope **+0.94**/month. Older looks better. | Adding `colony` flips the age slope to **−1.07**/month (true −1.0). |
| **Curvature** | Inverted U with an optimum at 20 °C | Correlation of temperature with regeneration is **r = −0.10**. Temperature looks irrelevant. | The quadratic term is **−0.29** (true −0.30). Temperature is a major driver, and R² jumps. |
| **Interaction** | Fgf8 and Shh only help when *both* are high. Their main effects are zero. | Marginal correlations are **−0.05 and 0.00**. Neither gene looks relevant. | `fgf8 × shh` coefficient is **+5.97**. This is the biggest single effect in the model. |
| **Unequal noise** | Residual SD grows with age | Point estimates look fine. | Breusch–Pagan **p ≈ 3×10⁻¹⁴**. Residual SD is 4.6 in GFP_inbred vs 10.3 in Xochimilco. Use robust (HC3) SEs, WLS, or a variance model. |
| **Correlated predictors** | Colony age ranges barely overlap (GFP 3–16, Classic 5–25, Xochimilco 19–40 months). | — | The age SE goes from 0.064 to 0.108 once colony is added. Age and colony compete to explain the same variation, so the separate "age effect" is less precise than the fit suggests. |

Model fit: the naive linear model has **R² = 0.50**. The correct model (colony + age + temp² + fgf8×shh) has **R² = 0.77**.

**Verdict on the reviewer:** wrong on all three counts. Within each colony, older axolotls regenerate *slower*. Temperature matters a lot, but its effect is a curve rather than a straight line. The two genes matter jointly but not separately.

Biological plausibility: regeneration really does slow with age in axolotls, and Fgf8 and Shh really do act as a mutually dependent signalling loop in the regenerating limb (Nacu et al., *Nature* 2016). The effect sizes in this dataset are invented.

---

## Repository layout

This project holds **two datasets with different provenance**, kept in separate trees:

**Part 1 — generated here, sent to the partner** (axolotl limb regeneration)

| Path | Contents |
| --- | --- |
| `workflows/simulate_axolotl.py` | Generates the dataset (seed 2026). |
| `data/processed/part1-mine-axolotl/axolotl_regeneration.csv` | The 500-animal dataset sent to the partner. |
| `tests/verify.py` | Reproduces the answer-key checks above. |

**Part 2 — received from the partner** (agonist-antibody screen). Originals, not modified:

| Path | Contents |
| --- | --- |
| `data/raw/part2-partner-antibody/antibody_agonism.csv` | 500 candidate antibodies, as received. |
| `data/raw/part2-partner-antibody/README.md` | The partner's story, data dictionary and task. |
| `src/`, `workflows/`, `results/`, `figures/`, `reports/` | My analysis of the partner's data. |

`data/metadata/checksums.txt` records SHA-256 for both datasets and labels which is which.

Both scripts resolve paths relative to the project root, so they can be run from any
working directory:

```bash
python workflows/simulate_axolotl.py   # rewrites the CSV, byte-identical
python tests/verify.py                 # prints the answer-key numbers
```
