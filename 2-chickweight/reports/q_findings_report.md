# ChickWeight Analysis: Findings Report

## Dataset Overview

- **Source:** Classic dataset in R; 50 chicks over 21 days
- **Treatments:** 4 diets
- **Sample sizes:** Diet 1 (n=20), Diet 2 (n=10), Diet 3 (n=10), Diet 4 (n=9)
- **Measurements:** 12 time points (days 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 21)

---

## Q1: How do growth trajectories differ by diet?

**Key finding:** All diets show sigmoidal (S-shaped) growth curves, but the **rate and final outcome differ substantially**.

### Mean final weight ranking:
1. **Diet 3: 270.3g** — Fastest growth, highest final weight
2. **Diet 4: 238.6g** — Intermediate growth
3. **Diet 2: 214.7g** — Moderate growth
4. **Diet 1: 177.8g** — Slowest growth

### Total weight gain over 21 days:
- Diet 3: +229.5g (5.6× starting weight)
- Diet 4: +197.6g (4.8× starting weight)
- Diet 2: +174.0g (4.3× starting weight)
- Diet 1: +136.3g (3.3× starting weight)

**Interpretation:** Diet 3 is superior for overall growth. The trajectories diverge most in weeks 2–3, suggesting diet-dependent growth acceleration in later development.

---

## Q2: Which diet achieves the best balance of high final weight and low variability?

**Key finding:** **Diet 4 provides the best consistency** (lowest variability), while **Diet 3 achieves highest final weight at moderate cost to variability**.

### Variability metrics (final weight, day 21):

| Diet | Mean ± SD | Coefficient of Variation | Range (min–max) |
|------|-----------|-------------------------|-----------------|
| 1 | 177.8 ± 58.7g | 33% | 96–305g |
| 2 | 214.7 ± 78.1g | 36% | 74–331g |
| 3 | 270.3 ± 71.6g | 27% | 147–373g |
| 4 | 238.6 ± 43.3g | **18%** | 196–322g |

**Outliers:** Only Diet 1 produced one extreme outlier (Chick 7: 305g, +127g above mean). All other diets show consistent responses.

**Trade-off analysis:**
- **Best for high weight:** Diet 3 (270g mean, 27% CV)
- **Best for consistency:** Diet 4 (18% CV, stable 196–322g range)
- **Worst:** Diet 2 shows high variability (36% CV) without compensatory high final weight

---

## Q3: Do chicks grow at constant rates, or do growth patterns change over time?

**Key finding:** **Growth accelerates in the mid-period (8–14 days) for all diets**, then either plateaus (Diet 1) or continues accelerating (Diets 2, 3, 4).

### Daily weight gain by period:

#### Early (0–8 days):
- Diet 4: 8.07 ± 0.36 g/day (fastest start)
- Diet 3: 7.20 ± 0.50 g/day
- Diet 2: 6.38 ± 0.61 g/day
- Diet 1: 4.51 ± 0.45 g/day (slowest start)

#### Mid (8–14 days): **Major acceleration zone**
- Diet 3: 11.02 ± 1.20 g/day (↑54% from early)
- Diet 4: 9.37 ± 0.60 g/day (↑16% from early)
- Diet 2: 8.37 ± 1.57 g/day (↑31% from early)
- Diet 1: 6.68 ± 1.12 g/day (↑48% from early)

#### Late (14–21 days): **Divergent patterns**
- Diet 3: 15.11 ± 2.19 g/day (↑37% from mid — **continued acceleration**)
- Diet 2: 10.40 ± 2.01 g/day (↑24% from mid)
- Diet 4: 9.64 ± 1.95 g/day (↑3% from mid — stable growth)
- Diet 1: 6.48 ± 1.20 g/day (↓3% from mid — **plateau/deceleration**)

**Interpretation:**
- **Diet 1** reaches a growth plateau by mid-period, limiting final size
- **Diet 3** shows sustained acceleration throughout, driving highest final weight
- **Diet 4** shows strong acceleration early and mid, then stabilizes at high level
- **Diet 2** shows variability in growth rates (high SEM in all periods)

---

## Q4: Which individual chicks are outliers, and do certain diets show more variability?

**Key finding:** **Diet 1 and Diet 2 show high individual variability**, while **Diet 4 produces highly consistent responses**.

### Outlier count (>2 SD from diet mean):
- **Diet 1:** 1 outlier (Chick 7: 305g, +71% above mean) — "failure to thrive" and one exceptional grower
- **Diet 2:** 0 outliers but high overall spread (74–331g range)
- **Diet 3:** 0 outliers, tight consistency at high mean weight
- **Diet 4:** 0 outliers, tightest range (196–322g)

### Interpretation by diet:
- **Diet 1:** Produces most variable responses; includes both slow and exceptional growers
- **Diet 2:** Supports growth but with unpredictable individual outcomes (widest range despite no outliers)
- **Diet 3:** High-quality diet producing consistent, superior growth across individuals
- **Diet 4:** High-quality diet producing consistent, moderate growth (reliable if extreme size not needed)

---

## Summary: Answering All 4 Questions

| Question | Answer |
|----------|--------|
| **Q1: Trajectory differences?** | Diet 3 grows fastest; Diet 1 slowest. Divergence accelerates in weeks 2–3. |
| **Q2: Best balance?** | Diet 4 (consistency, 18% CV) vs. Diet 3 (high weight, 27% CV). Diet 2 poor. |
| **Q3: Constant growth rate?** | No—growth accelerates in mid-period. Diet 3 sustains acceleration; Diet 1 plateaus. |
| **Q4: Individual variation?** | Diet 4 & 3 consistent. Diet 1 & 2 variable. Only Diet 1 produces outliers. |

### Recommendation:
- **For maximum growth:** Use Diet 3
- **For reliable, consistent growth:** Use Diet 4
- **Avoid:** Diet 2 (high variability, no advantage)
- **Note:** Diet 1 baseline is poor quality; variable responses suggest nutrient imbalance

---

**Analysis date:** October 2026
**Visualizations:** 4 figures (q1–q4)
**Methods:** Descriptive statistics, group-level means with SEM, individual trajectory overlay, outlier detection (>2 SD)
