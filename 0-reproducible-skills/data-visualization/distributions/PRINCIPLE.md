# Distributions: "How Does It Vary?"

When your question is about the distribution or variability of a single variable across groups.

## Core Principle: Identical Means ≠ Similar Observations

**Showing only the mean and standard deviation can hide multimodality, skewness, and outliers.** Always show the data or a distribution summary.

### Example: Simpson's Paradox in Distributions
```
Three datasets with identical means and SDs:
Dataset A: [8.04, 6.95, 7.58, 8.81, 8.33, ...]  → bell curve
Dataset B: [9.14, 8.14, 8.74, 8.77, 9.26, ...]  → y = 3 + 0.5x
Dataset C: [7.46, 6.77, 12.74, 7.11, 8.84, ...]  → with one outlier

All have mean ≈ 7.5, SD ≈ 2.03, but the distributions are very different!
```

**Solution:** Plot the distribution, not just the summary statistic.

## Chart Types for This Question

### Exploratory (Show All Data)
Best when n is small to moderate (<500 observations per group).

#### Points (Strip Plot, Jitter Plot)
**Best for:** Raw observations; showing count and spread simultaneously

**Principles:**
- Each point is an observation
- Add jitter (random small displacement) to reduce overplotting
- Can add a summary (mean/median line) on top
- Color or shape can distinguish subgroups

**Checklist:**
- [ ] Overplotting is visible and not hiding observations
- [ ] Jitter amount is reasonable (not obscuring pattern)
- [ ] Sample size is shown in caption or on axis
- [ ] If overlapping, consider transparency

#### Histogram
**Best for:** Understanding bin structure and modality

**Principles:**
- **Bin width matters.** Different widths reveal different patterns:
  - Too narrow: noisy (every bin has 0–1 observations)
  - Too wide: smooth but loses detail
- Start with Freedman-Diaconis rule or Sturges' rule; then adjust
- Show count or density on y-axis (specify which)
- Include axis labels with units

**Checklist:**
- [ ] Bin width is specified or justified (caption or code)
- [ ] Y-axis is labeled (count vs. density)
- [ ] All data are included (no silent trimming)
- [ ] Color choice: single color or gradient if showing groups separately

#### Empirical Cumulative Distribution Function (ECDF)
**Best for:** Non-parametric comparison across groups

**Principles:**
- X-axis: variable of interest
- Y-axis: proportion of data ≤ x
- One line per group (different colors)
- Useful for comparing distributions without assuming normal distribution

**Checklist:**
- [ ] Lines are distinct (different colors or line styles)
- [ ] Legend identifies each group
- [ ] X-axis range includes all data
- [ ] Y-axis ranges from 0 to 1 (or 0–100%)

### Summary Shapes (Show Distribution without All Points)
Best when n is large (>500 observations per group) or when a few summary statistics are sufficient.

#### Box Plot
**Best for:** Quick comparison of medians, quartiles, and outliers

**Principles:**
- Box: IQR (25th to 75th percentile)
- Line in box: median (50th percentile)
- Whiskers: extend to ~1.5 × IQR
- Points beyond whiskers: outliers
- **Limitation:** Hides multimodality; two bimodal distributions can look identical

**Checklist:**
- [ ] Outlier threshold is standard (1.5 × IQR) or specified
- [ ] Multiple groups use distinct colors
- [ ] Y-axis includes all data (no trimming)
- [ ] Sample size per group is shown (if imbalanced)

#### Violin Plot
**Best for:** Showing full distribution shape without showing every point

**Principles:**
- Mirrored kernel density estimate
- Width at each y-value ∝ density at that value
- Often combined with box plot or points for detail
- **Caveat:** Smoothing bandwidth choice affects appearance

**Checklist:**
- [ ] Bandwidth (smoothing) is documented or justified
- [ ] Not misinterpreted as continuous (it's a density estimate)
- [ ] Paired with box plot or points for quantiles
- [ ] Sample size per group shown (if imbalanced)

#### Combination: Points + Box + Violin
**Best for:** Maximum information with minimal clutter (when n is moderate, ~50–500)

**Principles:**
- Violin: overall shape
- Box: quartiles and median
- Points: individual observations (with transparency or jitter if overlapping)
- Layers convey different information at different visual weights

**Checklist:**
- [ ] All three layers add information (remove if redundant)
- [ ] Transparency on points allows seeing overlap
- [ ] Colors distinguish groups, not layers
- [ ] Not overly complex (save complex figures for exploratory analysis)

## Common Pitfalls

### Pitfall 1: Bimodal Distribution Masked as Normal
**Problem:** Two groups (e.g., expression in two cell types) pooled and shown as single distribution.
**Solution:** Facet by group, or use color/shape to distinguish subgroups.

### Pitfall 2: Bin Width Hiding the Pattern
**Problem:** Histogram with default bins obscures bimodality; violin with default bandwidth over-smooths.
**Solution:** Explore multiple bin widths or bandwidths. Document your choice.

### Pitfall 3: Outliers Dominating the Axis
**Problem:** One extreme outlier stretches axis, making the bulk of the distribution hard to read.
**Solution:**
- Trim x-axis to show 95th percentile, and note outliers in caption
- Or use log scale if appropriate
- Or use box plot (automatically handles outliers)

### Pitfall 4: Showing Only Mean ± SD
**Problem:** Hides multimodality, skewness, and non-normal tails.
**Solution:** Show the actual distribution (histogram, ECDF, violin, or points).

## Data Inspection Checklist

Before creating any distribution visualization:

```python
from utilities.data_inspection import inspect_variable

# 1. Understand the variable
inspect_variable(data['expression'], groups=data['cell_type'])
# → Shows distribution, outliers, missingness, by group

# 2. Check for bimodality or skewness
import numpy as np
assert not (data['expression'].skew() > 2), "Extreme skewness; consider log transform"

# 3. Decide on transformation
# Log scale useful for: count data, right-skewed, fold-changes
# Square root: for count data
# Box-Cox: for automatic transform selection
```

## When to Use Each Chart

| Situation | Chart Type | Why |
|-----------|-----------|-----|
| Few observations (n<100), need to see each point | Points + jitter | Preserves all data; no information loss |
| Moderate observations (n=100–500) | Points + box + violin | Shows distribution without overwhelming |
| Many observations (n>500), publication | Violin + box | Clean summary; avoids overplotting |
| Comparing distributions non-parametrically | ECDF | No assumption of normality |
| Need to show outliers explicitly | Box plot | Automatically marks outliers |
| Quick exploratory look at multiple groups | Multiple histograms | Shows modality, skewness per group |
| Distribution is heavy-tailed or skewed | Box plot or ECDF | Better than assuming normal |

## Example Workflow

```python
import matplotlib.pyplot as plt
from scipy import stats
from utilities.data_inspection import inspect_variable
from utilities.color_scales import categorical_palette

# Step 1: Inspect data by group
inspect_variable(data['cd8_expression'], groups=data['response'])
# → Shows distribution, n per group, range

# Step 2: Choose chart type based on n and question
n_responders = len(data[data['response'] == 'responder'])
n_non = len(data[data['response'] != 'responder'])
# Both > 200 → violin + box; show summary not individual points

# Step 3: Create visualization
fig, ax = plt.subplots(figsize=(5, 6))
colors = categorical_palette(n_categories=2)

groups = ['responder', 'non-responder']
data_by_group = [data[data['response'] == g]['cd8_expression'] for g in groups]

# Box + violin
parts = ax.violinplot(data_by_group, positions=[0, 1], showmeans=False, showmedians=False)
bp = ax.boxplot(data_by_group, positions=[0, 1], widths=0.3, patch_artist=True)

# Color
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

ax.set_ylabel("CD8 expression (log₂ TPM)")
ax.set_xlabel("Response status")
ax.set_xticklabels(groups)

# Add sample sizes
for i, g in enumerate(groups):
    n = len(data[data['response'] == g])
    ax.text(i, ax.get_ylim()[0] - 1, f"n={n}", ha='center', fontsize=9)

fig.savefig("cd8_distribution_by_response.png", dpi=300, bbox_inches='tight')
```

---

See `exploratory.py` and `summary_shapes.py` for utility functions.
