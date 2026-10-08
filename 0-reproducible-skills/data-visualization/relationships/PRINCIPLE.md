# Relationships: "How Do Variables Relate?"

When your question is about the association or correlation between two or more variables.

## Core Principle: Overplotting Hides Density

**When many observations overlap, you lose information about concentration.** A sparse-looking scatterplot might actually be dense where many points collide.

### Example
```
True data: 10,000 cells plotted as points
Visual appearance: Sparse cloud (looks like ~100 cells)
Reason: 95% of cells cluster in one region, overlapped

Solution: Use transparency, hexbin, or 2D density heatmap
Result: Reveals that most cells are in the dense region
```

## Chart Types for This Question

### Scatter Plot (Two Continuous Variables)
**Best for:** Showing relationship between two measurements

**Principles:**
- One point per observation
- X and Y axes are measurements (variables of interest)
- **Scale choice** (linear vs. log) affects interpretation:
  - Linear: shows absolute differences
  - Log: shows relative differences (fold-change)
- Choose scale that answers your question while retaining context

**Handling Overplotting:**
- **Transparency (alpha):** Points stack visually; darker = more overlapping
- **Jitter:** Small random displacement to break up exact overlaps
- **Hexbin or 2D histogram:** Count how many points fall in each region
- **Correlation heatmap:** Shows correlation without showing all points (see below)

**Checklist:**
- [ ] Overplotting strategy chosen (transparency, hexbin, or other)
- [ ] Axis ranges include all data or are justified (e.g., "showing 95th percentile")
- [ ] Units shown on axes
- [ ] If log scale used, justification clear in caption or axis label
- [ ] Subgroups (if any) are distinguishable by color or shape

### Binned Counts (Hexbin, 2D Histogram)
**Best for:** Large datasets (n>1000) where overplotting obscures patterns

**Principles:**
- Divide plot into hexagonal (or square) bins
- Color each bin by count or other summary (mean, median)
- Reveals density without plotting individual points
- **Bin size matters:** Too large = lost detail; too small = noise

**Checklist:**
- [ ] Bin size is specified or justified
- [ ] Color scale appropriate for the summary statistic
- [ ] Legend explains what color represents
- [ ] Sample size noted in caption

### Correlation Matrix
**Best for:** Exploring which variables are related (many variables at once)

**Principles:**
- Rows and columns are variables
- Each cell shows correlation between that row and column variable
- Use diverging color scale (e.g., blue-white-red) centered at 0
- Values range from -1 (perfect negative) to +1 (perfect positive)
- White or light center = no correlation

**Checklist:**
- [ ] Diverging color scale used (not sequential)
- [ ] Correlation method specified (Pearson, Spearman, etc.)
- [ ] If p-values shown, test specified and threshold stated
- [ ] Sample size noted (especially if computed within subgroups)

## Common Pitfalls

### Pitfall 1: Overplotting (Invisible Points)
**Problem:** Scatterplot with 10,000 points looks empty; relationships hard to see.
**Solution:** Use transparency, hexbin/2D histogram, or density contours.

```python
# Bad: opaque points, overplotting obscures data
ax.scatter(x, y)  

# Good: transparent points, density visible
ax.scatter(x, y, alpha=0.1)

# Also good: hexbin reveals density directly
ax.hexbin(x, y, cmap='YlOrRd')
```

### Pitfall 2: Inappropriate Scale
**Problem:** Linear plot when log scale is needed (or vice versa).
**Solution:** Choose scale that answers your question and shows data well.

Examples:
- Gene expression vs. dose: **log scale** (fold-change is meaningful)
- Age vs. outcome: **linear scale** (absolute difference is meaningful)
- Very skewed data: **log or sqrt scale** (to show the range)

### Pitfall 3: Ignoring Subgroups
**Problem:** Single correlation across all data when subgroups have different relationships.
**Solution:** Color or facet by subgroup. Or, report correlation within each group.

Example: CD8 vs. tumor burden might be negative in responders but positive in non-responders.

### Pitfall 4: Sample Size Imbalance Affecting Interpretation
**Problem:** Correlation driven by one subgroup (which happens to be large).
**Solution:** Show sample size for each subgroup. Color by subgroup. Report correlation separately per group if relevant.

## Data Inspection Checklist

Before creating any relationship visualization:

```python
from utilities.data_inspection import inspect_variable

# 1. Understand each variable separately
inspect_variable(data['var1'])  # → range, distribution, missingness
inspect_variable(data['var2'])

# 2. Check for outliers in each variable
# These can dominate correlation and distort visual appearance

# 3. Check for subgroups
# Correlation structure often differs by group

# 4. Assess sample size
# Large n needed for reliable correlation estimates
assert len(data) >= 50, "Sample too small for correlation analysis"
```

## When to Use Each Chart

| Situation | Chart Type | Why |
|-----------|-----------|-----|
| Two variables, n<500, want to see each point | Scatter with transparency or jitter | Shows all data; density visible via overlap |
| Two variables, n>1000, overplotting severe | Hexbin or 2D histogram | Clearly shows density without overplotting |
| Many variables (5+), want correlations | Correlation matrix heatmap | Shows all pairwise relationships at once |
| Two variables, clear subgroups | Scatter colored by subgroup | Reveals whether relationship changes by group |
| Want to show confidence/uncertainty in trend | Scatter + regression line + confidence band | Shows both data and fitted trend |
| Many variables, want to see all points AND patterns | Small multiples (faceted scatter plots) | Each variable pair in separate panel; shared axes |

## Example Workflow

```python
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from utilities.data_inspection import inspect_variable
from utilities.color_scales import categorical_palette

# Step 1: Inspect variables
inspect_variable(data['cd8_expression'])    # log₂(TPM), range 0–15
inspect_variable(data['tumor_burden'])      # cells/mm³, range 100–10000

# Step 2: Check for outliers and subgroups
# By group: responders vs. non-responders
n_resp = len(data[data['response'] == 'responder'])
n_nonresp = len(data[data['response'] != 'responder'])
# Both > 200 → can split plot by group

# Step 3: Create scatterplot with subgroup coloring
fig, ax = plt.subplots(figsize=(6, 5))
colors = categorical_palette(n_categories=2)

for i, (group, color) in enumerate(zip(['responder', 'non-responder'], colors)):
    subset = data[data['response'] == group]
    ax.scatter(subset['cd8_expression'], subset['tumor_burden'], 
               alpha=0.3, s=20, color=color, label=f"{group} (n={len(subset)})")

ax.set_xlabel("CD8 expression (log₂ TPM)")
ax.set_ylabel("Tumor burden (cells/mm³)")
ax.legend()

# Step 4: Add summary statistics
# Compute correlation within each group
for group, color in zip(['responder', 'non-responder'], colors):
    subset = data[data['response'] == group]
    r, p = spearmanr(subset['cd8_expression'], subset['tumor_burden'])
    print(f"{group}: ρ = {r:.2f}, p = {p:.3f}")

fig.savefig("cd8_vs_burden_by_response.png", dpi=300, bbox_inches='tight')
```

---

See `scatter_and_binned.py` and `correlation_matrices.py` for utility functions.
