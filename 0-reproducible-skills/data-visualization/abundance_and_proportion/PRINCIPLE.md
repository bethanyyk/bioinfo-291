# Abundance & Proportion: "How Many? What Fraction?"

When your question is about counts, frequencies, or proportions across groups.

## Core Principle: Proportional Ink

**When a shaded region is used to represent a numerical value, the area of that shaded region must be directly proportional to the corresponding value.**

This principle is fundamental to bar charts and proportion visualizations.

### Violation Example (Bad)
```
Actual values: [100, 90, 10]
Truncated bar chart (y-axis from 80 to 100):
[████ 100]  [████ 90]   [█ 10]  ← Visually implies 100 ≈ 90 ≈ 10 (false!)
```

### Correct Example (Good)
```
Actual values: [100, 90, 10]
Full bar chart (y-axis from 0 to 100):
[████████████ 100]  [███████████ 90]   [█ 10]  ← True proportional ink
```

## Chart Types for This Question

### Count Bar Charts
**Best for:** Comparing numbers of items in different categories

**Principles:**
- Y-axis should start at 0 (proportional ink)
- Bar length represents magnitude
- One bar per category
- Categories can be ordered (alphabetical, by size, by importance) or unordered

**Checklist:**
- [ ] Y-axis starts at 0
- [ ] Bar height matches the value
- [ ] Categories are clearly labeled
- [ ] If many categories (>10), consider sorting by value or faceting

### Proportion/Percentage Bars
**Best for:** Showing parts of a whole (stacked or side-by-side)

**Principles:**
- **Always show the denominator.** A percentage means nothing without context.
  - ✓ "30% of 1,000 cells"
  - ✗ "30% of cells"
- Bar represents 100% (or whatever the whole is)
- Colors distinguish categories (see `utilities/color_scales.py`)
- When stacked, hard-to-read categories go in the middle (easiest to read: left, right, then middle)

**Checklist:**
- [ ] Caption states denominator (total n, what the 100% includes)
- [ ] Total bar width/height is comparable across groups (if comparing multiple groups)
- [ ] Colors are distinct and appropriate for categorical data
- [ ] Ordering (stacked order) is logical or serves the narrative

### Mosaic Plots
**Best for:** Showing how group sizes affect proportions (preserves sample composition information)

**Principles:**
- Width of bar represents group size
- Height represents proportion within that group
- Area of rectangle represents the count
- Preserves information about group imbalance

**Checklist:**
- [ ] Group sizes are clear from bar widths
- [ ] Legend or annotation explains what is being shown
- [ ] Color represents categories being compared
- [ ] Caption explains sample composition (which groups, which n)

## Common Pitfalls

### Pitfall 1: Truncated Y-Axis (Proportional Ink Violation)
**Problem:** Bar chart with y-axis from 80 to 100 makes 90 look much smaller than 100.
**Solution:** Start y-axis at 0, or use a dot/slope chart instead (these don't have the proportional ink requirement).

### Pitfall 2: Missing Denominator
**Problem:** "CD8+ frequency increased 2×" — 2× from what baseline? In what sample?
**Solution:** Always state the denominator. "CD8+ frequency: 15% in pre-treatment (n=200) vs. 30% in on-treatment (n=250)"

### Pitfall 3: Stacked Bar with Too Many Categories
**Problem:** 15 categories stacked in one bar → middle categories are unreadable.
**Solution:** 
- Separate into multiple facets (one bar per panel)
- Use a legend with distinct colors
- Or, use grouped bars instead of stacked

### Pitfall 4: Sample Imbalance Ignored
**Problem:** Bar chart comparing "responders" (n=20) and "non-responders" (n=150) looks like equal bars.
**Solution:** Show sample sizes on the figure or in caption. Consider mosaic plot to visualize imbalance.

## Data Inspection Checklist

Before creating any count/proportion visualization:

```python
from utilities.data_inspection import sample_composition, check_proportions

# 1. Understand sample composition
sample_composition(data, group_col='treatment')
# → Shows group sizes, imbalance, nesting (if any)

# 2. Check proportions
check_proportions(data, numerator='cd8_positive', denominator='all_cells', by='treatment')
# → Shows proportions and confidence intervals for each group

# 3. Verify sample sizes are adequate
assert len(data) >= 100, "Sample too small for meaningful proportions"
```

## When to Use Each Chart

| Situation | Chart Type | Why |
|-----------|-----------|-----|
| Comparing counts across independent categories (no nesting) | Count bar chart | Proportional ink; easy to read |
| Comparing proportions (parts of a whole) | Stacked bar or proportion bar | Shows both parts and whole |
| Showing group sizes matter for interpretation | Mosaic plot | Visualizes sample composition |
| Many categories (>10) | Grouped bars or faceted bars | Avoids clutter; easier to read |
| Emphasizing that groups are NOT all the same size | Mosaic plot | Width variation shows imbalance |
| Want to minimize visual emphasis on axis baseline | Dot chart (not bar) | No proportional ink requirement |

## Example Workflow

```python
import matplotlib.pyplot as plt
from utilities.data_inspection import sample_composition
from utilities.color_scales import categorical_palette
from utilities.annotations import add_comparison_annotation

# Step 1: Inspect data
sample_composition(data, 'response')
# → responders: 24, non-responders: 24

# Step 2: Calculate abundances
cd8_counts = data.groupby('response')['cd8_positive'].sum()
total_counts = data.groupby('response').size()
proportions = (cd8_counts / total_counts * 100).round(1)

# Step 3: Create bar chart
fig, ax = plt.subplots(figsize=(4, 5))
colors = categorical_palette(n_categories=2)
bars = ax.bar(proportions.index, proportions.values, color=colors, width=0.6)

# Step 4: Annotate
ax.set_ylabel("CD8+ T cells (%)")
ax.set_xlabel("Response status")
ax.set_ylim(0, 100)
add_comparison_annotation(ax, groups=['responder', 'non-responder'], p_value=0.001, test_name="Fisher's exact")

# Step 5: Add sample sizes
for i, (group, prop) in enumerate(proportions.items()):
    n = total_counts[group]
    ax.text(i, prop + 2, f"n={n}", ha='center', fontsize=9)

fig.savefig("cd8_by_response.png", dpi=300, bbox_inches='tight')
```

---

See `bar_charts.py` and `proportion_and_mosaic.py` for utility functions.
