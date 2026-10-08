# Patterns Across Features: "What Pattern Spans Many Features?"

When your question involves many variables and you want to show how a pattern affects or manifests across all of them.

## Core Principle: Shared Scales Enable Comparison

**When showing multiple features in separate panels or as columns/rows, they must use the same scale (or differences in scale must be intentional and highlighted).**

### Example
```
Good: Heatmap showing all genes on same color scale
      → Can compare expression levels across genes directly

Bad: Each gene gets its own color scale (auto-scaled per gene)
     → Genes that are actually low can appear high (scaled to their local range)
     → Cannot compare across genes; visual is misleading
```

## Chart Types for This Question

### Heatmap (Matrix)
**Best for:** Many variables × few groups (e.g., genes × cell types, features × samples)

**Principles:**
- Rows: variables (genes, markers, etc.)
- Columns: observations or groups
- Color: represents value (or transformed value, e.g., log fold-change)
- **Color scale choice:**
  - Sequential lightness: for magnitude (0 to max)
  - Diverging: for centered data (e.g., log fold-change from -2 to +2)
- Often includes clustering (reorder rows/columns by similarity)
- Row or column names can be labels on the side; use smaller font

**Checklist:**
- [ ] Color scale is appropriate for data type (sequential or diverging)
- [ ] Color scale is shared across all rows/columns (not auto-scaled per row)
- [ ] Labels (gene names, sample names) are readable
- [ ] Color bar shows what values correspond to what colors
- [ ] Transformation (if any) is documented (e.g., "log₂ fold-change")
- [ ] Clustering method documented (if applicable)

### Small Multiples (Faceted Plots)
**Best for:** Same pattern shown across multiple groups or variables

**Principles:**
- Each panel shows one subset of data (one variable, one group, etc.)
- All panels use the same scale (x and y axes have same range)
- Panels arranged in a grid (rows and columns make sense semantically)
- Allows showing variability across features without overwhelming a single plot

**Checklist:**
- [ ] All panels share x and y axes ranges (or differences are justified)
- [ ] Panel labels clearly identify what each shows
- [ ] Grid layout is logical (don't randomize; order by value or importance)
- [ ] Enough panels to show pattern, but not so many as to require scrolling
- [ ] Legend (if needed) shared across all panels

## Common Pitfalls

### Pitfall 1: Auto-Scaled Per Feature
**Problem:** Each gene in a heatmap gets its own color scale → low-expression genes look high.
**Solution:** Use single color scale across all features.

```python
# Bad: each row auto-scaled
import seaborn as sns
sns.heatmap(data)  # Default behavior varies per row

# Good: shared color scale
vmin, vmax = data.min().min(), data.max().max()
sns.heatmap(data, vmin=vmin, vmax=vmax)
```

### Pitfall 2: Too Many Small Multiples
**Problem:** 50 small panels arranged in a grid; unreadable.
**Solution:**
- Limit to ~16 panels (4×4 grid)
- Or use interactive visualization (hover to see details)
- Or select subset of most important features to show

### Pitfall 3: Different Scales Across Panels
**Problem:** Panel A (0–100) and Panel B (0–10) look equivalent despite 10× difference.
**Solution:** Share scales across panels. If ranges must differ, use log scale or explicitly label each panel.

### Pitfall 4: Clustering Obscures Sample Order
**Problem:** Heatmap reorders samples by clustering; original sample order (e.g., time order) is lost.
**Solution:**
- Cluster features (rows) only, not samples (columns)
- Or, use a "hybrid" approach: cluster within each group, but preserve group order
- Document if ordering is changed from original

## Data Inspection Checklist

Before creating any multi-feature visualization:

```python
from utilities.data_inspection import inspect_variable

# 1. Understand each feature separately
for col in data.columns:
    inspect_variable(data[col])
    # → Range varies widely? May need log scale or different transform

# 2. Check correlations
corr = data.corr()
# High correlation → features are redundant; consider selecting subset

# 3. Assess clustering
from scipy.cluster.hierarchy import dendrogram, linkage
Z = linkage(data.T, method='ward')
# Dendrogram shows which features group together
```

## When to Use Each Chart

| Situation | Chart Type | Why |
|-----------|-----------|-----|
| Many genes/markers × few cell types or samples | Heatmap | Compresses data; easy to scan |
| Same type of plot across many subgroups | Small multiples | Shows whether pattern is consistent |
| Genes ordered by biological importance | Heatmap without clustering | Preserves user-chosen ordering |
| Want to find similar genes/samples | Heatmap with clustering | Dendrograms reveal structure |
| Continuous variable across many features | Small multiples of line plots or scatter | Shows trend consistency |
| Many categorical variables | Heatmap of proportions or Fisher's exact test | Shows which categories associate |

## Example Workflow: Heatmap

```python
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.cluster.hierarchy import linkage
from utilities.color_scales import diverging_palette

# Step 1: Prepare data (e.g., log fold-change matrix)
# Rows: genes
# Columns: cell types
# Values: log₂(fold-change in responders vs. non-responders)

# Step 2: Compute clustering (cluster genes, not cell types)
Z_genes = linkage(data, method='ward')  # Clusters rows (genes)
# Leave columns (cell types) in original order

# Step 3: Create heatmap
fig, ax = plt.subplots(figsize=(8, 12))

# Use diverging palette centered at 0
cmap, norm = diverging_palette(vmin=-2, vmax=2, center=0, palette="RdBu_r")

sns.heatmap(data, 
            cmap=cmap,
            norm=norm,
            xticklabels=True,
            yticklabels=True,  # Gene names
            cbar_kws={'label': 'log₂(fold-change)'},
            ax=ax)

ax.set_title("Immune marker expression\nResponders vs. non-responders")
ax.set_xlabel("Cell type")
ax.set_ylabel("Marker gene")

fig.savefig("immune_markers_heatmap.png", dpi=300, bbox_inches='tight')
```

## Example Workflow: Small Multiples

```python
import matplotlib.pyplot as plt
from utilities.color_scales import categorical_palette

# Step 1: Prepare data subsets
# One plot per cell type, showing expression distribution

fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharex=True, sharey=True)
axes = axes.flatten()

cell_types = data['cell_type'].unique()
colors = categorical_palette(n_categories=2)

# Step 2: Create subplots with shared scales
for idx, cell_type in enumerate(cell_types):
    ax = axes[idx]
    subset = data[data['cell_type'] == cell_type]
    
    for i, response in enumerate(['responder', 'non-responder']):
        resp_subset = subset[subset['response'] == response]
        ax.hist(resp_subset['expression'], alpha=0.5, bins=20, color=colors[i], label=response)
    
    ax.set_title(f"{cell_type} (n={len(subset)})")
    ax.set_ylabel("Frequency")
    ax.set_xlabel("Expression (log₂ TPM)")

# Shared legend
axes[0].legend(loc='upper right')

# Remove extra subplots
for idx in range(len(cell_types), len(axes)):
    fig.delaxes(axes[idx])

fig.suptitle("Gene expression distribution across cell types")
fig.tight_layout()
fig.savefig("expression_small_multiples.png", dpi=300, bbox_inches='tight')
```

---

See `heatmap.py` and `small_multiples.py` for utility functions.
