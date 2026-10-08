# Data Visualization Guidelines & Tools

A structured approach to choosing and designing scientific visualizations. Organized by the question your data answers, not just chart type.

## Organizing Principle: Start with the Question

Every visualization answers one of five types of questions:

| Question | Useful starting point | Key check | Module |
|----------|----------------------|-----------|--------|
| **How many? What fraction?** | Count bars, proportion bars, mosaic | Denominator and sample composition | `abundance_and_proportion/` |
| **How does it vary?** | Points, histogram, box, violin, ECDF | Raw n, bins or smoothing | `distributions/` |
| **How do variables relate?** | Scatter, binned counts, correlation matrix | Subgroups, outliers, units | `relationships/` |
| **How does it change?** | Time series, paired lines, change plot | Order, pairing, missing intervals | `time_and_change/` |
| **What pattern spans many features?** | Heatmap, small multiples | Transform, color scale, shared axes | `patterns_across_features/` |

**Start here:** Identify your question, then read the corresponding module.

## Universal Principles (All Charts)

Before choosing a specific chart type, apply these principles to every visualization:

### Data Inspection
1. **Understand your variables:** units, range, missingness, distribution
2. **Know your sample:** composition, groups, potential imbalances
3. **Identify the comparison:** What are you comparing? Can your data support it?

### Visual Encoding
1. **Axis titles** should name the measurement (not just "value" or "count")
2. **Captions** should identify the sample (what subset is shown?)
3. **Color** should match the structure of the variable:
   - **Hue** distinguishes categories (nominal)
   - **Ordered lightness** communicates magnitude (ordinal or continuous)
   - **Diverging scale** needs a meaningful center (diverging continuous)
4. **Non-color channels** (labels, shapes, line patterns, panels) preserve meaning without relying on hue alone
5. **Annotations** state the comparison directly — identify outcome, units, intended comparison, retain context

### Design Decisions
- Prefer **simpler, more direct visualizations** when the question and data allow
- **Complexity is a tradeoff:** between information available and what a viewer can read
- Some complex figures (e.g., scatterplot with marginal histograms) can work if you choose the right moment to show them
- **Iterate and refine:** Start with exploration, then move to publication-grade presentation

## Directory Structure

```
data-visualization/
├── README.md                          ← This file
├── universal_principles.md            ← Detailed universal principles
├── abundance_and_proportion/
│   ├── PRINCIPLE.md                   ← Principles for count & proportion questions
│   ├── bar_charts.py                  ← Bar chart utilities and checks
│   ├── proportion_and_mosaic.py       ← Proportion bars and mosaics
│   └── examples.md                    ← Example checklist and common pitfalls
├── distributions/
│   ├── PRINCIPLE.md
│   ├── exploratory.py                 ← Points, histograms, ECDF (exploration)
│   ├── summary_shapes.py              ← Box plots, violin plots (summaries)
│   └── examples.md
├── relationships/
│   ├── PRINCIPLE.md
│   ├── scatter_and_binned.py          ← Scatter, hexbin, binned heatmaps
│   ├── correlation_matrices.py        ← Correlation visualization
│   └── examples.md
├── time_and_change/
│   ├── PRINCIPLE.md
│   ├── time_series.py                 ← Time series and temporal trends
│   ├── paired_and_change.py           ← Paired lines and change plots
│   └── examples.md
├── patterns_across_features/
│   ├── PRINCIPLE.md
│   ├── heatmap.py                     ← Heatmaps and clustering
│   ├── small_multiples.py             ← Faceted plots, small multiples
│   └── examples.md
└── utilities/
    ├── __init__.py
    ├── color_scales.py                ← Color scales for different data types
    ├── data_inspection.py             ← Utilities to inspect data before plotting
    ├── annotations.py                 ← Utilities for annotations and labels
    └── publication_checks.py          ← Checklist for publication-ready figures
```

## Workflow

### 1. Start with Data Inspection

Before you design any visualization, use utilities to understand your data:

```python
from utilities.data_inspection import inspect_variable, sample_composition

inspect_variable(data['age'])        # → units, range, missingness, distribution
sample_composition(data, group_col)  # → group sizes, imbalances
```

### 2. Identify Your Question

Which of the five questions does your data answer?
- **"How many cells in each group?"** → abundance_and_proportion/
- **"What's the distribution of expression values?"** → distributions/
- **"Are these two genes correlated?"** → relationships/
- **"How does abundance change over treatment?"** → time_and_change/
- **"What cell types express these genes?"** → patterns_across_features/

### 3. Read the Relevant PRINCIPLE.md

Each category has a PRINCIPLE.md file that explains:
- Why certain chart types answer that question
- Common pitfalls for that question type
- Checklist of what to verify before publication

### 4. Choose a Specific Visualization

Use the module for your question type. Each module provides:
- **Utilities:** functions to create the visualization
- **Checks:** verify your visualization follows principles
- **Examples:** worked examples with common variations

### 5. Verify with Checklist

Use `utilities/publication_checks.py` to verify:
- Axes are labeled and units are clear
- Colors are appropriate for the variable type
- Annotations clarify the intended comparison
- No overlapping text or visual conflicts

## Design Philosophy

This toolkit reflects a few core principles:

1. **Form follows function** — The visualization should match your question, not the other way around
2. **Data inspection comes first** — You can't choose a good visualization without understanding your data
3. **Simplicity by default** — Start simple; add complexity only when the data and question demand it
4. **Annotation as explanation** — When a visualization needs explanation, annotations are better than captions
5. **Color is not mandatory** — Many visualizations can work without color or with minimal color

## Examples

### "How many immune cells in each treatment group?"

**Question:** Abundance (How many?)  
**Module:** `abundance_and_proportion/`  
**Starting point:** Bar chart with count, or proportion bar with sample size noted

### "Are CD8 expression values normally distributed?"

**Question:** Distribution (How does it vary?)  
**Module:** `distributions/`  
**Starting point:** Histogram or box plot + violin plot

### "Which genes correlate with patient outcome?"

**Question:** Relationship (How do variables relate?)  
**Module:** `relationships/`  
**Starting point:** Scatter (outcome vs expression) or correlation matrix

### "How does tumor burden change with treatment over time?"

**Question:** Change (How does it change?)  
**Module:** `time_and_change/`  
**Starting point:** Line plot (time on x-axis, burden on y-axis) with observed/missing/estimated treated differently

### "Which cell types are enriched in responders vs non-responders?"

**Question:** Pattern across features (What pattern spans many features?)  
**Module:** `patterns_across_features/`  
**Starting point:** Heatmap (cell types as rows, response groups as columns) or small multiples (one panel per cell type)

## Integration with Your Project

### Copy utilities into your project

```bash
cp -r 0-reproducible-skills/data-visualization/utilities/ your-project/src/
```

### Import in your analysis scripts

```python
from src.utilities.data_inspection import inspect_variable, sample_composition
from src.utilities.color_scales import categorical_palette, diverging_palette
from src.utilities.annotations import add_comparison_annotation
```

### Add to your config

```yaml
# configs/visualization.yaml
figures:
  style: "publication"  # or "exploratory"
  dpi: 300
  palette: "default"    # categorical palette
  diverging_palette: "RdBu_r"
  fontsize: 9
```

## What This Toolkit Is NOT

- **Not a style guide** — These are principles, not strict rules. Apply them with judgment.
- **Not a replacement for domain expertise** — You know your data; use that knowledge.
- **Not a complete graphics library** — Use matplotlib, seaborn, plotly, ggplot2, etc. as your graphics engine; this toolkit documents *principles* and provides *checklists*.

## Next Steps

1. **Read `universal_principles.md`** for detailed design guidance
2. **Choose your question type** and read the corresponding PRINCIPLE.md
3. **Explore the examples.md** in that module for worked cases
4. **Use the utilities** to inspect data and verify your visualization

---

**See also:** The `1-scRNA/` project demonstrates these principles in action. Check `figures/` for publication-ready figures and `src/icb_scrna/plots.py` for how they were implemented.
