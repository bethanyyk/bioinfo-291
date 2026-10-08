# Universal Visualization Principles

Apply these principles to **every visualization** before publication, regardless of chart type.

## Phase 1: Data Inspection

Before you design anything, understand what you're working with.

### Understand Your Variables

For each variable, ask:

- **What are the units?** (e.g., log₂(TPM), percentage, Kelvin, count)
- **What is the range?** (min, max, IQR; are there outliers?)
- **What fraction is missing?** (And is missingness informative?)
- **What is the distribution?** (Symmetric, skewed, multimodal, bounded?)
- **How many unique values?** (Continuous? Ordinal? Nominal?)

Example:
```python
from utilities.data_inspection import inspect_variable

inspect_variable(data['cd8_expression'])
# → Units: log₂(TPM)
#   Range: -2.1 to 15.8 (IQR: 0.5–8.2)
#   Missing: 0%
#   Distribution: bimodal (many zeros, then right-skewed)
#   Unique values: ~5000 (continuous)
```

### Know Your Sample

- **Total n:** How many observations?
- **Composition:** How are they divided into groups?
- **Imbalance:** Are groups very different sizes?
- **Nesting:** Are observations nested (e.g., cells within patients)?

Example:
```python
sample_composition(data, group_col='treatment')
# → Total n: 16,291 cells
#   Groups: (pre-treatment: 4,200 | on-treatment: 6,100 | post-treatment: 5,991)
#   Imbalance ratio: 1.45
#   Nesting: 48 patients, 200–450 cells per patient
```

### State Your Question Clearly

What comparison are your data meant to support?

- ✓ **Good:** "CD8 abundance differs between responders and non-responders"
- ✗ **Vague:** "Compare CD8 across groups"
- ✓ **Good:** "Is CD4 expression correlated with patient age?"
- ✗ **Vague:** "Scatter plot of CD4 and age"

## Phase 2: Visual Encoding

Now that you understand the data, choose how to represent it.

### Axis Titles: Name the Measurement

- ✓ "CD8+ T cell frequency (%)"
- ✓ "log₂(TPM)"
- ✗ "CD8"
- ✗ "Value"

**Include:** Variable name, units (if not obvious), transformation (if any)

### Captions: Identify the Sample

The caption should tell a reader what subset of data is shown.

- ✓ "CD8+ T cell frequency across 48 melanoma patients, stratified by treatment response (responders n=24, non-responders n=24)"
- ✓ "Correlation matrix of immune markers in 10,000 treatment-naive cells"
- ✗ "CD8 by response"

**Include:** What variables, which sample, sample size(s), filters or subsets applied

### Color: Match the Variable Structure

Choose a color encoding that matches what you're encoding:

#### Categorical (Nominal) Data
**Use:** Distinct hues, no implied order
- Good: ["#E41A1C", "#377EB8", "#4DAF4A"] (red, blue, green)
- Bad: ["light gray", "medium gray", "dark gray"] (implies order)
- Bad: Single hue at different saturations (implies magnitude)

```python
from utilities.color_scales import categorical_palette
colors = categorical_palette(n_categories=5, palette="default")
```

#### Ordered/Magnitude Data (Ordinal, Continuous)
**Use:** Sequential lightness (darker = higher value)
- Good: Sequential from light to dark
- Bad: Rainbow palette (implies categorical)
- Bad: Unordered hues

```python
colors = diverging_palette(vmin=data.min(), vmax=data.max(), 
                           palette="RdBu_r")
```

#### Diverging Data (Has a Meaningful Center)
**Use:** Diverging from a central color to two extremes
- Good: Light center, two colors diverging (e.g., blue ← gray → red)
- Bad: Single color (loses the center)
- Bad: Sequential scale (doesn't emphasize the center)

```python
colors = diverging_palette(vmin=-1, vmax=1, center=0, 
                           palette="RdBu_r")
```

### Non-Color Channels Preserve Meaning

When hue alone is insufficient (or for accessibility), use:

- **Line patterns (dashes, dots, solids):** Distinguish time series
- **Marker shapes (circles, triangles, squares):** Mark subgroups
- **Facets (separate panels):** Show categories, ordered or unordered
- **Labels and annotations:** Name points, highlight comparisons
- **Opacity:** Show density or confidence
- **Size:** Encode a fourth variable (but see scatter-plot section on bubble sizing)

Example: Time series with observed vs. estimated values
```
—— solid line: observed values
- - dashed line: imputed values (missing data)
... dotted line: forecast/estimated
```

### Annotations State the Comparison Directly

Good annotations:
- Name what is being compared
- State the outcome or finding
- Include units if not obvious in axes
- Retain context (show n, effect size, p-value as appropriate)

Examples:
```
✓ "p < 0.001 (Welch's t-test)"
✓ "Responders (n=24) vs. non-responders (n=24)"
✓ "Fold change: 2.3× (95% CI: 1.8–2.9)"
✗ "Different" (outcome not stated)
✗ "p < 0.05" (no context; which test?)
```

Use `utilities.annotations.add_comparison_annotation()` to streamline this.

## Phase 3: Chart Selection

Given your question and data, which chart type answers it best?

### How Many? What Fraction?
→ See `abundance_and_proportion/`

**Key principle: Proportional ink**  
When a shaded region represents a value, the area must be proportional to that value.

- **Problem:** Bar chart with truncated y-axis (area ≠ value)
- **Solution:** Start y-axis at 0, or use dot chart instead

**Denominator matters:**
- "30% of what?" (total sample, subgroup, treatment?)
- Always state the denominator in caption or on axis

### How Does It Vary?
→ See `distributions/`

**Key principle: Show the data, not just summaries**  
Identical means can hide different distributions.

- **Problem:** Reporting only mean ± SD
- **Solution:** Show points, histogram, or box + violin

**Bin width and smoothing matter:**
- Histogram bin width can reveal or hide modes
- Violin smoothing choice affects shape
- Specify these in methods or caption

### How Do Variables Relate?
→ See `relationships/`

**Key principle: Overplotting hides density**  
When many points overlap, you lose information about how many.

- **Problem:** 10,000 cells plotted as overlapping points (looks sparse)
- **Solution:** Use transparency, hexbin, or 2D density heatmap

**Scale choice serves the question:**
- Linear scale shows absolute differences
- Log scale shows relative differences
- Choose the scale that answers your question while retaining context

### How Does It Change?
→ See `time_and_change/`

**Key principle: A line means connection**  
Lines imply order or causation. Don't connect unrelated points.

- **Problem:** Bar plot with line connecting groups (misleads)
- **Solution:** Only use lines for time series, paired measurements, or trajectories

**Missing and estimated values need different treatment:**
- Observed: solid line
- Missing: gap or break in line
- Estimated/imputed: dashed or dotted line

### What Pattern Spans Many Features?
→ See `patterns_across_features/`

**Key principle: Shared scales make comparisons**  
If features are in separate panels, they must use the same scale to be comparable.

- **Problem:** Each panel autoscales (hides differences)
- **Solution:** Share scales across all panels, document if ranges differ

**Color scale choice for heatmaps:**
- Sequential lightness for magnitude
- Diverging for data with a center
- Specify transform (log, square root) if axes are transformed

## Phase 4: Publication Checklist

Before finalizing any figure, verify:

### Accessibility
- [ ] No text overlaps
- [ ] Labels are horizontal (no hard-to-read angles)
- [ ] Color is not the only way to distinguish categories (use shape, pattern, or label)
- [ ] Sufficient contrast (text vs. background, colors vs. background)
- [ ] Font size ≥ 8pt (check when printed or exported at intended size)

### Clarity
- [ ] Axis titles name the measurement (with units)
- [ ] Caption identifies the sample and comparison
- [ ] Legend or annotation identifies what colors/shapes represent
- [ ] Scale choice (linear/log) is justified or obvious
- [ ] n values shown for each group (especially if sample sizes vary)

### Correctness
- [ ] Axis ranges include all data
- [ ] Missing data is accounted for (shown as gaps, noted in caption)
- [ ] Statistical tests are named (e.g., "Welch's t-test, p < 0.001")
- [ ] Transformations are documented (e.g., "log₂(TPM)")
- [ ] If using custom colors, verify color-blind accessibility (check with simulator)

### Consistency
- [ ] All figures in the same paper use consistent fonts and sizes
- [ ] Same variables use the same colors across figures
- [ ] Same statistical test nomenclature across figures
- [ ] Same caption style and detail level

Use `utilities.publication_checks.verify_figure()` to automate some of these.

## Common Pitfalls

### Pitfall: "My data doesn't fit a single chart type"

**Solution:** Break it into multiple visualizations. Each figure answers one question clearly.

- Instead of: One dense, complex scatterplot with marginal distributions, sized points, and colors
- Try: Separate scatterplot (relationship question), separate histogram or ECDF (distribution question), separate bar chart (abundance question)

**Exception:** Scatterplot with marginal histograms can work well for two-variable exploration, but save it for exploratory analysis, not publication.

### Pitfall: "The y-axis doesn't start at zero"

**Context matters:**
- ✓ Bar charts for magnitudes: y-axis should start at 0 (proportional ink principle)
- ✓ Line plots for time series: y-axis can start at a non-zero value if context demands it (e.g., temperature from 35–40°C is more useful than 0–40°C)
- ✓ Dot charts or slope plots: y-axis can start at a non-zero value (no expectation of proportional ink)

**Rule:** If the reader expects proportional ink (area = value), start at 0. Otherwise, justify in caption.

### Pitfall: "I need to show too much information"

**Solution:** Divide into facets or separate figures.

- Large faceted plot with 20 panels: too much to read at once
- Five separate figures, each answering one specific question: clearer

**Complexity trade-off:** Complexity is acceptable when the viewer *needs* all the information *right now* to answer the question. Otherwise, simplify.

### Pitfall: "My colors are too similar / too distinct"

**Use utility checks:**
```python
from utilities.color_scales import check_color_blind
check_color_blind(my_palette)  # → See how colors appear to color-blind viewers
```

---

## Summary

1. **Inspect your data first** — Understand units, range, distribution, sample composition
2. **State your question clearly** — What comparison are you making?
3. **Use visual encoding that matches your variable types** — Hues for categories, lightness for magnitude
4. **Choose a chart type that answers your question** — See the question-specific modules
5. **Annotate to explain** — Axis titles, captions, and annotations carry meaning
6. **Check for accessibility and correctness** — Use the publication checklist before finalizing

---

**See also:** Each question type has a detailed PRINCIPLE.md file in its subdirectory. Start with your question type.
