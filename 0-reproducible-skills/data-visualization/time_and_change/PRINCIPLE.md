# Time & Change: "How Does It Change?"

When your question is about how a measurement changes over time, between treatments, or across a sequence.

## Core Principle: A Line Means Connection

**Lines imply order, causation, or continuity.** Do not draw lines connecting points unless that connection is meaningful (time series, paired measurements, or trajectories).

### Example (Bad vs. Good)
```
Bad: Bar plot with connecting line across treatment groups
     Implies causation or ordering; but groups are categorical
     
Good: Line plot with time on x-axis
      Line shows natural progression through time
      Meaningful connection between consecutive points
```

## Chart Types for This Question

### Time Series (One or More Lines)
**Best for:** Measurements at regular or irregular time points

**Principles:**
- X-axis: time (date, day, year, etc.)
- Y-axis: measurement of interest
- One line per group or subject (if following individuals)
- **Missing values need explicit treatment:**
  - Solid line: observed values
  - Gap (break in line): missing data in that interval
  - Dashed line: imputed or estimated values
- **Linear vs. logarithmic scale:**
  - Linear: shows absolute change (e.g., burden goes from 100 → 50 cells/mm³)
  - Log: shows relative/fold change (e.g., burden halves, regardless of baseline)

**Checklist:**
- [ ] X-axis is time (with units, e.g., "Days since start of treatment")
- [ ] Time points are equally spaced (or spacing justified)
- [ ] Missing values shown (gaps, not connected)
- [ ] Estimated/imputed values distinguished from observed (dashed line, different color)
- [ ] Legend shows what each line/color represents
- [ ] Y-axis includes 0 (if proportional ink expected) or range is justified

### Paired Lines (Connected Points Before & After)
**Best for:** Comparing individuals across two time points or conditions

**Principles:**
- Two vertical axes (or x-positions) representing "before" and "after"
- Lines connect each individual's before-and-after values
- Line slope shows direction of change
- Color or shape can code for outcome (e.g., responders vs. non-responders)
- Alternative: "spaghetti plot" if many time points

**Checklist:**
- [ ] Each line represents one subject (or observation)
- [ ] Lines are thin (don't dominate the plot) and semi-transparent (show overlaps)
- [ ] Categories (responders vs. non-responders) are color-coded
- [ ] Summary (mean or median) line overlaid, if helpful
- [ ] Sample size shown (n = number of lines)

### Change Plot (Slopes, Differences)
**Best for:** Emphasizing the direction and magnitude of change

**Principles:**
- X-axis: categories or time points
- Y-axis: the change (e.g., percentage point difference, fold-change)
- Can use slopes (one line per subject), bars, or points
- Useful when baseline varies widely (focus on change, not absolute value)

**Checklist:**
- [ ] Scale is change/difference (not absolute values)
- [ ] Zero line is shown and emphasized (changes above/below are clear)
- [ ] Direction of change is interpretable (which direction is "good"?)
- [ ] Units of change are labeled

## Common Pitfalls

### Pitfall 1: Lines Connecting Unrelated Points
**Problem:** Bar plot with lines across treatment groups (categorical, not ordered).
**Solution:** Use separate bars without lines, or use a line plot only for time series.

```python
# Bad: lines across categorical treatments
ax.plot(['pre', 'on', 'post'], [100, 150, 120])

# Good: separate bars (no connection implied)
ax.bar(['pre', 'on', 'post'], [100, 150, 120])
```

### Pitfall 2: Missing Data Obscured
**Problem:** Time series with gaps in observations, but line still connects (implies continuity).
**Solution:** Show gaps explicitly.

```python
# Bad: connects across missing data
ax.plot(times, values)

# Good: breaks line at missing values
ax.plot(times, values, drawstyle='steps-post')  # or use explicit NaN handling
```

### Pitfall 3: Spaghetti Plot Too Dense
**Problem:** Many individual trajectories plotted; center of the mass is invisible.
**Solution:** Add summary line (mean or median) on top; or use facets; or use 2D density.

### Pitfall 4: Linear Scale Obscures Relative Change
**Problem:** Baseline values differ widely; absolute change plot is dominated by high-baseline subjects.
**Solution:** Use log scale (to show fold-change) or fold-change on y-axis instead of absolute value.

Example: Drug response from baseline to week 4:
- Patient A: burden 1 million → 10,000 (100-fold reduction)
- Patient B: burden 1,000 → 100 (10-fold reduction)
- Linear plot: A dominates; hard to see that both improved
- Log plot: Both show as parallel downward slopes; equal fold-change visible

## Data Inspection Checklist

Before creating any time/change visualization:

```python
from utilities.data_inspection import inspect_variable

# 1. Understand time structure
# Are measurements at regular intervals? Irregular?
# How many time points per subject?

# 2. Check for missing time points
# Are there gaps? Are gaps informative (e.g., subjects dropped out)?
missing_counts = data.groupby('subject_id')['measurement'].apply(lambda x: x.isna().sum())
print(f"Missing measurements per subject: {missing_counts.describe()}")

# 3. Check baseline variation
baseline = data[data['timepoint'] == 'baseline']['measurement']
cv = baseline.std() / baseline.mean()  # coefficient of variation
if cv > 0.5:
    print("High baseline variation; consider fold-change visualization")
```

## When to Use Each Chart

| Situation | Chart Type | Why |
|-----------|-----------|-----|
| Regular time series, one group | Single line | Shows trend clearly |
| Multiple groups over time | Multiple lines (color-coded) | Shows whether trends differ by group |
| Following individuals over time | Spaghetti plot (thin, transparent lines) | Shows heterogeneity in trajectories |
| Before-after comparison, paired data | Paired lines or dot plot | Shows who changed, direction |
| Multiple time points, want summary | Line + confidence band | Shows trend and uncertainty |
| Large baseline variation, relative change | Log scale or fold-change y-axis | Shows fold-change (not absolute) |
| Many irregular time points | Loess smooth or GAM curve | Shows underlying trend without overemphasizing individual points |

## Example Workflow

```python
import matplotlib.pyplot as plt
import pandas as pd
from utilities.color_scales import categorical_palette

# Step 1: Prepare data
# Assume: data has columns [subject_id, timepoint, measurement, response]
data_sorted = data.sort_values(['subject_id', 'timepoint'])

# Step 2: Create time series (multiple groups)
fig, ax = plt.subplots(figsize=(8, 6))
colors = categorical_palette(n_categories=2)

for i, (response, color) in enumerate(zip(['responder', 'non-responder'], colors)):
    subset = data_sorted[data_sorted['response'] == response]
    
    # Plot each subject (thin, transparent lines)
    for subject in subset['subject_id'].unique():
        subj_data = subset[subset['subject_id'] == subject]
        ax.plot(subj_data['timepoint'], subj_data['measurement'], 
               color=color, alpha=0.2, linewidth=1)
    
    # Overlay group mean
    mean_by_time = subset.groupby('timepoint')['measurement'].mean()
    ax.plot(mean_by_time.index, mean_by_time.values, 
           color=color, linewidth=2.5, label=response, marker='o')

ax.set_xlabel("Timepoint")
ax.set_ylabel("Tumor burden (cells/mm³)")
ax.legend()
ax.set_ylim(bottom=0)

fig.savefig("tumor_burden_over_time.png", dpi=300, bbox_inches='tight')
```

---

See `time_series.py` and `paired_and_change.py` for utility functions.
