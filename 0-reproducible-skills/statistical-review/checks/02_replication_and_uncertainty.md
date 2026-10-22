# Stage 2 — Replication and Uncertainty

> SD describes the spread of individual observations. SE/SEM describes the
> uncertainty in the summary statistic, and it shrinks as n grows.

Two questions: **what is n**, and **what do the reported error bars describe**.
Both are answerable from the data plus one sentence of the methods, and both are
common sources of a p-value that is too small by a large factor.

## Checklist

### What is n

- [ ] **The paper states the unit of replication,** and the noun matters more
      than the number: "n = 3 independent experiments" and "n = 3000 cells from
      3 mice" can describe the same experiment, and only the first names the unit
      a test may treat as independent.
- [ ] **The n used by the test is the number of independently randomized
      units,** not the number of measurements.
- [ ] **Degrees of freedom are consistent with the stated n.** A paper reporting
      n = 3 and t(18) ran its test on 20 measurements, not on 3 units. This is
      checkable from the PDF alone and is one of the highest-yield checks in the
      whole review.
- [ ] **Technical and biological replicates are distinguished.**

### What the error bars describe

- [ ] **The figure legend says what the error bars are.** If it does not, the
      figure cannot be read, and that is the finding.
- [ ] **The stated error matches what was plotted.** SD and SEM differ by
      sqrt(n), so at n = 25 an SEM bar is one-fifth the length of an SD bar.
- [ ] **The error bar matches the claim being made.** A claim about the spread of
      a population needs SD; a claim about the precision of a mean needs SEM or a
      CI. A paper describing "variability between animals" with SEM bars is
      making a confidence statement about its own mean.
- [ ] **Reported means are attainable** for the stated n, where the measurement
      is integer-valued (GRIM).

## What to compute

```python
from kernel import describe_data
desc = describe_data(data, value_col="volume", group_col="arm", unit_col="mouse")

desc["independence"]
# {'n_observations': 300, 'n_units': 6, 'obs_per_unit_mean': 50.0,
#  'icc': 0.83, 'design_effect': 41.7, 'effective_n': 7.2,
#  'pseudoreplication_risk': True, 'note': '...'}
```

The design effect is `1 + (k-1) * ICC`, and the effective n is the row count
divided by it. When observations cluster within units, a test on rows reports a
p-value far smaller than the design supports.

```python
from utilities.recompute import error_bar_identity, grim_check

error_bar_identity(reported_error=2.0, sd=10.0, n=25)
# best_match='SEM', sd_to_sem_factor=5.0

grim_check(mean=3.15, n=7, decimals=2)
# consistent=False  → 3.15 is not a multiple of 1/7
```

## Figure

**F1, the design map** — observations per unit on the left, independent units per
group on the right, annotated with the ICC and the effective n. The gap between
the row count and the effective n *is* the finding; the figure exists so a reader
sees it without reading a paragraph.

**F2, the raw data** — every observation, with mean ± SD and mean ± SEM drawn
side by side on the same axis. This is the clearest way to show what a
bar-and-SEM plot was hiding.

## Findings template

> **Unit of replication.** The paper reports *n = [N] [unit]* for *[claim]*. The
> supplied data contain *[rows]* observations from *[units]* units (*[k]* per
> unit, ICC = *[icc]*), giving an effective n of about *[eff]*. The reported
> *[statistic]* with *[df]* degrees of freedom is consistent with a test run on
> *[rows]* rather than on *[units]*.
>
> **Resolution.** Re-running the comparison at the unit level gives *[p]* against
> the reported *[p]*.

## Verdict implications

| Observation | Verdict |
|-------------|---------|
| df consistent with measurements, not units; unit-level test crosses alpha | `not-reproducible` |
| Clustering present but the unit-level test still crosses alpha | `supported-with-caveat` |
| Error bars are SEM while the claim is about spread | `supported-with-caveat` |
| Error bar definition absent from the paper | note; qualifies any figure-based claim |
| GRIM inconsistency | report the arithmetic; ask which n or which exclusions apply |
