# Stage 4 — Assumptions and the Three Routes

> Three ways to make a probability statement: parametric, non-parametric and
> permutation. Don't treat a normality check as a mechanical pass/fail rule.

Run all three routes on the same data and report all three. Agreement makes the
claim robust to the modelling choice. Disagreement localises which assumption the
claim rests on — and *that* is the finding.

## The three routes, and what each buys

| Route | Mechanism | Assumes | Costs |
|-------|-----------|---------|-------|
| Parametric (t) | model: measurement = group mean + normal noise; t = difference / SE | a distribution family with fixed parameters; normal means | wrong family → wrong p |
| Non-parametric (ranks, signs) | convert to signs or ranks and reason with exact probability | independence + the test-specific null (symmetry for signed-rank; nothing about shape for the sign test) | discards magnitude; loses effect size in original units; ~5% efficiency when data really are normal |
| Permutation | act out the null story: shuffle labels, recompute, compare | exchangeability — the *entire* distribution is identical under the null, not just the mean | any statistic is allowed, and different statistics give different answers |

Rank tests lose little when the data are normal and can need less data when
tails are heavy. Permutation is the most flexible and the most demanding: its
null is a statement about whole distributions, so a difference in variance alone
can reject it.

## Checklist

- [ ] **Shape is described, not graded.** Q-Q plots per group, or on the paired
      differences for paired data. Report skewness, excess kurtosis, and the
      Shapiro-Wilk p-value where it is defined — labelled as a test whose power
      depends on n.
- [ ] **The n is taken into account when reading the diagnostic.** Below n ≈ 30
      the t-distribution's heavier tails matter and a normality test has almost
      no power, so a non-significant result is not evidence of normality. Above
      it, averages are close to normal by the central limit theorem, so a
      significant normality test need not invalidate a test about means.
- [ ] **For paired data, the diagnostics are on the differences,** which is
      where the assumption lives.
- [ ] **Variance is compared,** with Levene (median-centred) rather than
      Bartlett where normality is in doubt.
- [ ] **All three routes are reported,** including the ones that disagree.
- [ ] **The permutation statistic is named.** "A permutation test" is not a
      method description; "label permutation on the difference in means,
      two-sided, 10,000 permutations, seed 0" is.
- [ ] **Independence is flagged as uncheckable** — no diagnostic here can
      establish it, and neither can any amount of data.

## What to compute

```python
from utilities.resampling import three_routes
from utilities.assumption_checks import shape_table, variance_comparison, paired_difference_check

shape_table({"treated": a, "control": b})
variance_comparison(a, b)

routes = three_routes(a, b, n_perm=10000)          # or paired=True
routes[["route", "test", "p_value", "assumes", "min_achievable_p"]]
routes.attrs["unanimous"]                           # do they all agree at 0.05?
```

For paired designs:

```python
paired_difference_check(after, before)
# ... 'symmetry_concern': True,
#     'signed_rank_note': 'differences are visibly asymmetric (|skew| > 1): ...',
#     'direction_vs_magnitude_note': 'direction and magnitude disagree: ...'
```

## Figure

**F3, distribution shape** — one Q-Q panel per group (or on the differences),
each carrying n, skewness, excess kurtosis and the Shapiro p-value, with the
figure-level note stating that the diagnostic is descriptive. The figure title
says so explicitly, because a bare Q-Q plot in a review reads as an accusation.

**The per-claim panel** carries the three routes as a dot plot on a log p-axis,
with the unreachable region shaded where a route's floor exceeds alpha.

## Findings template

> **Shape.** *[Group]* is *[symmetric / right-skewed / heavy-tailed]* (n = *[n]*,
> skew = *[s]*, Shapiro p = *[p]*). At this n, *[what that does and does not
> imply]*.
>
> **Routes.** Welch's t gives *[p]*, Mann-Whitney *[p]*, label permutation on the
> difference in means *[p]* (10,000 permutations, seed 0). *[All three agree at
> alpha, so the claim does not depend on the modelling choice. / They disagree:
> the claim holds under [route] and not under [route], which means it rests on
> [the specific assumption].]*

## Verdict implications

| Observation | Verdict |
|-------------|---------|
| All routes agree at alpha | `supported` on this axis |
| Routes disagree at alpha | `assumption-dependent` — name which assumption |
| Signed-rank used on visibly asymmetric differences | `supported-with-caveat`; report sign test alongside |
| Permutation reported without naming the statistic or the shuffling scheme | note; the tested null is unknown |

## What not to write

Never "normality was violated, so the test is invalid". Write what the shape is,
what it implies at this n, and whether the routes disagree. A t-test on skewed
data with n = 200 is usually fine; a t-test on 5 observations is fragile whatever
the Shapiro p-value says.
