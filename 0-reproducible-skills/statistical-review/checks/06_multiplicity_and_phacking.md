# Stage 6 — Multiplicity and Analytic Flexibility

> A p-value only says how likely data like yours are under certain assumptions.
> If you run 100 tests, about 5 will come out "significant" by chance.

The arithmetic first: at alpha = 0.05 over m tests you expect 0.05 × m
significant results even when every null is true. Put that number next to the
paper's hit count before anything else.

## The corrections answer different questions

| Correction | Controls | Data-dependent? | Effect |
|------------|----------|-----------------|--------|
| None | the per-test error rate only | — | 0.05 × m expected false positives |
| Bonferroni | FWER — the probability of **any** false positive | No | very conservative; threshold alpha/m regardless of the results |
| Benjamini-Hochberg | FDR — the expected **fraction** of false positives among the calls | Yes | more permissive; the effective threshold rises with the number of small p-values, so up to ~5% of the calls may be false |
| Max-statistic permutation | FWER, respecting dependence between tests | Yes | the right choice for correlated features; Bonferroni is both conservative and wrong about the dependence there |

Neither Bonferroni nor BH is "correct" in general. They answer different
questions, and the report states which question the paper's claim needs.

## Checklist

- [ ] **Every test is counted,** including those visible only as figure panels,
      supplementary tables, or subgroup analyses. The denominator for a
      correction is the number of tests *run*, not the number reported.
- [ ] **0.05 × m is stated next to the hit count.**
- [ ] **A correction is applied where the claim is "these features differ",**
      as distinct from one pre-specified primary comparison.
- [ ] **The correction matches the dependence structure.** Bonferroni across
      thousands of correlated neighbouring features is the wrong tool; the
      max-statistic permutation approach carries the correlation in the
      permutation itself.
- [ ] **The paper says which correction it used, over how many tests.**
- [ ] **Subgroup and post-hoc analyses are identified as such.**
- [ ] **The p-value distribution is inspected** for shape: a left spike means
      real signal, a flat profile means mostly nulls, and a spike just under the
      threshold means something about the analysis rather than the biology.

## What to compute

```python
from kernel import audit_reported_statistics
audit = audit_reported_statistics(paper["inventory"])

audit["multiplicity"]["summary"]
#            correction                                     controls  n_significant  threshold_on_raw_p
#                  none                nothing (per-test error rate)             14               0.050
#            Bonferroni      FWER: probability of any false positive              2               0.001
#    Benjamini-Hochberg   FDR: expected fraction of false calls <= 0.05           7               0.012

audit["multiplicity"]["note"]
# '50 tests at alpha=0.05 yield about 2.5 significant results by chance alone. ...'

audit["just_below"]
# {'n_in_suspect_window': 9, 'n_in_reference_window': 2, 'binomial_p_excess': 0.033,
#  'interpretable': True, ...}
```

For correlated features, where the paper supplies or you can recompute a
permutation null per feature:

```python
from utilities.multiplicity import max_statistic_null
max_statistic_null(observed_stats, null_matrix, alternative="two-sided")
```

## Figure

**F6, multiplicity** — the p-value histogram against the uniform expectation,
with the just-below-threshold window highlighted; beside it, the hit count under
no correction, BH and Bonferroni, with 0.05 × m drawn as a reference line. The
figure shows at a glance how much of the paper's hit list depends on no
correction having been applied.

## The p-hacking diagnostic, and its limits

Published p-values cluster just below 0.05 across the literature, and that
clustering is measurable. `just_below_rate` compares the count in [0.045, 0.05)
with the equally wide [0.04, 0.045). Under any continuous distribution — null or
not — the two should be similar, and the left window slightly larger where there
is real signal.

Three limits the report must carry every time this is used:

1. It is a property of a **collection** of p-values and says nothing about any
   individual result.
2. It needs a reasonable number of p-values. Below 20 the function returns
   `interpretable: False`, and the report prints the counts without a conclusion.
3. Rounding conventions alone can produce apparent pile-up at exactly 0.05.

And the thing it cannot do: **it cannot distinguish p-hacking from legitimate
exploratory analysis.** The line between them is genuinely fine. Report the
counts and what they are consistent with. Never attribute intent, and never
write "the authors p-hacked" — write that the distribution of reported p-values
shows *[counts]*, which is *[consistent with / not distinguishable from]*
*[explanation]*, and that a pre-registered analysis plan or the full set of
tests run would resolve it.

## Findings template

> **Multiplicity.** The paper reports *[m]* tests at alpha = 0.05, where about
> *[0.05m]* significant results are expected by chance alone. *[k]* are reported
> as significant uncorrected; *[k_bh]* survive Benjamini-Hochberg and *[k_bonf]*
> survive Bonferroni. The claims that do not survive correction are *[ids]*.
>
> **Correction choice.** The paper used *[correction]* over *[stated m]*. *[The
> features tested are correlated, so a max-statistic permutation correction
> would be less conservative and would respect the dependence. / This matches the
> claim being made.]*
>
> **Distribution of reported p-values.** *[n]* fall in [0.045, 0.05) against
> *[n]* in [0.04, 0.045) (one-sided binomial p = *[p]*). *[Interpretation with
> the three limits stated. / Fewer than 20 p-values: not interpretable at this
> scale.]*

## Verdict implications

| Observation | Verdict |
|-------------|---------|
| Claim does not survive a correction the analysis needs | `supported-with-caveat`, or `not-reproducible` where the paper asserted significance after correction |
| No correction applied across an explicitly exploratory screen | `supported-with-caveat`; report the BH recount |
| Correction applied over fewer tests than were run | report the corrected count at the true m |
| Bonferroni over correlated features | note; propose max-statistic permutation |
