# Stage 3 — Test Selection

> Key test terms: one-sample, two-sample or paired; one-tailed or two-tailed;
> equal or unequal variance.

Check that the named test matches the design on each of those axes, then check
the harder thing: that the test's null hypothesis is the one the scientific
question needs.

## Checklist

### Does the test match the design

- [ ] **One-sample, two-sample or paired.** Paired data analysed as two
      independent samples throws away the pairing and usually loses power;
      independent data analysed as paired is invalid.
- [ ] **One-tailed or two-tailed.** A one-tailed test must be justified by a
      directional hypothesis stated *before* the data were seen. Halving a
      two-tailed p after the fact is p-hacking with extra steps. Where a paper
      reports a one-tailed p without pre-registration, report the two-tailed
      value alongside.
- [ ] **Equal or unequal variance.** Welch's is the safe default; it costs
      negligible power when variances are equal. Student's equal-variance t with
      a variance ratio above ~2 is not appropriate.
- [ ] **The test handles the number of groups.** Three groups compared by three
      pairwise t-tests is three tests, and stage 6 will count them.
- [ ] **Repeated measures are modelled, not ignored** — see stage 2.

### Does the null match the question

- [ ] **Mann-Whitney is not a test of equal medians.** It tests which group
      tends to be larger across cross-group pairs. Two groups can share a median
      and still differ strongly by this test. A paper claiming "medians did not
      differ (Mann-Whitney p = …)" has mismatched its null to its wording.
- [ ] **Wilcoxon signed-rank assumes symmetric differences.** A zero median is
      not sufficient. Where the differences are visibly skewed, the sign test or
      a permutation test on a stated statistic is the honest route.
- [ ] **A test of means does not answer "do most subjects improve".** The
      scientific question decides the summary statistic. Where 9 of 10 subjects
      improve and the mean change is negative, a t-test and a sign test answer
      two different, both legitimate, questions — and the paper must say which it
      asked.
- [ ] **"Parametric" does not mean "Gaussian".** It means a distribution family
      with a fixed set of parameters. A paper that justifies a Poisson or
      negative-binomial model has not thereby assumed normality.
- [ ] **"Non-parametric" does not mean assumption-free.** Independence and the
      test-specific null still apply.

## What to compute

```python
paper["inventory"]["named_tests"]    # which tests the paper says it used, with context
paper["inventory"]["tests"]          # the statistics it reports, with df

from utilities.assumption_checks import variance_comparison
variance_comparison(a, b)["recommended_t_variant"]
```

Cross-check the named test against the degrees of freedom: a two-sample t on
n1 = n2 = 10 has df = 18 under equal variance and a non-integer df under Welch.
Integer df where the paper claims Welch, or vice versa, tells you which test was
actually run.

## Findings template

> **Test and design.** The paper analyses *[design]* with *[named test]*. The
> design is *[paired / independent / one-sample]* because *[reason]*, so
> *[appropriate test]* is the matching choice. Re-running with *[appropriate
> test]* gives *[p]* against the reported *[p]*.
>
> **Null and question.** The claim is *[wording]*, which is a claim about
> *[means / medians / stochastic ordering / direction of change]*. *[Named test]*
> tests *[its actual null]*. *[These match / The claim should be reworded, or a
> test of [X] reported instead.]*

## Verdict implications

| Observation | Verdict |
|-------------|---------|
| Wrong pairing structure, and the corrected test crosses alpha | `not-reproducible` |
| Equal-variance t with variance ratio > 2, corrected test still significant | `supported-with-caveat` |
| One-tailed without pre-stated direction, two-tailed crosses alpha | `not-reproducible` |
| Null does not match the claim's wording | `supported-with-caveat`; propose the wording or the test |
| Test cannot address the claim at all | `not-assessable`; state what analysis would |
