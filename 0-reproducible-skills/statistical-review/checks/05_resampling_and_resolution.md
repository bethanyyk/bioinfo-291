# Stage 5 — Resampling and Resolution

> Small samples limit the possible permutations. With n = 3 vs 3, the minimum p
> is 0.10 even when the t-test gives p = 0.0004.

Before interpreting any small p-value, ask whether it could have been small. A
design has a hard floor set by how many relabellings exist, and a procedure has a
resolution set by how many iterations were run. Both are computable from the
design alone, before looking at any effect.

## The floors

| Design | Floor | Function |
|--------|-------|----------|
| Two-sample permutation, n1 vs n2 | 2 / C(n1+n2, n1) two-sided | `min_achievable_p(n1, n2)` |
| 3 vs 3 | **0.10** | `min_achievable_p(3, 3)` |
| 4 vs 4 | 0.029 | `min_achievable_p(4, 4)` |
| Paired sign-flip, n pairs | 2 / 2**n | `min_achievable_p(n, paired=True)` |
| Sign test, 5 non-zero pairs | **0.0625** | `min_sign_test_p(5)` |
| Any permutation with N iterations | 1 / (N+1) | `p_value_resolution(N)` |

A p-value reported below its design's floor did not come from the data's own
combinatorics; it came from the parametric model. That is not automatically
wrong — the model may be right — but it must be stated, and the review reports
it as `underpowered-resolution` when the claim's whole weight rests on it.

## Checklist

### Resolution

- [ ] **The design's floor is above or below alpha.** Compute it before anything
      else; a design that cannot reach 0.05 cannot support a claim of
      significance at 0.05, whatever test is applied.
- [ ] **The reported permutation p is achievable** at the stated number of
      permutations: p = (k+1)/(N+1), so p < 1e-6 needs about 1e6 permutations. A
      paper reporting p < 1e-8 from 10,000 permutations reports something its
      procedure cannot produce.
- [ ] **The number of permutations is stated,** along with the seed.
- [ ] **Where many features were tested, an adaptive procedure is used** rather
      than a fixed small N for every feature.

### How the permutation was done

- [ ] **The shuffling scheme is stated, and it matches the null.** Shuffling
      across the whole dataset encodes "the groups are identical"; shuffling
      within patient encodes "the label is arbitrary within patient"; for a
      spatial expression gradient, shuffling within cell type. These are
      different hypotheses.
- [ ] **Paired data are permuted within pairs,** not across the whole sample.
- [ ] **The statistic is named** — permutation allows any statistic, and
      different statistics give different answers on the same data. That freedom
      is the method's cost, and a review must know which one was used.

### Bootstrap and intervals

- [ ] **The bootstrap's purpose is stated**: estimating the distribution of a
      summary statistic, typically to get a confidence interval. It does not test
      a hypothesis.
- [ ] **The key assumption is acknowledged**: the sample is representative. The
      bootstrap propagates a biased sample's bias; it does not repair it.
- [ ] **Enough iterations for the tails.** Interval endpoints need far more
      resamples than the centre; below ~2000 the endpoints are unstable.
- [ ] **Bias is considered** where the statistic is known to be biased under
      resampling (variance estimates, for instance), and a bias-corrected method
      used — `bootstrap_ci(..., method="BCa")` is the default here.

## What to compute

```python
from utilities.resampling import (permutation_test, bootstrap_ci,
                                  min_achievable_p, min_sign_test_p, p_value_resolution)

min_achievable_p(3, 3)          # 0.10 — before looking at the data at all
perm = permutation_test(a, b, statistic="mean_diff", n_perm=10000, random_state=0)
perm["p_value"], perm["min_achievable_p"], perm["floor_binds"], perm["p_resolution"]

# a stratified null, where the paper's design demands one
permutation_test(a, b, strata=cell_type_labels, n_perm=10000)["stratified"]   # True

boot = bootstrap_ci(a, np.mean, n_boot=10000, method="BCa", random_state=0)
boot["ci_low"], boot["ci_high"], boot["bootstrap_bias"]
```

## Figure

**F5, the null distribution** — the permutation null with the observed statistic
marked, annotated with the p-value, the resolution limit, and the design floor
where it binds (the annotation box turns amber when it does). Beside it, the
bootstrap distribution with its interval, answering the separate question of how
precise the estimate is.

The floor is *drawn*, not merely described: a reader should be able to see that
the p-value could not have been smaller.

## Findings template

> **Resolution.** The design is *[n1]* vs *[n2]*, giving a two-sided permutation
> floor of *[floor]*. The paper reports *p = [p]*. *[This is below the floor, so
> it rests entirely on the parametric model rather than on the data's own
> combinatorics. / The floor does not bind.]*
>
> **Null.** The paper *[states / does not state]* how it permuted. Permuting
> *[scheme]* gives *[p]* over *[N]* permutations (seed *[s]*, resolution
> *[1/(N+1)]*).
>
> **Interval.** The bootstrap *[statistic]* is *[point]*, 95% BCa CI
> *[lo, hi]* over *[N]* resamples.

## Verdict implications

| Observation | Verdict |
|-------------|---------|
| Design floor above alpha | `underpowered-resolution` |
| Reported p finer than the stated permutation count allows | `not-reproducible` on that number |
| Permutation scheme unstated | note; the tested null is unknown |
| Permutation scheme does not match the design (unpaired shuffle of paired data) | `not-reproducible` if the corrected scheme crosses alpha |
| Interval from fewer than ~2000 resamples | note; endpoints are unstable |
