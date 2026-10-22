# Statistical Review — *[Paper short title]*

*Reviewer: [name/agent] · Date: [YYYY-MM-DD] · Protocol: statistical-review v1*

---

## 1. Inputs and coverage

| Input | Value |
|-------|-------|
| Manuscript | `[path]` — [N] pages, text layer [usable / not usable] |
| Raw data | `[path]` — [rows] rows × [cols] columns, sha256 `[…]` |
| Data corresponds to | [which figures/tables of the paper] |
| Claims extracted automatically | [n] |
| Claims added by hand | [n] |
| Claims with no supporting data | [ids] — reported as `not-assessable` |
| Seeds | `random_state=[s]`, `n_perm=[N]`, `n_boot=[N]` |

> State here, before any finding: what fraction of the paper's claims this review
> actually covers. A review of 3 of 14 claims is a review of 3 claims.

---

## 2. Summary of findings

| Verdict | Claims |
|---------|--------|
| supported | [n] |
| supported-with-caveat | [n] |
| assumption-dependent | [n] |
| underpowered-resolution | [n] |
| not-reproducible | [n] |
| contradicted | [n] |
| not-assessable | [n] |

**At most three findings** — the ones that change how the paper reads. Each one
sentence: the claim, the discrepancy as two numbers, the stage it came from.

1. *[C2: the paper reports p = 0.033 for …; the unit-level test on the 3
   independent experiments gives p = 0.41 (stage 2).]*
2. *[…]*
3. *[…]*

---

## 3. Claim ledger

![F4 — claim ledger](figures/F4_claim_ledger.png)

*Recomputed effect (filled marker) against the reported effect (open marker)
with 95% bootstrap CI; right panel, reported against recomputed p-values on log
axes with alpha = 0.05 marked. Points off the diagonal across 0.05 are
decision-changing.*

| id | claim | verdict | effect rep. | effect recomp. | 95% CI | p rep. | p recomp. | reason |
|----|-------|---------|------------|----------------|--------|--------|-----------|--------|
| C1 | | | | | | | | |

---

## 4. Design and sampling — stages 1–2

![F1 — design map](figures/F1_design_map.png)

*Left: observations per experimental unit. Right: independent units per group,
the n available to a test. Annotated with ICC and effective n.*

![F2 — raw data behind the claim](figures/F2_raw_behind_claim.png)

*Every observation with its box summary; right panel shows mean ± SD against
mean ± SEM on the same axis.*

**Population.** *[What population the claim is about; whether the sample
plausibly represents it; what the methods do and do not state.]*

**Unit of replication.** *[Reported n and unit; rows vs units in the data; ICC,
design effect, effective n; whether the reported df are consistent with the
stated n.]*

**Uncertainty reporting.** *[What the paper says its error bars are; what they
match; whether that matches the claim being made.]*

---

## 5. Method and assumptions — stages 3–4

![F3 — distribution shape](figures/F3_assumption_shape.png)

*Q-Q per group (or on the paired differences), with n, skewness, excess kurtosis
and the Shapiro-Wilk p-value. Descriptive, not a pass/fail gate.*

**Test selection.** *[Named test; whether it matches one-sample/two-sample/
paired, tail, and variance assumption; whether its null matches the claim's
wording.]*

**Three routes.**

| route | test | p | assumes | reachable at 0.05 |
|-------|------|---|---------|-------------------|
| parametric | | | | |
| non-parametric (ranks) | | | | |
| permutation | | | | |

*[Whether they agree; if not, which assumption the claim rests on.]*

---

## 6. Resolution and evidence — stage 5

![F5 — null distribution](figures/F5_null_distribution.png)

*Permutation null with the observed statistic marked, annotated with the
p-value, the resolution limit and the design floor where it binds; right,
bootstrap distribution of the estimate with its interval.*

**Resolution.** *[Design floor; whether it binds; whether the reported p is
achievable at the stated permutation count; what shuffling scheme was used and
what null that encodes.]*

**Interval.** *[Bootstrap point estimate, CI, method, n_boot, bias.]*

---

## 7. Multiplicity — stage 6

![F6 — multiplicity](figures/F6_multiplicity.png)

*p-value histogram against the uniform expectation with the just-below-threshold
window highlighted; right, hits surviving each correction against the number
expected by chance.*

| correction | controls | n significant | threshold on raw p |
|------------|----------|---------------|--------------------|
| none | per-test error rate only | | |
| Benjamini-Hochberg | FDR | | |
| Bonferroni | FWER | | |

**Counting.** *[m tests run, including panels and supplementary; 0.05 × m;
which claims fail which correction.]*

**Distribution of reported p-values.** *[Counts in the two windows, binomial p,
and the three interpretation limits. No attribution of intent.]*

---

## 8. Per-claim supplement

![C1](figures/panel_C1.png)

*Raw observations, permutation null, and the three routes for C1. Shaded region
marks p-values unreachable for this design.*

*[One per headline claim, and per claim whose verdict is not `supported`.]*

---

## 9. What was not assessed

*[Paste `assumption_checks.uncheckable_assumptions_note()` verbatim.]*

Additionally:

- Claims with no supporting data: *[ids]*
- Analyses the supplied data could not reconstruct: *[what]*
- Anything subsampled, filtered or approximated during this review: *[what, and
  why]*

---

## 10. Reproduction

```bash
python workflows/run_review.py --pdf data/raw/manuscript.pdf \
                               --data data/raw/source.csv \
                               --out results/ --figures figures/
```

| | |
|--|--|
| Seeds | `random_state=[s]` |
| Permutations / resamples | `n_perm=[N]`, `n_boot=[N]` |
| Package versions | numpy [v], scipy [v], pandas [v], statsmodels [v] |
| Data checksum | sha256 `[…]` |
| Skill version | statistical-review v1 |

---

*Findings in this report are arithmetic statements about reported and recomputed
values. Where they differ, the cause is often an input the review could not see;
the inputs used are named above so that discrepancies can be resolved.*
