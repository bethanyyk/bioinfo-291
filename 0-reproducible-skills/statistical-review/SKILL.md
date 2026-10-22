---
name: statistical-review
description: "Evaluate the statistical claims of a paper against its own raw data. Takes a paper PDF plus the raw data, extracts every reported test, p-value, interval and sample size, re-derives them, runs parametric / rank-based / permutation routes side by side, checks resolution floors and multiple-testing, and produces a figure-led review report with a per-claim verdict ledger. Use when reviewing a manuscript, auditing a published analysis, or checking your own paper before submission."
---

# Statistical Review

Review the statistical claims of a paper. **Inputs: the paper PDF and the raw
data.** Output: a figure-led review report with a verdict for every claim.

## What this skill decides

For each claim the paper makes, three separate questions, which are easy to
conflate and must be answered separately:

1. **Is the number right?** Does the reported p-value follow from the reported
   statistic? Does the test, re-run on the authors' own data, return what the
   paper says?
2. **Is the method right for the design?** Does the test's null hypothesis match
   the scientific question? Is the n the test used the n that was randomized?
   Does the conclusion survive the parametric, rank-based and permutation
   routes, or does it depend on which was chosen?
3. **Could the number have been anything else?** Is the design's resolution
   floor above 0.05? How many tests were run, and how many hits survive
   correction?

A claim passes only if all three hold. Most review findings come from the second
and third, not the first.

## Quickstart

```python
# From the skill root, or after copying statistical-review/ into your project
from kernel import (load_paper, load_data, describe_data,
                    audit_reported_statistics, review_claim,
                    build_review_figures, review_summary_table, verdict_counts)

paper = load_paper("data/raw/manuscript.pdf")
print(paper["summary"])                  # how many claims of each kind were found
print(paper["text_quality"])             # is the text layer usable at all?

audit = audit_reported_statistics(paper["inventory"])
print(audit["n_checked"], audit["n_inconsistent"], audit["n_decision_errors"])

data = load_data("data/raw/figure2_source.csv")
print(describe_data(data, value_col="volume", group_col="arm", unit_col="mouse"))

c1 = review_claim(data, value_col="volume", group_col="arm", unit_col="mouse",
                  claim_id="C1", claim_short="treated tumours are smaller",
                  p_reported=0.033, figure_path="figures/panel_C1.png")
print(c1["verdict"], c1["reasons"])

figs = build_review_figures(data, "mouse", "arm", "volume",
                            claims=[c1], audit=audit, outdir="figures",
                            value_label="tumour volume (mm3)")
```

## The review protocol

Six stages, in order. `REVIEW_PROTOCOL.md` has the full procedure; each stage
has a checklist in `checks/`.

| Stage | Question | Checklist |
|-------|----------|-----------|
| 1 | What population is the claim about, and is the sample representative of it? | `checks/01_population_and_sampling.md` |
| 2 | What is n, and does the reported uncertainty describe the right thing? | `checks/02_replication_and_uncertainty.md` |
| 3 | Does the test match the design and the question? | `checks/03_test_selection.md` |
| 4 | Do the parametric, rank and permutation routes agree? | `checks/04_assumptions_and_three_routes.md` |
| 5 | Could the p-value have been smaller, and is the interval trustworthy? | `checks/05_resampling_and_resolution.md` |
| 6 | How many tests were run, and what survives correction? | `checks/06_multiplicity_and_phacking.md` |

## Report format

The report is figure-led: findings are presented as figures with the numbers
written on them, and prose explains what each figure shows. `REPORT_FORMAT.md`
specifies the structure; `templates/REVIEW_REPORT_TEMPLATE.md` is the skeleton
to fill in.

Six standard figures, **always attempted** — a slot that cannot be built is
rendered as a panel naming the reason, so a missing diagnostic stays visible:

| Figure | Shows | data-visualization module |
|--------|-------|---------------------------|
| F1 design map | observations per unit; units per group; effective n | `abundance_and_proportion/` |
| F2 raw behind claim | every observation; SD against SEM on one axis | `distributions/` |
| F3 assumption shape | Q-Q per group, or on the paired differences | `distributions/` |
| F4 claim ledger | reported vs recomputed effect and p, per claim | `relationships/` |
| F5 null distribution | permutation null, observed value, bootstrap interval | `distributions/` |
| F6 multiplicity | p-value histogram; hits surviving each correction | `patterns_across_features/` |

Plus one per-claim panel (raw data, null, three routes) for each headline claim.

## Verdict vocabulary

Every claim gets exactly one verdict. `propose_verdict` assigns the first that
applies; the reviewer confirms or overrides it after reading the methods.

| Verdict | Meaning |
|---------|---------|
| `supported` | reproduces, and every route agrees at alpha |
| `supported-with-caveat` | conclusion stands, but the design or the reporting qualifies it |
| `assumption-dependent` | parametric, rank and permutation routes disagree at alpha |
| `underpowered-resolution` | the design's p-value floor is above alpha; no result from it could have been significant |
| `not-reproducible` | recomputed p falls on the other side of alpha from the reported one |
| `contradicted` | the data point the other way |
| `not-assessable` | the supplied data cannot address the claim; say what is missing |

## What this skill will not do

- **It will not stamp assumptions as met.** Normality diagnostics return shape
  descriptions, never a pass/fail. A Shapiro-Wilk p-value is reported with the
  note that its power depends on n.
- **It cannot check independence, representativeness or confounding.** Those are
  design questions, readable only from the methods section.
  `assumption_checks.uncheckable_assumptions_note()` returns the disclaimer to
  paste into every report.
- **It does not pick the smallest p-value.** The three routes are reported
  together. Choosing among them after seeing the answers is the thing the review
  exists to detect.

## Files

```
statistical-review/
├── SKILL.md                  ← this file
├── README.md                 ← design rationale and worked example
├── REVIEW_PROTOCOL.md        ← the six-stage procedure
├── REPORT_FORMAT.md          ← report structure and figure spine
├── kernel.py                 ← load_paper, review_claim, build_review_figures
├── checks/                   ← one checklist per stage
├── templates/                ← report skeleton and ledger schema
├── utilities/                ← the computational core (see below)
└── tests/                    ← 38 tests, including the lecture's exact numbers
```

| Module | Provides |
|--------|----------|
| `utilities/claim_extraction.py` | `read_pdf_text`, `build_claim_inventory`, `find_reported_tests` |
| `utilities/recompute.py` | `statcheck_table`, `grim_check`, `error_bar_identity`, `effect_from_summary` |
| `utilities/assumption_checks.py` | `distribution_shape`, `qq_points`, `variance_comparison`, `design_independence` |
| `utilities/resampling.py` | `permutation_test`, `sign_test`, `bootstrap_ci`, `three_routes`, `min_achievable_p` |
| `utilities/multiplicity.py` | `discovery_recount`, `just_below_rate`, `max_statistic_null` |
| `utilities/review_figures.py` | the six figures, the per-claim panel, `not_assessable_figure` |

## Testing

```bash
cd 0-reproducible-skills/statistical-review
python -m pytest -c pytest.ini --rootdir . tests/
```

The suite pins the lecture's specific numerical facts, so the module cannot
drift from them: 3 vs 3 permutation has a two-sided floor of exactly 0.10; a
sign test on 5 pairs cannot go below 0.0625; identical means can still differ by
a count-based statistic; and 9 of 10 subjects can improve while the mean change
is negative.

---

See `../data-visualization/` for the figure principles this skill builds on, and
`../README.md` for the full reproducible-skills toolkit.
