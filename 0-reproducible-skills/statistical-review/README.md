# Statistical Review

Evaluate the statistical claims of a paper against its own raw data.

**Inputs:** a paper PDF and the raw data behind its figures.
**Output:** a figure-led review report with a verdict for every claim.

## Organizing principle: three questions per claim

A statistical claim can fail in three independent ways, and conflating them
produces reviews that argue about decimal places while missing the design.

| Question | Needs | Stage |
|----------|-------|-------|
| Is the number right? | PDF, then data | 0, 2 |
| Is the method right for the design? | methods section + data | 1, 3, 4 |
| Could the number have been anything else? | data + the full set of tests run | 5, 6 |

Most findings come from the second and third. A paper whose arithmetic is
flawless can still report a p-value its design could never have produced, or run
forty tests and report the four that crossed 0.05.

## Directory structure

```
statistical-review/
├── README.md                              ← this file
├── SKILL.md                               ← loadable skill definition
├── REVIEW_PROTOCOL.md                     ← the six-stage procedure
├── REPORT_FORMAT.md                       ← report structure and figure rules
├── kernel.py                              ← the workflow, as a few calls
├── pytest.ini
├── checks/
│   ├── 01_population_and_sampling.md      ← what population? representative?
│   ├── 02_replication_and_uncertainty.md  ← what is n? SD or SEM?
│   ├── 03_test_selection.md               ← does the test match the design?
│   ├── 04_assumptions_and_three_routes.md ← parametric / rank / permutation
│   ├── 05_resampling_and_resolution.md    ← floors, resolution, bootstrap
│   └── 06_multiplicity_and_phacking.md    ← corrections, p-value distribution
├── templates/
│   ├── REVIEW_REPORT_TEMPLATE.md          ← report skeleton to fill in
│   └── claim_ledger_schema.md             ← ledger columns and verdict rules
├── utilities/
│   ├── claim_extraction.py                ← parse claims out of the PDF
│   ├── recompute.py                       ← statcheck, GRIM, error-bar identity
│   ├── assumption_checks.py               ← shape, spread, clustering
│   ├── resampling.py                      ← permutation, sign, bootstrap, floors
│   ├── multiplicity.py                    ← BH/Bonferroni, max-statistic, p-curve
│   └── review_figures.py                  ← the six figures + per-claim panel
└── tests/
    └── test_statistical_review.py         ← 38 tests
```

## Workflow

### 1. Establish what you were given

```python
from kernel import load_paper, load_data, describe_data

paper = load_paper("data/raw/manuscript.pdf")
paper["text_quality"]      # {'usable_text_layer': True, 'chars_per_page': 3184, ...}
paper["summary"]           # rows found per claim type
```

A scanned PDF yields an empty claim table. That is a reportable condition, not a
finding of "no claims" — `pdf_text_quality` exists so the two cannot be confused.

### 2. Audit the paper's own arithmetic

```python
from kernel import audit_reported_statistics

audit = audit_reported_statistics(paper["inventory"])
audit["n_checked"], audit["n_inconsistent"], audit["n_decision_errors"]
```

`statcheck_table` recomputes each p-value from its reported statistic and df.
Inconsistencies that cross 0.05 are flagged separately as `decision_error` —
those are the ones that can change a conclusion. Rows that cannot be checked
come back as NA and are counted, never dropped.

### 3. Look at the data

```python
desc = describe_data(data, value_col="volume", group_col="arm", unit_col="mouse")
desc["independence"]["effective_n"]     # 6.3, where the paper's test used n=300
```

### 4. Review each claim

```python
from kernel import review_claim

c1 = review_claim(data, value_col="volume", group_col="arm", unit_col="mouse",
                  claim_id="C1", claim_short="treated tumours are smaller",
                  p_reported=0.033, figure_path="figures/panel_C1.png")
c1["verdict"]    # 'assumption-dependent'
c1["reasons"]    # ['the routes disagree at alpha: Welch's t-test p = 0.031; ...']
```

`review_claim` runs the three routes, the permutation null with its resolution
floor, a bootstrap interval, and the assumption diagnostics, then proposes a
verdict from the evidence. The proposal is bookkeeping over computed
diagnostics — it cannot read the methods section, so the reviewer confirms it.

### 5. Build the report

```python
from kernel import build_review_figures, review_summary_table, verdict_counts

figs = build_review_figures(data, "mouse", "arm", "volume",
                            claims=[c1, c2], audit=audit, outdir="figures",
                            value_label="tumour volume (mm3)")
review_summary_table([c1, c2]).to_csv("results/claim_ledger.csv", index=False)
```

Then fill `templates/REVIEW_REPORT_TEMPLATE.md`.

## Design decisions

**Normality checks return descriptions, not verdicts.** There is no function
that returns "assumptions met". `distribution_shape` returns skewness, excess
kurtosis, a Shapiro-Wilk p-value where it is defined, and a note stating what
that p-value can mean at this n — at n < 30 a non-significant result is not
evidence of normality, and at large n a significant one need not invalidate a
test about means.

**Resolution floors are first-class outputs.** Every p-value this module
produces arrives with the smallest p its design could have produced.
`min_achievable_p(3, 3)` is 0.10; `min_sign_test_p(5)` is 0.0625. A claim whose
design cannot reach alpha gets its own verdict, `underpowered-resolution`,
rather than being argued about as a borderline result.

**The three routes are reported together.** `three_routes` returns one row per
route with its p-value and what it assumes. The module deliberately offers no
"choose the best test" function: selecting a route after seeing the p-values is
the behaviour the review exists to detect.

**How you permute is recorded.** `permutation_test(..., strata=...)` makes the
null explicit. Shuffling across the whole dataset tests a different hypothesis
from shuffling within patient, and the result dict carries `stratified` so the
report can state which null was tested.

**Multiplicity is counted, not asserted.** `discovery_recount` prints the hit
count under no correction, BH and Bonferroni side by side, with 0.05 × m as the
chance expectation. `just_below_rate` measures pile-up below the threshold and
refuses to interpret itself below 20 p-values.

**The figure spine is fixed.** All six figures are attempted on every review.
A slot that cannot be built renders as a panel naming the reason, because in a
report that simply omits a figure, an impossible diagnostic and an unconsidered
one look identical.

## What it cannot do

Independence, representativeness and confounding cannot be checked from the data
— only from the methods section. `uncheckable_assumptions_note()` returns the
disclaimer, and `REPORT_FORMAT.md` requires it verbatim in section 9 of every
report. A review that does not say what it could not check overstates itself.

## Testing

```bash
python -m pytest -c pytest.ini --rootdir . tests/
# 38 passed
```

The suite pins the lecture's exact numbers: the 3 vs 3 permutation floor of
0.10, the five-pair sign-test floor of 0.0625, the dataset with identical means
that differs by a count statistic, and the paired design where 9 of 10 subjects
improve while the mean change is negative.

---

**See also:** `../data-visualization/` for the figure principles, and
`../MODULES.md` for the rest of the toolkit.
