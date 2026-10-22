# Review Protocol

Six stages, in this order. The order matters: a design problem found in stage 1
or 2 changes how every later number should be read, and running the arithmetic
first invites you to argue about decimal places while the real finding is that
n = 3 means something different from what the paper implies.

**Inputs:** the paper PDF and the raw data. Nothing else is assumed.

**Before stage 1:** look at the data. Load it, plot it, count the rows, find the
missing values. The lecture's first instruction — always look at the data before
running numbers — applies to the reviewer as much as to the author.

```python
data = load_data("data/raw/source.csv")
print(describe_data(data, value_col=..., group_col=..., unit_col=...))
```

---

## Stage 0 — Establish what you were given

```python
paper = load_paper("data/raw/manuscript.pdf")
print(paper["text_quality"])
print(paper["summary"])
```

Record, in the report header:

- Whether the PDF has a usable text layer. If not, say so: claims could not be
  extracted automatically, and the claim table was built by hand (or was not
  built at all). An empty claim table from a scanned PDF is **not** "no claims
  found".
- What the raw data actually contains, and which figures or tables of the paper
  it corresponds to. Data that does not correspond to a specific claim cannot
  review that claim.
- Which claims have no supporting data at all. These get `not-assessable`, and
  the report names what would be needed. This is a finding in its own right.

---

## Stage 1 — Population and sampling

→ `checks/01_population_and_sampling.md`

Before any test: what population is the claim about? The paper must define it,
and the sample must plausibly be representative of it.

Nothing in `utilities/` can check this. It is read from the methods section, and
if the methods do not say, the report says they do not say.

**Finding template:** "The paper claims X about population P. The sample was
drawn as S, which [does / does not] plausibly represent P because …"

---

## Stage 2 — Replication and uncertainty

→ `checks/02_replication_and_uncertainty.md`

```python
desc = describe_data(data, value_col, group_col, unit_col="mouse")
desc["independence"]          # ICC, design effect, effective n
```

Two questions:

- **What is n?** The n a test may use is the number of independently randomized
  units, not the number of measurements. `design_independence` returns the ICC
  and the implied effective n; when rows are cells within animals and the ICC is
  non-trivial, the effective n is far below the row count.
- **What do the error bars describe?** SD describes the spread of observations;
  SEM describes uncertainty in the mean and shrinks as n grows. They differ by
  sqrt(n). `error_bar_identity(reported_error, sd, n)` says which the paper's
  number actually is.

---

## Stage 3 — Test selection

→ `checks/03_test_selection.md`

```python
paper["inventory"]["named_tests"]     # which tests the paper says it used
```

Check that the named test matches the design on every axis the lecture lists:
one-sample / two-sample / paired; one-tailed or two-tailed; equal or unequal
variance. Then check that its null hypothesis is the one the scientific question
needs — Mann-Whitney is not a test of equal medians, and a test of means does
not answer "do most subjects improve".

---

## Stage 4 — Assumptions and the three routes

→ `checks/04_assumptions_and_three_routes.md`

```python
routes = three_routes(a, b)                 # or paired=True
shape  = shape_table({"treated": a, "control": b})
var    = variance_comparison(a, b)
```

Run the parametric, rank-based and permutation routes on the same data and
report all of them. Agreement makes the claim robust to the modelling choice.
Disagreement localises exactly which assumption the claim depends on — that is
the finding, and it is reported as `assumption-dependent`, not as "the paper used
the wrong test".

Q-Q plots describe shape. They are not a gate. Do not write "normality was
violated, therefore the test is invalid"; write what the shape is, what it
implies at this n, and whether the routes disagree.

---

## Stage 5 — Resampling and resolution

→ `checks/05_resampling_and_resolution.md`

```python
perm = permutation_test(a, b, n_perm=10000)
perm["min_achievable_p"], perm["floor_binds"]
boot = bootstrap_ci(a, np.mean, n_boot=10000)
```

Before interpreting any small p-value, ask whether it could have been small:

- The permutation floor for the design — `min_achievable_p(n1, n2)`. For 3 vs 3
  it is 0.10, so a reported p = 0.0004 from such a design comes entirely from the
  parametric model, not from anything the data can demonstrate on their own.
- The sign-test floor — `min_sign_test_p(n_pairs)`. Five pairs cannot go below
  0.0625.
- The resolution of the reported permutation analysis — p = (k+1)/(N+1), so
  resolving p < 1e-6 needs about 1e6 permutations. A paper reporting p < 1e-8
  from 10,000 permutations is reporting something its procedure cannot produce.

Also check **how** the paper permuted. Shuffling across the whole dataset
encodes a different null from shuffling within patient or within cell type. If
the paper does not say how it permuted, the null it tested is unknown.

---

## Stage 6 — Multiplicity and analytic flexibility

→ `checks/06_multiplicity_and_phacking.md`

```python
audit = audit_reported_statistics(paper["inventory"])
audit["multiplicity"]["summary"]
audit["just_below"]
```

Count every test the paper ran, including those only visible as panels or
supplementary tables. Put 0.05 × m next to the paper's hit count, then recount
under BH and Bonferroni and report what changes.

Where the paper ran many correlated tests (a genome-wide scan, a panel of
markers), `max_statistic_null` is the correction that respects the dependence,
and is what the report should ask for if the paper used Bonferroni on
neighbouring correlated features.

The just-below-threshold diagnostic is a property of a *collection* of p-values.
Report it as a count with its binomial p-value and the explicit caveat that it
says nothing about any individual result. Never attribute intent.

---

## Stage 7 — Assemble

```python
figs = build_review_figures(data, unit_col, group_col, value_col,
                            claims=[c1, c2, ...], audit=audit, outdir="figures")
ledger = review_summary_table(claims)
counts = verdict_counts(claims)
```

Fill `templates/REVIEW_REPORT_TEMPLATE.md`. The ledger is the spine of the
report; the figures are the evidence; the prose explains what each figure shows
and what it does not.

**Every report ends with `uncheckable_assumptions_note()` verbatim.** A review
that does not state what it could not check overstates its own coverage.

---

## How to write a finding

Four parts, in this order:

1. **What the paper claims**, quoted or closely paraphrased, with the page.
2. **What was computed**, with the inputs named: the test, the n, the statistic,
   the seed.
3. **The discrepancy**, as both numbers, never as an adjective.
4. **What would resolve it** — a different test, a stated error definition, the
   per-unit data, a correction across the panel.

Write "the reported p = 0.033 does not follow from t(18) = 2.31; recomputing
gives p = 0.0329, which agrees, but the same sentence reports n = 3 independent
experiments, so the 18 degrees of freedom imply the test was run on 20
measurements rather than on 3 units" — not "the statistics are sloppy".

A discrepancy is an arithmetic fact; its cause is usually something the authors
can explain in one sentence. Leave them that room.
