# Claim Ledger Schema

One row per claim. Produced by `kernel.review_summary_table(claims)`, rendered
as figure F4, and the spine of the report.

## Columns

| column | type | content |
|--------|------|---------|
| `claim_id` | str | `C1`, `C2`, … Stable throughout the report; every figure and finding refers to claims by this id. |
| `claim_short` | str | The claim in the paper's own terms, under 60 characters. Not your paraphrase of what it should have said. |
| `verdict` | str | Exactly one of the seven below. Never blank, never two. |
| `effect_reported` | float | In the paper's units. NaN where the paper reports no effect size — which is itself worth noting. |
| `effect_recomputed` | float | The same quantity from the raw data. |
| `ci_low`, `ci_high` | float | Bootstrap interval on the recomputed effect. |
| `p_reported` | float | As printed. Where the paper gives a bound (`p < 0.001`), record the bound and note the relation. |
| `p_recomputed` | float | The permutation p-value — the route that assumes least. |
| `reasons` | str | Why this verdict and not the adjacent one. Pipe-separated. |

Optional columns worth adding per project: `page`, `figure_panel`,
`n_reported`, `n_effective`, `survives_correction`.

## Verdicts

Assigned by `kernel.propose_verdict`, first match wins. The reviewer confirms or
overrides after reading the methods — the function cannot see them.

| verdict | assigned when | colour |
|---------|---------------|--------|
| `supported` | recomputation reproduces the result and every route agrees at alpha | green `#009E73` |
| `supported-with-caveat` | the conclusion stands, but clustering, an error-bar mismatch, or a missing correction qualifies it | light blue `#56B4E9` |
| `assumption-dependent` | the parametric, rank and permutation routes disagree at alpha | orange `#E69F00` |
| `underpowered-resolution` | the design's permutation floor is above alpha, so no result from it could have been significant | pink `#CC79A7` |
| `not-reproducible` | the recomputed p falls on the other side of alpha from the reported one | vermillion `#D55E00` |
| `contradicted` | the recomputed effect has the opposite sign to the reported one | black `#000000` |
| `not-assessable` | the supplied data cannot address the claim | grey `#999999` |

Palette is Okabe-Ito; verdict is always also written as text, so colour is never
the only channel.

## Rules

1. **One verdict per claim.** A claim that is both underpowered and
   assumption-dependent takes the earlier verdict in the cascade; the other goes
   in `reasons`.
2. **`not-assessable` is a real verdict, not a gap.** Use it, and say in
   `reasons` what data would change it. Omitting the row instead makes the review
   look like it covered more than it did.
3. **Both p-values always appear,** even when they agree. A ledger showing only
   the recomputed value cannot be checked by the reader.
4. **`reasons` carries numbers, not adjectives.** "routes disagree at alpha:
   Welch p = 0.031; Mann-Whitney p = 0.08; permutation p = 0.07", not "weak
   evidence".
5. **A claim with no data gets a row, not silence.**

## Example

```python
import pandas as pd
from kernel import review_claim, review_summary_table

c1 = review_claim(data, "volume", "arm", unit_col="mouse", claim_id="C1",
                  claim_short="treated tumours are smaller", p_reported=0.033)
c2 = review_claim(data2, "score", "arm", claim_id="C2",
                  claim_short="marker predicts response", p_reported=0.029)

ledger = review_summary_table([c1, c2])
ledger.to_csv("results/claim_ledger.csv", index=False)
```

Claims the data cannot address are appended by hand:

```python
ledger = pd.concat([ledger, pd.DataFrame([{
    "claim_id": "C3",
    "claim_short": "effect persists at 12 months",
    "verdict": "not-assessable",
    "p_reported": 0.041,
    "reasons": "no longitudinal data supplied; source data cover the 4-week endpoint only",
}])], ignore_index=True)
```
