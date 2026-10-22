# Report Format

The review report is **figure-led**. A statistical finding is nearly always a
statement about a distribution, a sample, or a set of tests, and all three are
read faster from a figure with the numbers written on it than from a paragraph.
Prose explains what each figure shows and what it does not; it does not restate
the numbers already on the axes.

Figures follow the `data-visualization` skill. Load its `universal_principles.md`
before building any figure here, and the question-specific `PRINCIPLE.md` for the
module each figure draws on.

## Structure

```
1. Header                  inputs, what was and was not assessable
2. Summary of findings     verdict counts + the three findings that matter
3. Claim ledger            F4 + the ledger table: one row per claim
4. Design and sampling     F1 + stage 1-2 findings
5. Method and assumptions  F3 + stage 3-4 findings
6. Resolution and evidence F5 + stage 5 findings
7. Multiplicity            F6 + stage 6 findings
8. Per-claim supplement    one panel per headline claim
9. What was not assessed   uncheckable_assumptions_note(), verbatim
10. Reproduction           code, seeds, versions, data checksums
```

## 1. Header

State the inputs and their condition before anything else:

- PDF path, page count, whether the text layer was usable.
- Data path(s), row counts, checksums, and **which figure or table of the paper
  each file corresponds to**.
- Claims extracted automatically vs added by hand.
- Claims with no supporting data, listed by id. These are `not-assessable`, and
  saying so up front stops the report from reading as if it covered the paper.

## 2. Summary of findings

A verdict-count table, then **at most three** findings — the ones that would
change a reader's interpretation of the paper. Each is one sentence naming the
claim, the discrepancy as two numbers, and the stage it came from.

Everything else goes in the ledger. A summary that lists fifteen findings has not
decided which ones matter.

## 3. Claim ledger — F4

The spine of the report. One row per claim:

| column | content |
|--------|---------|
| `claim_id` | C1, C2, … — used everywhere else in the report |
| `claim_short` | the claim in the paper's own terms, under 60 characters |
| `verdict` | one of the seven; never blank, never two |
| `effect_reported` / `effect_recomputed` | in the paper's units |
| `ci_low` / `ci_high` | recomputed interval |
| `p_reported` / `p_recomputed` | both, always, even when they agree |
| `reasons` | why this verdict and not the adjacent one |

F4 renders the same content: a forest plot of recomputed effects with the
reported point overlaid, and reported against recomputed p-values on log axes
with the 0.05 lines drawn.

## 4–7. The body

One section per stage group, each led by its figure. The prose answers three
questions about the figure and stops:

1. What does this figure show?
2. What does it imply for the claims in the ledger, by `claim_id`?
3. What would resolve the discrepancy?

## 8. Per-claim supplement

One `fig_claim_panel` per headline claim: raw observations, the permutation
null, and the three routes side by side. For a paper with many claims, include
panels for the claims in the summary and the ones whose verdict is not
`supported`; list the rest in the ledger only.

## 9. What was not assessed

`assumption_checks.uncheckable_assumptions_note()`, verbatim, plus:

- claims with no data
- analyses the data could not reconstruct (e.g. the paper's own model fit)
- anything subsampled, filtered or approximated during the review, and why

## 10. Reproduction

Seeds, `n_perm`, `n_boot`, package versions, data checksums, and the script that
regenerates every figure. A review that cannot be re-run is an opinion.

---

## Figure rules specific to this report

On top of the universal principles:

- **Both numbers, always.** Any figure comparing reported to recomputed shows
  both, never the difference alone.
- **The floor is drawn, not described.** Where a design's p-value floor is above
  alpha, the figure marks the unreachable region. A reader should see that the
  p-value could not have been smaller.
- **Verdict colour is never the only channel.** Verdict is written as text next
  to every marker. The palette is Okabe-Ito and survives greyscale printing.
- **Error bars are named on the panel.** A figure showing error bars states what
  they are. Where the paper did not say, the panel says that instead.
- **n on every panel.** Group n, and unit n where they differ.
- **A slot that cannot be filled is still rendered**, as a panel naming the
  reason. The spine is fixed so a missing diagnostic cannot be mistaken for one
  that was considered and passed.

## Tone

Findings are arithmetic, not character assessment. Write what was reported, what
was recomputed, and what inputs were used; let the authors explain the gap.
Never infer intent from a p-value distribution — p-hacking diagnostics describe
a *set* of results and cannot speak to any one of them, nor to why it looks the
way it does.
