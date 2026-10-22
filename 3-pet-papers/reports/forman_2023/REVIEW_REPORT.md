# Statistical review: Forman, Renner & Leavens (2023), *Fetching felines*

Forman J., Renner E. & Leavens D.A. (2023) Fetching felines: a survey of cat owners
on the diversity of cat (*Felis catus*) fetching behaviour. *Scientific Reports*
**13**:20456. doi:10.1038/s41598-023-47409-w

Reviewed with the `statistical-review` skill (six-stage protocol); figures follow the
`data-visualization` skill. Every number below was recomputed from the supplied
Supplementary S3 dataset; nothing is simulated, reconstructed or approximated.

---

## 1. Header: what was supplied, and what was assessable

| input | path | size | sha256 (first 16) | corresponds to |
|---|---|---|---|---|
| paper text | `v9629532c_forman_2023_paper_fulltext.txt` | 58,930 characters | `e08b1b2a7444df9d` | whole article (text layer) |
| raw data | `v2210645b_forman_2023_dataset_S3.csv` | 1154 rows x 75 columns | `c9451c33a9675e7b` | Supplementary S3; Table 1-2, Figures 2-3 and every test in Results |
| authors' analysis code | `vf76dcde0_forman_2023_supplementary_S2_Rcode.txt` | 122,185 characters | `8ba14b0bc577eac8` | Supplementary S2, extracted from the .docx: R source with console output interleaved |

The paper was supplied as an extracted text layer, not a PDF, so `load_paper` was bypassed and the claim extractors were run on the text directly. Page numbers are therefore unavailable and claims are located by section.

**Supplementary S2.** S2 was supplied after the first pass of this review and this
revision incorporates it. Word collapsed many in-paragraph line breaks, so it is not a
runnable script, but the printed console output survives intact: 23 test statistics
with their degrees of freedom and p-values, 8 pairwise-Wilcoxon comparison
matrices, and the call that produced each one. Every recomputed test in this report is now
checked against the authors' own call and printed result (section 5), and the one
discrepancy the first pass could not resolve -- the Mann-Whitney Z -- is settled.

**Claim extraction.** The automatic extractors found
15 test statistics,
30 p-values and
5 named tests.
Four of the extracted p-values are artefacts of the text layer (three numbers inside the
Table 1 demographic block, and the alpha threshold quoted in the Levene sentence); they
are excluded, leaving **28 reported p-values** from
**17 named tests plus 11 pairwise Wilcoxon comparisons**. The
15 reported test
statistics were each recomputed by hand against the raw data.

**Coverage.** 25 claims are in the ledger. They cover all
15 reported test
statistics, all 28 reported p-values (the 11 pairwise Wilcoxons are
recomputed inside C12-C17), the four headline descriptives in the abstract, and ten
Table 1 descriptives cross-checked in section 4. 3 claims are
`not-assessable` and are listed by id in section 9.

**Not covered.** The object-category tallies and thematic map (the paper's Figures 1-4
and the object groupings) are descriptive or qualitative; no statistical claim rests on
them and none was re-derived beyond the breed and training counts in section 4. S2 cannot
be executed as supplied (collapsed line breaks), so each test is still re-implemented
rather than re-run -- but it is now re-implemented against a known call, and the result is
compared with the authors' printed output test by test.

---

## 2. Summary of findings

| verdict | n_claims |
|---|---|
| supported | 17 |
| supported-with-caveat | 5 |
| not-assessable | 3 |

Three findings would change how a reader reads the paper.

1. **The unit of analysis is the cat; the sampling unit is the owner; and S3 cannot
   connect them (C23, stage 1-2).** 924 owners reported on 1154 cats
   (1.249 cats per owner), 269 rows are
   flagged as repeating a previous respondent's demographics, and 186 cats
   are recorded as sharing a household with at least one other fetching cat -- but S3
   contains no respondent, household or submission identifier, so no cluster-robust or
   mixed model can be fitted. What *can* be computed is the bound: the worst-case design
   effect (intra-cluster correlation = 1) is 1.249 from the owner/cat
   ratio and 1.188 from the recorded fetching-cat household sizes.
   Dividing every reported chi-square by 1.249 leaves all of them on the same
   side of alpha -- the sex-ratio test moves from p = 0.0185 to
   p = 0.0351. Clustering is therefore a reporting gap here, not a
   result-changing one.

2. **The outlier sensitivity analysis changes one interpretation, and the Methods say it
   changes none (C20, stage 5-6).** The Methods state that analyses were run with and
   without extreme outliers and that "there were no differences in the interpretation of
   the findings". For the sex difference in retrievals that is not the case. With outliers
   retained -- the reported analysis -- S2 prints W = 156213,
   p = 0.09177, which recomputes exactly
   (W = 156213, p = 0.09177).
   With extreme outliers excluded by the authors' own rule S2 prints
   W = 149449, p = 0.02865, which also recomputes
   exactly (W = 149449,
   p = 0.02865; 14 extreme values removed,
   leaving 526 female and 614
   male cats). That crosses alpha: the comparison the paper reports as non-significant is
   significant in its own sensitivity run, and S2's code comment on that line reads
   "significant differences between groups". The same check on monthly sessions does not
   change (p = 0.6226 retained,
   p = 0.5858 excluded).

3. **"Significantly fewer retrievals if the owner usually initiates" is quoted with two
   identical medians (C13, stage 3-4).** The sentence reads "(median = 5, IQR = 4)
   compared to the cat (median = 5, IQR = 6, Pairwise Wilcoxon p = 0.009)". The adjusted
   p reproduces (0.0093) and all five routes agree at alpha, so
   the result stands -- but it is a statement about stochastic ordering and the mean, not
   about medians: the recomputed mean difference is
   -1.31 retrievals
   (95% CI [-2.02,
   -0.58]), rank-biserial
   -0.129 with the owner group first.

Everything else is in the ledger. The headline results of the paper -- that fetching is
overwhelmingly untrained, that it is first seen in young cats, and that cats initiate and
terminate fetching more often than owners do -- reproduce exactly from S3.

---

## 3. Claim ledger -- F4

![F4 claim ledger]({{artifact:art_0cef1b0b-5e22-4539-a312-dc04eae63950}})

Left: the recomputed effect with its 95% bootstrap interval, one row per claim. Every
open marker is absent because **the paper reports no effect size for any of its tests**;
the only reported effects are the percentages in C1, C2 and C25. Right: reported against
recomputed p-values on log axes. No point lies in an off-diagonal quadrant across 0.05,
so no recomputation changes a decision. The points sitting above the diagonal are the
claims whose reported p is a `< 0.001` bound recorded at 0.001 while the recomputed value
is smaller, and C20/C21, where the paper's Mann-Whitney p and the review's rank
permutation p differ in the third decimal.

| claim_id | claim_short | verdict | effect_metric | effect_reported | effect_recomputed | ci_low | ci_high | p_reported | p_reported_relation | p_recomputed | p_smallest_possible | survives_bonferroni_m28 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C1 | 94.4% of fetching emerged without training | supported-with-caveat | proportion untrained | 0.944 | 0.9437 | 0.9288 | 0.9563 |  |  |  |  |  |
| C2 | fetching first noticed <1 yr (701) or 1-7 yr (415) | supported | proportion <1 yr |  | 0.6075 | 0.5786 | 0.6358 |  |  |  |  |  |
| C3 | unequal distribution of who initiates fetching | supported | Cramer's V |  | 0.2321 | 0.1916 | 0.2735 | 0.001 | < | 1e-05 | 1e-05 | yes |
| C4 | cats initiate more often than owners | supported | Cramer's V |  | 0.3791 | 0.3167 | 0.4414 | 0.001 | < | 1e-05 | 1e-05 | yes |
| C5 | cats initiate more often than 'about equal' | supported | Cramer's V |  | 0.2221 | 0.158 | 0.2862 | 0.001 | < | 1e-05 | 1e-05 | yes |
| C6 | owners initiate less often than 'about equal' | supported | Cramer's V |  | 0.1714 | 0.09484 | 0.2479 | 0.001 | < | 2e-05 | 1e-05 | yes |
| C7 | unequal distribution of who ends fetching | supported | Cramer's V |  | 0.3776 | 0.3367 | 0.4216 | 0.001 | < | 1e-05 | 1e-05 | yes |
| C8 | cats terminate more often than owners | supported | Cramer's V |  | 0.4626 | 0.4063 | 0.519 | 0.001 | < | 1e-05 | 1e-05 | yes |
| C9 | cats terminate more often than 'about equal' | supported | Cramer's V |  | 0.4901 | 0.4327 | 0.5453 | 0.001 | < | 1e-05 | 1e-05 | yes |
| C10 | owner vs 'about equal' termination: no difference | supported | Cramer's V |  | 0.03549 | 0.002088 | 0.1232 | 0.437 | = | 0.464 | 1e-05 | no |
| C11 | more male (617) than female (537) fetching cats | supported-with-caveat | Cramer's V |  | 0.06932 | 0.01386 | 0.1265 | 0.018 | = | 0.01958 | 1e-05 | no |
| C12 | retrievals differ by who initiates | supported | epsilon-squared |  | 0.0104 | 0.001822 | 0.0273 | 0.001 | < | 0.0005999 | 9.999e-05 | yes |
| C14 | retrievals differ by who ends | supported | epsilon-squared |  | 0.09431 | 0.06493 | 0.1298 | 0.001 | < | 9.999e-05 | 9.999e-05 | yes |
| C16 | monthly sessions differ by who initiates | supported | epsilon-squared |  | 0.04426 | 0.02418 | 0.07191 | 0.001 | < | 9.999e-05 | 9.999e-05 | yes |
| C17 | monthly sessions differ by who ends | supported | epsilon-squared |  | 0.05297 | 0.03118 | 0.08188 | 0.001 | < | 9.999e-05 | 9.999e-05 | yes |
| C13 | fewer retrievals when the owner initiates | supported-with-caveat | rank-biserial (claim direction positive) |  | 0.1293 |  |  | 0.009 | = | 0.0031 | 9.999e-05 | no |
| C15 | more retrievals when the owner ends the session | supported | rank-biserial (claim direction positive) |  | 0.4317 |  |  | 0.001 | < | 1.06e-23 | 3.568e-232 | yes |
| C18 | no sex difference in who initiates | supported | Cramer's V |  | 0.03967 |  |  | 0.403 | = | 0.3978 | 9.999e-05 | no |
| C19 | no sex difference in who ends | supported | Cramer's V |  | 0.04128 |  |  | 0.374 | = | 0.3708 | 9.999e-05 | no |
| C20 | no sex difference in retrievals | supported-with-caveat | rank-biserial (male - female) |  | 0.05705 |  |  | 0.091 | = | 0.09639 | 9.999e-05 | no |
| C21 | no sex difference in monthly sessions | supported | rank-biserial (male - female) |  | 0.01674 |  |  | 0.622 | = | 0.6274 | 9.999e-05 | no |
| C22 | cats show independent, co-ordinated agency | not-assessable |  |  |  |  |  |  |  |  |  |  |
| C23 | 1154 cats analysed as independent observations | not-assessable |  |  |  |  |  |  |  |  |  |  |
| C24 | thematic analysis yields three main themes | not-assessable |  |  |  |  |  |  |  |  |  |  |
| C25 | 994 (86.1%) of fetching cats were mixed breed | supported-with-caveat | proportion mixed-or-unknown | 0.861 | 0.8614 | 0.8401 | 0.8808 |  |  |  |  |  |

**Sign convention.** `effect_recomputed` is oriented so that a positive value means "in
the direction the paper claims" wherever `effect_metric` says so (C13, C15); the raw
rank-biserial in the `reasons` column carries the sign of the group order the test was
run in, and both signs of the same statistic are given there so the two can be matched.
For C20 and C21 -- null claims, where there is no claimed direction -- the raw
male-minus-female sign is kept.

`reasons` for every row is in `claim_ledger.csv`; `p_smallest_possible` is the resolution
of the recomputed reference route (Monte-Carlo or permutation), i.e. the smallest p-value
that procedure could have returned.

---

## 4. Design and sampling -- F1

![F1 design map]({{artifact:art_80ca1406-8ba6-4188-ae48-11fb86352188}})

**Population.** The paper claims diversity of fetching behaviour among cats that fetch.
The sample was recruited through Twitter, Facebook and Reddit from owners who currently
or previously owned a fetching cat, so the population is self-selected owners of fetching
cats who use those platforms; 724 of 1154 cats
are in the USA. The paper states this limit explicitly and does not claim prevalence.

**Unit of replication.** The left panel is the sampling ladder: 1,258 owners consented,
924 completed (334 excluded as incomplete), reporting 1,184 cats, of
which 1154 survive the stated exclusions. Two things follow. First, the "924
owners / 1154 cats" pairing in the abstract combines a pre-exclusion owner count with
a post-exclusion cat count; the number of owners contributing to the final 1154 cats
is not reported and cannot be derived from S3. Second, cats are nested within owners and
no identifier links them (C23): the columns searched for an id are
`id/owner/respondent/participant/response: none found`.
The design-effect bound in the panel is the strongest statement the supplied data allow.

**Exclusions** are stated and counted in the paper (334 incomplete, 2 extreme, 23
ambiguous, 5 duplicate) and the arithmetic closes: 1,184 - 30 = 1154.

**Confounding.** Not checkable. Nothing in S3 identifies when a response was collected,
and the grouping variables of interest (who initiates, who ends) are themselves
owner-reported, so owner interpretation is shared between predictor and outcome.

**Descriptive cross-check.** Ten figures from the abstract and Table 1, recomputed:

| reported in paper | as printed | recomputed from S3 |
|---|---|---|
| untrained cats | 1089 (94.4%) | 1089 (94.4%) |
| trained cats | 65 (5.6%) | 65 (5.6%) |
| fetching first noticed: kitten | 701 | 701 |
| fetching first noticed: young adult | 415 | 415 |
| current cat / previous cat | 853 (73.9%) / 301 (26.1%) | 853 / 301 |
| female spayed / female intact | 501 (43.4%) / 36 (3.1%) | 501 / 36 |
| male neutered / male intact / male unknown | 601 (52.1%) / 15 (1.3%) / 1 (0.1%) | 601 / 15 / 1 |
| mixed breed | 994 (86.1%) | 493 mixed + 501 unknown = 994 (86.1%) |
| purebred | 160 (13.9%) | 160 (13.9%) |
| respondent gender: female / male / non-binary | 748 (64.8%) / 266 (23.1%) / 109 (9.4%) | 748 / 266 / 109 |

Nine of the ten reproduce exactly. The tenth is C25: "994 (86.1%) mixed breed" is
493 cats recorded as mixed plus 501 recorded
as unknown. The Discussion quotes the 994 as "mixed breeds" in one sentence and notes
unknown breed as a limitation in the next.

---

## 5. Method and assumptions -- F2, F3

![F2 raw data behind the claim]({{artifact:art_6b79a4ca-b2eb-4518-89a6-3d345fd7d27c}})

![F3 distribution shape]({{artifact:art_3436b3c3-09e1-4126-9003-50bde438ea94}})

**What the paper did.** Normality was screened with Shapiro-Wilk tests and Q-Q plots, and
because "the Shapiro-Wilk tests determined the distributions of the data were
significantly different from a normal distribution", non-parametric tests were used
throughout: chi-square for frequencies, Kruskal-Wallis with Bonferroni-adjusted pairwise
Wilcoxon tests for the two count outcomes, Mann-Whitney for sex contrasts. Extreme
outliers (beyond Q3 + 3xIQR) were identified and the analyses were run with and without
them, with the outlier-retaining version reported.

**What the screen can and cannot establish at this n.** Shapiro-Wilk is a test of a null
of exact normality whose power rises with n; at the group sizes here
(Cat: 675; About equal: 231; Myself: 248) it will
reject departures far too small to matter for any test about means, and a non-significant
result would not have been evidence of normality either. The diagnostic that *is*
informative is the shape itself, and F3 shows it: Cat: n = 675, median = 5, skew = 3.32, excess kurtosis = 19.6, Shapiro p = 3.6e-32; About equal: n = 231, median = 6, skew = 4.20, excess kurtosis = 30.5, Shapiro p = 4.9e-21; Myself: n = 248, median = 8, skew = 2.29, excess kurtosis = 8.5, Shapiro p = 5.2e-17. These are bounded,
heavily right-skewed integer counts with a floor at zero
(4 cats retrieved 0 times;
35 had fewer than one monthly session) and
a long upper tail (maximum 68 retrievals,
180 monthly sessions). Rank-based tests are a
defensible choice for such data -- but because the outcome is a skewed count and the
question is about which group tends to be higher, not because a normality test returned
p < 0.05. This review does not treat the Shapiro result as a pass/fail gate in either
direction, and no verdict in the ledger turns on it.

**Whether the gate mattered.** It did not. Running the parametric route anyway, the
one-way ANOVA agrees with Kruskal-Wallis at alpha on all four omnibus comparisons
(C12 p = 0.0048 vs 0.000923;
C14 3.98e-18 vs 9.89e-25;
C16 5.77e-08 vs 3.18e-12;
C17 7.21e-10 vs 2.12e-14), and the label
permutation of H agrees with both. The conclusions do not depend on the modelling choice.

**Null versus question.** Three mismatches are worth the authors' attention, none of them
decision-changing. (i) Kruskal-Wallis and Mann-Whitney test stochastic ordering, not
equality of medians; the Results sentences quote medians either side of a Wilcoxon
p-value, which reads as a median comparison (C13 is the clearest case: both medians 5).
(ii) The chi-square goodness-of-fit tests in C3-C11 test a uniform null over the response
categories; "cats initiate more than owners" is a claim about that sample's composition,
which is what the test delivers, but the sex-ratio claim (C11) is phrased as a claim about
cats in general. (iii) The paper reports Levene's tests showing unequal variances as the
reason for choosing the pairwise adjustment; Levene on heavily skewed counts is itself
shape-sensitive, and the choice of a Bonferroni-adjusted rank test does not depend on it.

**Cross-check against the authors' own calls (S2).** Every test the paper reports can now
be matched to the call that produced it and to its printed output:

| claim | call in S2 | authors' printed statistic | recomputed here | authors' p | recomputed p | agree |
|---|---|---|---|---|---|---|
| C3 | X-squared | 124.29 | 124.286 | < 2.2e-16 | 1.027e-27 | yes |
| C4 | X-squared | 115.23 | 115.232 | < 2.2e-16 | 7.001e-27 | yes |
| C5 | X-squared | 44.642 | 44.642 | = 2.366e-11 | 2.366e-11 | yes |
| C6 | X-squared | 17.652 | 17.6522 | = 2.652e-05 | 2.652e-05 | yes |
| C7 | X-squared | 329.08 | 329.076 | < 2.2e-16 | 3.483e-72 | yes |
| C8 | X-squared | 197.54 | 197.54 | < 2.2e-16 | 7.191e-45 | yes |
| C9 | X-squared | 217.59 | 217.589 | < 2.2e-16 | 3.035e-49 | yes |
| C10 | X-squared | 0.60334 | 0.60334 | = 0.4373 | 0.4373 | yes |
| C11 | X-squared | 5.5459 | 5.54593 | = 0.01852 | 0.01852 | yes |
| C12 | Kruskal-Wallis chi-squared | 13.976 | 13.9761 | = 0.0009229 | 0.0009229 | yes |
| C14 | Kruskal-Wallis chi-squared | 110.55 | 110.546 | < 2.2e-16 | 9.889e-25 | yes |
| C16 | Kruskal-Wallis chi-squared | 52.947 | 52.9472 | = 3.182e-12 | 3.182e-12 | yes |
| C17 | Kruskal-Wallis chi-squared | 62.972 | 62.972 | = 2.117e-14 | 2.117e-14 | yes |
| C13 | pairwise Wilcoxon (Bonferroni) | - | - | 0.00934 | 0.00934 | yes |
| C15 | pairwise Wilcoxon (Bonferroni) | - | - | < 2e-16 | 1.06e-23 | yes |
| C18 | X-squared | 1.8163 | 1.81631 | = 0.4033 | 0.4033 | yes |
| C19 | X-squared | 1.9667 | 1.96669 | = 0.3741 | 0.3741 | yes |
| C20 | W (Wilcoxon rank sum) | 156213 | 156213 | = 0.09177 | 0.09177 | yes |
| C21 | W (Wilcoxon rank sum) | 162892 | 162892 | = 0.6226 | 0.6226 | yes |

All 19 agree. Statcheck gains from S2 too: the manuscript prints degrees of
freedom for only 2 of its 15 test statistics, so
most triples were not checkable from the text, whereas all 19 printed (statistic, df, p) triples in S2 are internally consistent, and none is a decision error. The calls
are `chisq.test` on one-dimensional tables for C3-C11
("Chi-squared test for given probabilities", i.e. goodness-of-fit against equal
probabilities -- which is what this review implemented), `chisq.test` on a 2x3 table for
C18-C19, `kruskal.test` for C12-C17, `pairwise.wilcox.test(p.adjust.method =
"bonferroni")` for the pairwise comparisons, and
`wilcox.test(paired = FALSE, exact = FALSE, conf.int = TRUE)` for C20-C21. Three
consequences:

- **C11 is a transcription slip, not an analysis error.** S2 prints X-squared =
  5.5459, df = 1, p-value = 0.01852,
  identical to the recomputed 5.5459. The manuscript's "X2 = 5.50" is
  5.5459 re-rounded. Continuity correction never entered it: `chisq.test` applies Yates
  only to 2x2 contingency tables, so the Yates alternative floated in the first pass of
  this review is not a candidate and has been withdrawn.
- **The Mann-Whitney Z is resolved.** R's `wilcox.test` prints W, not Z; the reported
  Z = -1.69 and Z = -0.49 are the normal deviates of the printed p-values
  (qnorm(p/2) = -1.6861 and
  -0.4922). The tie-corrected |Z| computed in the first pass
  (1.686, 0.492) is the matching quantity; the
  uncorrected one is not. There is no discrepancy.
- **The statcheck flag on C19 was a rounding artefact.** S2 prints X-squared =
  1.9667, p-value = 0.3741; the manuscript's 1.97
  and 0.374 are both correct roundings of it.

**Error bars.** None are reported anywhere: the paper gives medians with IQRs and means
with SDs for demographics, and no interval for any estimate. F2's right panel shows
mean +/- SD against mean +/- SEM for the headline outcome so the difference between the
two is visible; neither is what the paper plots, which is why the panel says so.

---

## 6. Resolution and evidence -- F5

![F5 null distribution]({{artifact:art_0b051274-94b9-4c7a-8f1e-11875a23f6e1}})

No design in this paper is resolution-limited: with n in the hundreds, the combinatorial
floor of every test is far below any reportable p-value (for C15, 675 vs 248 cats, the
two-sided permutation floor is 3.57e-232). What *is*
limited is the review's own resampling: the Monte-Carlo chi-square nulls resolve to
1.0e-05 (100,000 replicates) and the permutation routes to
1.0e-04 (10,000 replicates), so every recomputed p-value at those values
is a resolution limit, not an estimate -- which is why the ledger carries
`p_smallest_possible` next to `p_recomputed`.

The left panel makes the point for the strongest claim in the paper (C7, who ends
sessions): the observed chi-square 329.1 lies far outside the entire
Monte-Carlo null (maximum 24.9 over 100,000 draws). The
right panel answers the separate question the paper does not: how precise the effect is.
Cramer's V = 0.378
[0.337, 0.422] -- a moderate
association, reported here because the paper reports none.

The paper reports no permutation or bootstrap analysis, so there is no shuffling scheme
or iteration count of its own to check.

---

## 7. Multiplicity -- F6

![F6 multiplicity]({{artifact:art_07dc4d3c-81be-48b5-8e7e-8bc352ed56a2}})

The paper runs 28 reported tests, where 1.4
significant results are expected by chance alone if every null were true.
21 are significant uncorrected; 21 survive Benjamini-Hochberg and 19 survive
Bonferroni at the full m = 28 (threshold 0.00179). The two that
do not survive Bonferroni are C11 (sex ratio, p = 0.018) and C13 (pairwise retrievals,
p = 0.009) -- both already flagged `supported-with-caveat` for other reasons.

**The paper did correct, within families.** Post-hoc chi-squares were interpreted at a
corrected alpha of 0.016 (0.05/3) and the pairwise Wilcoxon tests carry Bonferroni
adjustment, which is the right correction for those families. What is not stated is a
global count: the recount above is the view across the whole Results section, and it is
offered as context rather than as a demand, since the comparisons are pre-specified by
the survey's structure rather than selected after the fact.

**What S2 shows about the true denominator.** The script contains
11 `chisq.test`, 8 `kruskal.test`,
8 `pairwise.wilcox.test` (three comparisons each) and
4 `wilcox.test` calls -- 47 inferential p-values in
total, against the 28 the paper prints. The gap is almost entirely the
outlier-excluded duplicate of each count analysis (4 Kruskal-Wallis,
12 pairwise and 2
Wilcoxon tests), which is a robustness check rather than a second pass at the hypothesis,
and the paper summarises it in one Methods sentence -- a sentence that is inaccurate for
C20 (section 2, finding 2). On top of those, S2 runs 24
Shapiro-Wilk and 8 Levene screens, whose printed results confirm the
Methods description (for instance Levene F = 3.6907, p = 0.02525 on retrievals by
initiator).

**One pairwise comparison is run but not reported.** For monthly sessions by initiator the
paper reports two of the three pairwise Wilcoxons; the third, cat versus about-equal, is
printed in S2 as p = 0.66 and recomputes as
p = 0.659 (Bonferroni-adjusted).
Reporting all three keeps the correction's denominator visible.

**Distribution of reported p-values.** 19 of the 28 p-values are
printed as `p < 0.001` bounds and are recorded here at 0.001, so the left-hand spike in
F6 is partly a reporting convention rather than a measured distribution. The
just-below-threshold diagnostic returns
0 p-values in [0.045, 0.05) against
0 in [0.040, 0.045) -- both windows are empty, so the
diagnostic is uninformative here and no binomial comparison is defined. It is in any case
a property of a *collection* of p-values, says nothing about any individual result, and
cannot distinguish analytic flexibility from a pre-specified analysis; no inference about
intent is drawn from it.

---

## 8. Per-claim supplement

One panel per headline claim and per claim whose verdict is not `supported`.

**C1 -- 94.4% of fetching emerged without training** (`supported-with-caveat`)

![C1 panel]({{artifact:art_bdc97d9e-83ef-46de-8598-7d4380f8e3f9}})

1089 of 1154 cats are recorded `Untrained` (94.37%, exact 95%
CI [92.88, 95.63]%), reproducing the reported
94.4% exactly. The caveat is the denominator: the abstract says "94.4% of these *owners*
report", and the figure is a proportion of *cats*. With no respondent identifier the
owner-level proportion cannot be computed from S3. The centre panel's null (50:50
training) is shown because the panel template tests it, not because the paper claims it.

**C3 / C7 -- who initiates and who ends** (`supported`)

![C3 panel]({{artifact:art_e2671478-a26b-4d27-b23b-138307798c71}})

![C7 panel]({{artifact:art_86b20ddd-e6e6-45a5-9565-ff4ad2460b68}})

Both omnibus statistics reproduce exactly (initiation 124.29 against
the reported 124.29; termination 329.08 against 329.08), as do all
six post-hoc comparisons, and every route agrees. Cramer's V is
0.232 for initiation and 0.378 for
termination: the asymmetry in who ends a session is the stronger of the two, which the
paper's narrative also emphasises.

**C11 -- more male than female fetching cats** (`supported-with-caveat`)

![C11 panel]({{artifact:art_0bb264e5-0aaf-454e-87f6-531229a0a801}})

The routes agree with each other (asymptotic 0.0185, exact binomial
0.0200, Yates 0.0200, Monte-Carlo
0.0196) and the asymptotic value is exactly what S2 prints
(5.5459, p = 0.01852). The Yates row is shown for
completeness only: `chisq.test` does not apply it to a one-dimensional table, so it was
never part of the authors' result. What remains is the scope of the claim -- a 50:50 null
inside a sample of fetching cats describes composition, and with no non-fetching
denominator it cannot support a claim about sex and fetching -- plus the manuscript's
mis-transcription of 5.5459 as 5.50 and the fact that p = 0.018 does not survive
Bonferroni across the 28 reported tests.

**C13 -- fewer retrievals when the owner initiates** (`supported-with-caveat`)

![C13 panel]({{artifact:art_9ef3c08e-41a8-4d71-b9de-ee5246f5d206}})

See finding 3. The right-hand panel is the evidence that the result is robust to the
modelling choice: Welch, Student, Mann-Whitney, mean-difference permutation and rank-sum
permutation all land below alpha.

**C20 -- no sex difference in retrievals** (`supported-with-caveat`)

![C20 panel]({{artifact:art_31a9f047-2dff-4250-86d2-49d42b1f16c2}})

The reported Mann-Whitney reproduces exactly, including the W statistic S2 prints
(W = 156213, p = 0.09177; recomputed
W = 156213, p = 0.09177), and
the reported Z = -1.69 is the normal deviate of that p-value
(-1.6861). The verdict is `supported-with-caveat` for a
different reason -- the paper's own outlier-excluded run of this same comparison gives
p = 0.02865, on the other side of alpha (finding 2). This panel is
also worth reading because the routes, while unanimous at
alpha, differ by an order of magnitude in p: the rank routes sit near
0.09 while the mean-based routes sit near
0.67. They are answering different questions --
whether male cats tend to rank higher, versus whether their mean differs -- and the mean
difference is 0.141 retrievals
(95% CI [-0.519,
0.793]). A null result at this n is compatible both with
no association and with an association too small for 1154 cats to resolve.

**C22, C23, C24 -- not assessable**

![C22 panel]({{artifact:art_29f4f64b-2234-49ce-9361-26a27ef321ef}})

![C23 panel]({{artifact:art_31178810-acab-4c34-9a1b-af87e8083cc5}})

![C24 panel]({{artifact:art_57a6bd1f-eaa1-409f-845a-bd7a6f1ccd32}})

---

## 9. What was not assessed

Diagnostics in this report describe distribution shape and spread only. They cannot evaluate:
  - independence of observations
  - the sample is representative of the stated population
  - absence of unmeasured confounding
  - measurements were not selected after seeing the outcome
These are design questions. Assess them from the methods section, and say so explicitly when the methods do not contain enough information to assess them.

**Claims the supplied inputs cannot address** (verdict `not-assessable`, with what would
change each):

- **C22 - "cats demonstrate independent and co-ordinated agency".** S3 records, once per
  cat, who *usually* initiates and ends sessions on a three-point scale. That locates a
  frequency asymmetry; it carries no within-bout sequence, no latency, and no observation
  independent of the owner's interpretation. Bout-level observational or video-coded data,
  or an experimental manipulation of owner availability, would address it.
- **C23 - 1,154 cats analysed as independent observations.** No respondent, household or
  submission identifier exists in S3. A respondent id column, or owner-level counts, would
  allow a cluster-robust or mixed model. The worst-case bound computed in section 4 is the
  substitute, not a replacement.
- **C24 - the thematic analysis.** Qualitative coding is outside the scope of a
  statistical review, and the paper reports 268 initial codes with 35 child codes but no
  inter-coder agreement statistic. The code frame plus a second coder's assignments would
  make it checkable.

**Not checked for any claim:** independence of observations (see C23), representativeness
of the self-selected sample, and confounding between owner interpretation and the
reported behaviour. All three are design questions readable only from the methods, and
the methods here are explicit about the recruitment route.

**Review-side limitations.** Supplementary S2 is available to this revision but cannot be
executed as supplied -- the .docx extraction collapsed in-paragraph line breaks -- so each
test is re-implemented in Python against the authors' printed call rather than re-run in R.
The printed statistics, degrees of freedom and p-values are intact, and section 5 compares
every one of them with the recomputed value. The outlier-excluded sensitivity analysis is
now reproduced for the two cat-sex comparisons (the only ones where it could change an
interpretation, and it does for one); the outlier-excluded Kruskal-Wallis and pairwise runs
are read from S2's printed output and not independently recomputed, since all four remain
significant in both versions. Nothing was subsampled, filtered or approximated.

---

## 10. Reproduction

- Script: `review_forman_2023.py` (this review, end to end: recomputation, figures, ledger,
  report).
- Seeds: `SEED = 0` for every permutation, Monte-Carlo and bootstrap draw.
- Resampling: `N_PERM = 10,000` (resolution 1.0e-04),
  `N_MC = 100,000` (resolution 1.0e-05),
  `N_BOOT = 10,000`; alpha = 0.05.
- Inputs: `v9629532c_forman_2023_paper_fulltext.txt` sha256 `e08b1b2a7444df9db13e4800d49eac8297323a062740ab802ebf1728e7f9e04b`;
  `v2210645b_forman_2023_dataset_S3.csv` sha256 `c9451c33a9675e7b9ebbdc2484f6ef1c87d77dff01de8649c0d09b9c98ba9f8e`;
  `vf76dcde0_forman_2023_supplementary_S2_Rcode.txt` sha256 `8ba14b0bc577eac851f2eae46c6dca30ad449ad4c4003f6fefd9b774cf683384`.
- Outputs: `claim_ledger.csv` (25 claims, now with the authors' printed statistic per
  row), `statcheck.csv` (manuscript-text triples plus the 23 S2 console triples),
  `s2_crosscheck.csv` (the section 5 table), `reported_pvalues.csv`.
- Versions: python 3.11.16, numpy 2.4.6,
  pandas 2.3.3, scipy 1.17.1,
  matplotlib 3.11.2.
- One review-side fix: `resampling.min_achievable_p` evaluates 2/C(n1+n2, n1) directly
  and overflows a float at these sample sizes, so it is wrapped with a log-space version
  that underflows to 0. The skill itself is unmodified.
- Figure references in this report resolve to artifact markers when
  `figure_artifact_map.json` is present and to relative paths otherwise.
