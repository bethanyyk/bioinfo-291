# Statistical review — Wang & Horowitz (2026), *Dogs (Canis familiaris) distinguish conspecific emotional chemosignals*

*Sci Rep* **16**:11176 · doi:10.1038/s41598-026-41426-1 · review slug `wang_horowitz_2026`
Method: the six-stage `statistical-review` protocol; figures on `data-visualization`.
**Revision 2** — the authors' raw data arrived after the first pass. Every claim that was
`not-assessable` in revision 1 has been re-derived against their own files and re-issued.
Section 2.1 states which revision-1 inferences the data confirmed and which it overturned.

## 1. Header — what was supplied

| input | bytes | sha256 (first 16) | corresponds to |
|---|---|---|---|
| `paper_fulltext.txt` | 87,062 | `eda57f3291c4e127` | the article text (PMC extraction; figure images still absent) |
| `Part_1.xlsx` | 37,214 | `3ef960171f504113` | Figs 3–5: investigation times, 3 sheets |
| `Part_2.xlsx` | 79,458 | `fb594cabc647751d` | Fig 6 and the sex analyses: 8 behavioural sheets + 10 per-sex sheets |
| `statistical_analyses.R` | 7,892 | `bf63d0a9271a09a8` | the analysis as run |
| `README.md` | 5,933 | `c39250e094ef2b68` | variable dictionary |

**What the files contain.** `Part_1.xlsx` holds three sheets of the same trials at two different
units: `all` (122 rows, one per odour-pairing trial, **no subject column**),
`emotionalpairings` (43 rows, one per dog) and `specificorders` (43 rows, one
per dog). Dog identity was restored for the 122 pooled rows by rebuilding them from
`specificorders`; the reconstruction reproduces all four columns of `all` exactly as multisets
(`habit1`, `habit3`, `habit4`, `novel`), which is what licenses the dog-level re-analysis below.
Each sheet of `Part_2.xlsx` is three column blocks of different lengths side by side, not one table:
24 rows of per-dog 5-s scores, 36 rows of
`Condition`/`Frequency`, and 72 rows of
`Condition_MF`/`Frequency_MF`/`Sex_MF`. The `Frequency` column is the mean over dogs for each
sample point — verified to machine precision for all eight outcomes — so the three blocks are three
different units of analysis, and the script picks a different one in each section.

**Coverage.** All 39 numerically reported claims are in the ledger with a verdict, and
all were re-derived from the authors' data. Six comparisons the paper describes but never reports
numerically were recovered from the script and are reported here for the first time. Still out of
reach: the figure images (not supplied, so plotted values and significance brackets in Figs 3–6 were
not checked), and the design questions in section 9.

## 2. Summary of findings

| verdict | n |
|---|---|
| supported | 11 |
| supported-with-caveat | 6 |
| assumption-dependent | 4 |
| underpowered-resolution | 0 |
| not-reproducible | 18 |
| contradicted | 0 |
| not-assessable | 0 |

Every reported statistic reproduces. `statcheck.csv` checks 39 reported
(statistic, p) pairs against recomputation: 0 decision errors,
2 p-values outside their rounding window, 1 statistic
mismatches. The three exceptions are all last-digit matters with no bearing on any conclusion:
W3 prints *V* = 2061.5 against a recomputed 2061; W2 prints
*p* = 0.008 for 0.00747; W16 prints *p* = 0.045
for 0.04679. What fails is not the arithmetic but the unit the arithmetic was done on.

**Finding 1 — Part 2's battery is a test on time bins, and most of it does not survive the dogs
(stages 2–3).** Every Kruskal–Wallis and Dunn comparison in Part 2 ran on the 12 per-condition
sample-point means, each itself a mean over 6/9/9
dogs. Re-running the same comparisons on the 24 dogs' own means:
5 of
6 significant omnibus results and
8 of 11 Dunn comparisons no longer
cross alpha, including the paper's central attachment result — proximity to owner, reported
*p* < 0.001, dog-level Kruskal–Wallis *p* = 0.231 (ANOVA
0.190, quasi-Poisson on per-dog counts 0.212). Joy
behaviours is the one outcome that survives (0.038,
0.019, 0.034). Resolution: analyse the per-dog scores
that are already in the deposit, or model the 5-s bins with dog as a random effect.

**Finding 2 — the four sex claims cite an interaction but are worded as a simple effect
(stage 3).** Each reported coefficient is the `Condition × Sex` term, a difference of differences
against the blank condition, while each sentence claims a sex difference *within* one odour
condition. For two of the four the claimed contrast is not significant:
proximity to owner under stress, -0.092, *p* = 0.193
(reported *p* = 0.003), and proximity to stranger under joy,
+0.158, *p* = 0.071 (reported
0.007). At the dog level none of the four survives: Wilcoxon rank-sum
*p* = 0.70, 0.27, 0.71, 0.37 on
4–5 dogs per cell. Resolution: report the simple effects
the sentences describe, with the dog as the unit.

**Finding 3 — Part 1 holds, but its own negative control flips (stages 2 and 5).** The pooled
Part 1 tests mix units with the rest of the section (122 trials against
43 dogs), yet the intraclass correlation of the discrimination difference is
0.002 (effective n 122 of
122) — so for the discrimination comparison pooling barely inflates precision, though
for the dishabituation difference it does (ICC 0.30, effective n
78). Habituation, discrimination and dishabituation all
survive on the 43 dog means
(*p* = 2.5e-11, 0.026, 3.9e-05). The
exception is the format control: "no dishabituation due to the change in presentation format" is
*p* = 0.087 on trials but *p* = 0.013 on dogs, with
investigation rising from the third to the fourth presentation. That control is what separates
dishabituation to the novel odour from a response to seeing two boxes. Resolution: report the
dog-level control, or rest the dishabituation claim on the novel-vs-simultaneous comparison (W2)
alone.

## 2.1 Revision 1 inferences, checked against the data

| revision-1 inference (text only) | outcome |
|---|---|
| Pooled Part 1 tests ran on ~119–122 trial-level pairs, not 43 dogs (from *V* ≤ n(n+1)/2) | **confirmed** — sheet `all` is exactly 122 trial rows with no subject column |
| Part 2 tests ran on 12 sample-point means per condition (from the H ceiling of 23.00 for 9/9/6 dogs) | **confirmed** — `Frequency` is the per-sample-point mean over dogs, 12 per condition |
| GLM residual df ≈ 46–79 for G1, inconsistent with 24 dogs or 36 means | **confirmed and pinned** — residual df = 66 on 72 sex × condition × sample-point means |
| The four sex p-values lie below the 0.0159 dog-level permutation floor | **confirmed** — dog-level tests give *p* ≥ 0.27 |
| "Unreported comparisons showed no statistical significance" | **confirmed** — the 6 recovered tests give *p* = 0.39, 0.99, 0.84, 0.71, 0.38, 0.76 |
| Posture counts fall 5 short of 24 × 12 = 288 | **refined** — 284 cells were codable (4 missing 60-s cells), and the data give upright = 170 against 169 printed, so the gap is 4 uncoded cells plus a one-count difference in the paper |
| The age dispersion "± 40.1 months" must be an SD, not an SEM | **confirmed** — `Age_months` gives mean 63.8, SD 40.1, SEM 6.1 over n = 43 |
| Pooling would inflate the evidence for the Part 1 claims | **overturned** — ICC ≈ 0.00; effective n is 122 of 122, and all three claims survive at the dog level. The unit problem that matters is in Part 2, where averaging over dogs *created* the significance |
| The missing `family` argument makes these count models Gaussian-on-counts | **overturned** — `Frequency_MF` is a mean of per-dog counts over 3–5 dogs, not a count, so the family is not the issue; the unit is |
| F2, F3, F5 and all six per-claim panels are not assessable | **overturned** — all are built from the authors' data in this revision |

## 3. Claim ledger — F4

![F4 claim ledger: reported vs recomputed p-values and recomputed effects]({{artifact:art_23a62007-ba46-4e35-9cd8-2ebd7b7d7e16}})

| claim_id | claim_short | verdict | unit_tested | p_relation | p_reported | p_recomputed | p_permutation | p_floor_design |
|---|---|---|---|---|---|---|---|---|
| W1 | 1st vs 3rd habituation presentation | supported-with-caveat | pooled trial | < | 0.001 | 9.95e-18 | 0.0001 | 3.76e-37 |
| W2 | novel vs 4th habituated (discrimination) | supported-with-caveat | pooled trial | = | 0.008 | 0.00747 | 0.0243 | 3.76e-37 |
| W3 | novel vs 3rd habituated (dishabituation) | supported-with-caveat | pooled trial | < | 0.001 | 9.97e-05 | 0.0022 | 6.02e-36 |
| W4 | 3rd vs 4th habituated (format control) | not-reproducible | pooled trial | = | 0.087 | 0.0875 | 0.138 | 1.5e-36 |
| W5 | joy-baseline: novel vs 4th | assumption-dependent | dog (pairing) | = | 0.013 | 0.0125 | 0.0053 | 1.82e-12 |
| W6 | joy-baseline: novel vs 3rd | supported | dog (pairing) | < | 0.001 | 0.000588 | 0.0006 | 7.28e-12 |
| W7 | joy habituated / baseline novel: novel vs 4th | supported | dog (order) | = | 0.004 | 0.00392 | 0.0019 | 3.81e-06 |
| W8 | baseline habituated / joy novel: novel vs 4th | supported | dog (order) | = | 0.47 | 0.473 | 0.241 | 9.54e-07 |
| W9 | joy habituated / baseline novel: novel vs 3rd | supported | dog (order) | = | 0.003 | 0.00334 | 0.0041 | 3.81e-06 |
| W10 | baseline habituated / joy novel: novel vs 3rd | assumption-dependent | dog (order) | = | 0.042 | 0.0421 | 0.0813 | 3.81e-06 |
| W11 | stress-baseline: novel vs 4th | supported | dog (pairing) | = | 0.82 | 0.818 | 0.721 | 9.09e-13 |
| W12 | stress-joy: novel vs 4th | supported | dog (pairing) | = | 0.095 | 0.0949 | 0.18 | 9.09e-13 |
| W13 | stress habituated / joy novel: novel vs 4th | supported | dog (order) | = | 0.001 | 0.00121 | 0.0019 | 1.91e-06 |
| W14 | stress-joy: novel vs 3rd | supported | dog (pairing) | < | 0.001 | 0.000727 | 0.0058 | 1.46e-11 |
| W15 | stress habituated / joy novel: novel vs 3rd | supported | dog (order) | = | 0.003 | 0.00315 | 0.0015 | 1.91e-06 |
| W16 | joy habituated / stress novel: novel vs 3rd | assumption-dependent | dog (order) | = | 0.045 | 0.0468 | 0.189 | 1.53e-05 |
| K1 | stress behaviours ~ odour condition | not-reproducible | 12 sample-point means per condition | < | 0.001 | 2.8e-05 | — | 1.53e-10 |
| K2 | joy behaviours ~ odour condition | supported-with-caveat | 12 sample-point means per condition | < | 0.001 | 6.21e-06 | — | 1.53e-10 |
| K3 | proximity to owner ~ odour condition | not-reproducible | 12 sample-point means per condition | < | 0.001 | 9.61e-06 | — | 1.53e-10 |
| K4 | proximity to odour box ~ odour condition | supported | 12 sample-point means per condition | = | 0.37 | 0.367 | — | 1.53e-10 |
| K5 | proximity to stranger ~ odour condition | not-reproducible | 12 sample-point means per condition | = | 0.016 | 0.0159 | — | 1.53e-10 |
| K6 | upright posture ~ odour condition | not-reproducible | 12 sample-point means per condition | < | 0.001 | 0.000215 | — | 1.53e-10 |
| K7 | sitting posture ~ odour condition | supported | 12 sample-point means per condition | = | 0.079 | 0.0789 | — | 1.53e-10 |
| K8 | prone posture ~ odour condition | not-reproducible | 12 sample-point means per condition | < | 0.001 | 0.000403 | — | 1.53e-10 |
| D1 | stress behaviours: blank vs joy | assumption-dependent | 12 vs 12 sample-point means | < | 0.001 | 0.000109 | 0.0587 | 0.0004 |
| D2 | stress behaviours: blank vs stress | not-reproducible | 12 vs 12 sample-point means | < | 0.001 | 1.41e-05 | 0.112 | 0.0004 |
| D3 | joy behaviours: blank vs joy | supported-with-caveat | 12 vs 12 sample-point means | < | 0.001 | 2.91e-06 | 0.0166 | 0.0004 |
| D4 | joy behaviours: blank vs stress | supported-with-caveat | 12 vs 12 sample-point means | < | 0.001 | 5.41e-05 | 0.0243 | 0.0004 |
| D5 | proximity to owner: blank vs stress | not-reproducible | 12 vs 12 sample-point means | < | 0.001 | 1.25e-06 | 0.0898 | 0.0004 |
| D6 | proximity to owner: joy vs stress | not-reproducible | 12 vs 12 sample-point means | < | 0.001 | 0.000708 | 0.217 | 4.11e-05 |
| D7 | proximity to stranger: blank vs stress | not-reproducible | 12 vs 12 sample-point means | = | 0.002 | 0.00227 | 0.18 | 0.0004 |
| D8 | upright posture: blank vs joy | not-reproducible | 12 vs 12 sample-point means | < | 0.001 | 2.74e-05 | 0.0739 | 0.0004 |
| D9 | upright posture: blank vs stress | not-reproducible | 12 vs 12 sample-point means | = | 0.004 | 0.00351 | 0.122 | 0.0004 |
| D10 | prone posture: blank vs joy | not-reproducible | 12 vs 12 sample-point means | < | 0.001 | 5.24e-05 | 0.15 | 0.0004 |
| D11 | prone posture: blank vs stress | not-reproducible | 12 vs 12 sample-point means | = | 0.005 | 0.00463 | 0.218 | 0.0004 |
| G1 | stress behaviours: F > M under stress odour | not-reproducible | 72 sex x condition x sample-point means | = | 0.007 | 0.00677 | 0.344 | 0.0159 |
| G2 | joy behaviours: F < M under stress odour | not-reproducible | 72 sex x condition x sample-point means | < | 0.001 | 2.1e-06 | 0.413 | 0.0159 |
| G3 | proximity to owner: F > M under stress odour | not-reproducible | 72 sex x condition x sample-point means | = | 0.003 | 0.00336 | 0.723 | 0.0159 |
| G4 | proximity to stranger: M > F under joy odour | not-reproducible | 72 sex x condition x sample-point means | = | 0.007 | 0.00708 | 0.531 | 0.0159 |

`p_recomputed` is the authors' own test re-run on their data; `p_permutation` is the label- or
sign-flip permutation at the unit that was randomised. Full reasons, effect sizes, bootstrap
intervals and the floor basis per row are in [`claim_ledger.csv`](claim_ledger.csv); the
statistic-by-statistic check is in [`statcheck.csv`](statcheck.csv).

## 4. Design and sampling — F1, F2

![F1 design map: observations per dog, units per pairing, effective n]({{artifact:art_447eb9ed-4efd-42ca-8867-6075b14a79ce}})

**Population.** 43 pet dogs recruited from one laboratory's owner mailing list in New York
City, screened on six owner-reported criteria including being comfortable with new people and
treat-motivated; mean age 63.8 months (SD 40.1, range
10–166), collected in three periods
(2025-11-24 00:00:00 n=17, 2025-04-24 00:00:00 n=14, Feb-Mar24 n=12). Odour stimuli came from **one** donor dog, so
every discrimination claim is a claim about discriminating that donor's three sample types. The
claims are supported for socially confident, food-motivated urban pet dogs whose owners volunteer
for cognition studies; generalising beyond that needs a sample not filtered on sociability.

**What is n, and of what.** F1 shows the Part 1 structure: 122 trials from
43 dogs (37 dogs x 3, 5 dogs x 2, 1 dogs x 1
trials). The ICC of the discrimination difference is 0.002 and of the
dishabituation difference 0.304; of the novel-odour investigation time
itself 0.371. Clustering is therefore mild for the differences the claims rest
on, which is why the pooled and dog-level answers agree for W1–W3. The reporting inconsistency
stands on its own: the three pooled tests use trial rows while every by-pairing and by-order test in
the same section uses dog rows, and the paper does not say so.

Part 2 nests three units in every sheet — 24 dogs (blank 6,
joy 9, stress 9), 36
condition-by-sample-point means, and 72 condition-by-sex-by-sample-point means
— and the script tests the second for the omnibus comparisons and the third for the sex models. The
dogs never appear as the unit. A further wrinkle: the number of dogs behind each sample-point mean
is not constant. For the 60-s bin it is
blank 5, joy 7, stress 8
instead of 6/9/9,
so the 12 values entering each test are means over different denominators and are not equally
precise.

![F2 raw data behind the owner-proximity claim: every dog, SD against SEM]({{artifact:art_927248c0-b9a8-4dba-837d-276d4f43af6d}})

**What the error bars describe.** F2 is the figure the paper does not contain. Left: every dog's
own score for the claim in H4. Right: mean ± SD beside mean ± SEM. Figs 3–6 of the paper plot
S.E.M., and Fig 6 states its unit explicitly as "means of the 12 sample points", so the published
interval describes the precision of a time-bin average, not variability between dogs — which is the
quantity the claims are about. The subject-age dispersion ("± 40.1 months") names neither quantity;
the deposit settles it as the SD (40.1 months; the SEM is
6.1).

**Internal arithmetic.** One posture is scored per coded sample point. The deposit contains
284 coded cells (4 of the 24 × 12 possible are missing at
60 s) and the posture counts total
upright 170 + sit 20 + prone
94 = 284, exactly the number of coded cells.
The Discussion prints upright 169, one fewer. The subject accounting
closes: 24 analysed in Part 2 + 2 excluded for not investigating + 3 discontinued + 14 not
offered Part 2 = 43.

**Confounders the paper tests.** Two rooms, an owner-distraction task added after the first 14
owners, and three collection periods, entered into six side-bias regressions reported only as
*p* > 0.05. Re-running them: 0 of
18 coefficients across the six models reach *p* < 0.05,
and the largest model R² is 0.097. A non-significant test is not a
demonstration of no effect, and with 17/14/12 dogs per period
these models have little resolution; the coefficients and intervals belong in the paper.

## 5. Method and assumptions — F3

![F3 distribution shape: Q-Q per condition and on the paired differences]({{artifact:art_f36e3291-cdca-4228-84db-fa9ef6f84ae3}})

Shape is described here, never graded, and no test in this paper is called invalid on the strength
of a normality diagnostic.

**The route was chosen by a normality test, and that choice is invisible in the paper.**
`statistical_analyses.R` runs `shapiro.test` on all eight Part 2 outcomes and branches on the
result: `box` (*W* = 0.961, *p* = 0.229) goes to
`aov`, the other seven to `kruskal.test`. Two of those seven sit just below the threshold — owner
*p* = 0.046 and stranger *p* = 0.040
— so which route those outcomes took was decided by a diagnostic within 0.01 of its cut-off. Two
observations about the diagnostic itself, both descriptive: it is computed on the 36 values pooled
across the three conditions, so it describes a three-group mixture rather than any group or the
model residuals, and a mixture of three groups with different locations departs from normality
*because* the groups differ; and at n = 36 the test has modest power either way. The consequence
worth reporting is procedural: a reader of the paper cannot tell which comparisons were parametric.

**The Dunn p-values are halved and unadjusted.** `dunn.test` prints one-tailed values by default,
and those are the numbers in the paper: for proximity to stranger the paper reports
*p* = 0.002 against a halved value of 0.0023 and a
two-sided value of 0.0045; for prone posture, 0.005
against 0.0046 and 0.0093. All eleven reported Dunn values
match the halved convention and none matches the two-sided one. No decision changes — the largest
two-sided value is 0.009 — and no adjustment was applied
either, though `dunn.test` can apply one. Resolution: state the convention, or report two-sided
adjusted values.

**Interaction against simple effect.** The sex models are
`glm(Frequency_MF ~ Condition_MF * Sex_MF)` with no family argument, i.e. Gaussian least squares.
The outcome is a mean of per-dog counts over 3–5 dogs, so it is not a count and the missing family
argument is not the defect; the unit is, and the contrast is. Each reported coefficient is an
interaction term; the matching within-condition contrasts are
G1 -0.308 (p = 4.9e-05); G2 +0.429 (p = 0.00135); G3 -0.092 (p = 0.193); G4 +0.158 (p = 0.0707).
G1 and G2 agree with their sentences; G3 and G4 do not.

**The three routes.** Every claim now carries parametric, rank-based and permutation routes at the
unit that was randomised; they appear per row in the ledger and on each per-claim panel. They agree
for most claims. Where they disagree the verdict is `assumption-dependent` and the disagreement is
named:
W5 (joy-baseline: novel vs 4th); W10 (baseline habituated / joy novel: novel vs 3rd); W16 (joy habituated / stress novel: novel vs 3rd); D1 (stress behaviours: blank vs joy).
The pattern in those rows is the one worth stating: the signed-rank and t routes reach alpha while
the **sign test does not**, because the effect is carried by a large difference in a minority of
dogs rather than by most dogs moving the same way. For W5, 26
of 40 dogs investigated the novel odour longer (sign test
*p* = 0.081) against Wilcoxon *p* = 0.0125. "Dogs investigate the
novel odour longer" and "most dogs investigate the novel odour longer" are different claims, and the
data support the first more strongly than the second.

## 6. Resolution and evidence — F5

![F5 permutation null and bootstrap interval for the pooled discrimination claim]({{artifact:art_e55f1f54-b7d7-45f6-a19a-1b869331a425}})

| claim | design | n | smallest_possible_p | reported_p |
|---|---|---|---|---|
| W2 | two-sided sign-flip floor 2/2^n | 43 | 3.76e-37 | 0.008 |
| W13 | two-sided sign-flip floor 2/2^n | 20 | 1.91e-06 | 0.001 |
| D5 | two-sided label-permutation floor 2/C(15,6) | 15 | 0.0004 | 0.001 |
| D7 | two-sided label-permutation floor 2/C(15,6) | 15 | 0.0004 | 0.002 |
| G1 | two-sided label-permutation floor 2/C(9,4) | 9 | 0.0159 | 0.007 |

No design floor in this paper exceeds α = 0.05, so no claim is `underpowered-resolution`: each
design could in principle have produced a significant result. Two resolution facts still shape how
the numbers read. The four sex-by-condition p-values are smaller than the smallest value a
relabelling of the dogs in that condition could produce (0.0159 for
4 vs 5 dogs), which is the arithmetic signature of a
p-value that comes from the model and the time bins rather than from the animals. And every Dunn
comparison involving the blank condition rests on 6 dogs, whose
dog-level floor is 0.0004.

**Effect sizes and intervals, which the paper does not report.** The pooled discrimination effect is
+1.00 s per dog (95% BCa CI
[0.23, 1.98], 43 dogs,
10,000 resamples), matched-pairs rank-biserial +0.28 at trial level.
Dishabituation is +1.77 s
[0.63, 2.79]. The
owner-proximity difference behind H4 is +0.37 score points between
6 blank and 9 stress dogs. F5 shows the permutation null and the
bootstrap distribution for the discrimination claim side by side; the permutation resolution is
1/10,001 = 1.0e-04, so no p-value from it should be read below that.

## 7. Multiplicity — F6

![F6 multiplicity: the full p-value inventory and hits under each correction]({{artifact:art_01999879-67e6-47dd-ac83-d393c1ddcd9f}})

| analysis | tests_run | reported_in_paper | basis |
|---|---|---|---|
| Wilcoxon signed-rank (Part 1) | 22 | 16 | 3 pooled + 6 by pairing + 12 by presentation order + 1 format control |
| Side-bias / season regression coefficients | 18 | 0 | 6 lm models x (side + 2 season contrasts); reported only as p > 0.05 |
| Shapiro-Wilk route selection | 8 | 0 | one per Part 2 outcome; decides aov vs kruskal.test |
| Part 2 omnibus | 8 | 8 | 7 kruskal.test + 1 aov |
| Dunn pairwise comparisons | 18 | 11 | 6 dunn.test families x 3 pairwise comparisons |
| GLM coefficients | 25 | 4 | 5 models x (2 condition + 1 sex + 2 interaction) coefficients |

The script produces **99 p-values**; the paper reports 39. At α = 0.05 over
99 tests about 5.0 significant results are expected by chance alone. Recounting
the 99 recomputed values that are defined:

| correction | controls | n_significant | threshold_on_raw_p |
|---|---|---|---|
| none | nothing (per-test error rate) | 53 | 0.05 |
| Bonferroni | FWER: probability of any false positive | 22 | 0.0005051 |
| Benjamini-Hochberg | FDR: expected fraction of false calls <= 0.05 | 46 | 0.01647 |

*99 tests at alpha=0.05 yield about 5.0 significant results by chance alone. 53 uncorrected hits become 46 under BH and 22 under Bonferroni.*

Within families, under Benjamini–Hochberg:

| family | tests | survive_bh |
|---|---|---|
| Part 1 order-level (12 tests) | 12 | 4 |
| Part 2 omnibus (8 tests) | 8 | 6 |
| Dunn (18 comparisons) | 18 | 11 |
| GLM coefficients (25) | 25 | 14 |

The paper states no correction anywhere. Which claims need one is a judgement about the claim: the
three pooled Part 1 comparisons are pre-specified by the paradigm and read as primary, while the
Part 2 battery — two behaviour categories, three proximity targets, three postures, then a sex
interaction on five of them — is a screen, and "dogs exhibited distinct behavioural responses
depending on the emotional content of the odor" is a claim about that family.

**Distribution of the p-values.** 2 of 99 fall in
[0.045, 0.05) against 2 in [0.04, 0.045) (one-sided binomial
*p* = 0.688; the function reports this as
interpretable at this scale). Three limits
apply and apply here: the diagnostic is a property of a *collection* of p-values and says nothing
about any individual result; it needs a reasonable number of them; and it cannot distinguish
selective reporting from legitimate exploratory analysis. **No intent is attributed here, and none
can be.** The inventory above — the full list of tests run, now visible because the script was
deposited — is what a reader needs, and depositing it is what made this recount possible.

## 8. Per-claim supplement

One panel per headline claim: the raw observations at the unit that was randomised, the permutation
null, and the three routes with the design floor shaded.

![H1 dogs investigate a novel odour longer than a repeated one]({{artifact:art_04a335cc-3936-4083-bdfa-a0695b768952}})

![H2 dogs distinguish joy from baseline odour]({{artifact:art_83df6f59-37cc-4fc5-93f5-57218007cccf}})

![H3 dogs distinguish joy from stress, not stress from baseline]({{artifact:art_4c78d884-08db-4425-800e-efffde35920c}})

![H4 stress odour moves dogs toward the owner and away from the stranger]({{artifact:art_910ceac8-424f-4339-ad41-bc942e74dc4a}})

![H5 conspecific odours shift behaviour and posture relative to blank]({{artifact:art_58116637-55cf-44c9-8daa-2c752f25bd4e}})

![H6 male and female dogs respond differently to the odour conditions]({{artifact:art_9399168a-c8ec-4a07-ae1a-9411c6c2d9a5}})


## 9. What was not assessed

Diagnostics in this report describe distribution shape and spread only. They cannot evaluate:
  - independence of observations
  - the sample is representative of the stated population
  - absence of unmeasured confounding
  - measurements were not selected after seeing the outcome
These are design questions. Assess them from the methods section, and say so explicitly when the methods do not contain enough information to assess them.

Specific to this review:

- **The figure images.** Figs 1–6 were not supplied, so no plotted value, error bar or significance
  bracket was checked against the text or the data.
- **Emotion induction.** Whether the three collection events produced the intended emotional states
  in the donor dog, and whether the samples differ chemically, are not statistical questions. The
  paper's validation is the owner's behavioural judgement of one donor.
- **Blinding and randomisation as executed.** Both are described in the methods (coder blind to
  condition, randomised pairing order and side); nothing in the deposit records the randomisation
  sequence, so the description was taken as given.
- **Video coding reliability.** No second coder, and no inter-rater agreement statistic, is reported
  or deposited; the per-dog scores were taken as given.
- **Why the 60-s sample point is missing for 4 dogs**, and whether missingness
  relates to the outcome. The README attributes "null" cells to early trial termination.
- **Nothing was subsampled or approximated.** Every claim was recomputed on the full deposited data.
  Where a reported p-value is a bound ("p < 0.001"), the recomputed value is given exactly and the
  bound is kept in `p_relation`.

## 10. Reproduction

| item | value |
|---|---|
| step 1 | `Rscript handoff/recompute_authors_tests.R` — the authors' own routines on their data, plus the dog-level and count-family alternatives |
| step 2 | `python review_wang_horowitz_2026.py` — routes, floors, intervals, figures, ledger, this report |
| seeds | `SEED = 0`, `n_perm = 10,000`, `n_boot = 10,000` (BCa) |
| permutation scheme | sign-flip within dog for paired Part 1 claims; label permutation between conditions for Part 2 claims; statistic = difference in means, two-sided |
| Part_1.xlsx sha256 | `3ef960171f504113ed0aab0e5841250e07dfaf42cc46c350f7fa2ee709f20d90` |
| Part_2.xlsx sha256 | `fb594cabc647751d5b4a9944a20700230785ec24a3a1dfa466718edcc98eee1e` |
| statistical_analyses.R sha256 | `bf63d0a9271a09a8b256c69a24e65571d9c304875b26504a406f2bc64317366f` |
| paper text sha256 | `eda57f3291c4e12793af05ae9e0d390f4d3fd07368039bf46de41822c134f9ab` |
| R | R 4.5.3 |
| python | 3.11.16 |
| numpy / pandas / scipy / matplotlib | 2.4.6 / 2.3.3 / 1.17.1 / 3.11.2 |

One defect in the deposited script, verified against the columns actually present:
`glm(Frequenc_MFy ~ Condition_MF * Sex_MF, data = box)` misspells the outcome
(present in the file; `box` contains `Frequency_MF`, not
`Frequenc_MFy`). The line immediately after it summarises the correctly spelled model, so no
reported number depends on it, but `source()`-ing the script stops there. Reproduction also needs
`install.packages` lines removed and the `dunn.test` dependency available; this review recomputes
Dunn's z directly instead, which is why both p-value conventions could be compared.

Derivations used above, each reproducible by hand: *V* ≤ n(n+1)/2 for a signed-rank statistic over
n non-zero pairs; the exact two-sided sign-flip floor 2/2ⁿ; the two-sided label-permutation floor
2/C(n₁+n₂, n₁); design effect 1 + (k−1)·ICC and effective n = rows / design effect.
