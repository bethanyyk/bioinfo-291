# Statistical Review — *Lonardo, Volter & Huber (2026), Cognition 274:106586*

*Reviewer: Statistical Reviewer (automated audit) · Date: 2026-10-19 · Protocol: statistical-review v1*

**Paper.** Lonardo L., Volter C.J. & Huber L. (2026) Pupil size changes reveal that dogs are
sensitive to the social conventions behind human gestures. *Cognition* **274**:106586.
doi:10.1016/j.cognition.2026.106586

---

## 1. Inputs and coverage

**This review was done in two passes.** The first pass had only the abstract (Elsevier paywall);
the published PDF was supplied afterwards, and this version is the revision against it. The first
pass's limitations are gone: statcheck now runs, `p_reported` is populated across the ledger, and
F4 and F6 have both axes.

**What the full text confirmed, and what it overturned.** Two of the five abstract-based readings
did not survive contact with the paper, and both errors were the reviewer's:

| First-pass reading | Full text | Outcome |
|--------------------|-----------|---------|
| **C1** — the 28-vs-39 sample split is an undisclosed dependence; `supported-with-caveat` | Methods: "we only retained dogs with four valid test trials, one per condition (N = 28)"; Results report the N = 39 re-analysis and its null (t = 0.54, p = 0.592) | **Overturned.** No disclosure gap. C1 is now a reproduction result and the verdict is `supported` |
| **C5** — "context-specific" rests on a difference of significance between the two gestures; `assumption-dependent` | The paper never splits the pupil analysis by gesture. "Context-relative" is an argument from the design: the same action is congruent after one gesture and incongruent after the other, so the action cancels in the trial-type contrast | **Overturned.** The reviewer read a per-gesture significance split into a sentence that restates the design. Verdict is now `supported-with-caveat` |
| **C2, C3** — the nulls reproduce but the intervals do not establish equivalence | Reported as chi2 = 0.12, df = 1, p = 0.731 and chi2 < 0.01, df = 1, p = 0.951; no equivalence test, no smallest effect of interest | **Confirmed**, and now evidenced from the text rather than inferred |
| **C6, C7** — argued from design and offered as discussion, with no test to audit | Confirmed verbatim in the discussion and hypothesis sections | **Confirmed**; both stay `not-assessable` |
| **C4** — N = 39 reproduces; the per-analysis n differ | The paper states n separately for every analysis, and its demographics reproduce exactly | **Mostly confirmed**; one trial-accounting gap remains (section 4) |

**How well the blind recomputation did.** The first pass derived its numbers from the raw data
without seeing any of the paper's. The subject-level Wilcoxon reproduces the paper's reported
statistic exactly — reported T+ = 313, N = 28,
p = 0.011 against recomputed T+ = 313,
N = 28, p = 0.0110 — and the cluster-based permutation test,
run as a stand-in for a GAMM that could not be refitted, placed the difference at
7800–8700 ms against the paper's reported significant window of
7697–9111 ms. The reviewer's errors were in reading the
prose, not in the arithmetic.

| Input | File | Size | sha256 (first 16) |
|-------|------|------|-------------------|
| pupil sample report (ds) | `dog_nonverbal_communication_samplereport_ds.parquet` | 5781 x 13 | `f046f8539c420125` |
| pupil sample report (ALL_DATA) | `dog_nonverbal_communication_samplereport_ds_ALL_DATA.parquet` | 6232 x 13 | `af83a389a31c34b0` |
| IA report, reaction IP | `IA_report_IP_from_reaction.txt` | 618 x 143 | `8e0615cf2227711f` |
| IA report, reaction IP from frame 579 | `IA_report_IP_from_reaction_frame579.txt` | 618 x 143 | `4357ed482bbd456d` |
| IA report, whole-video IP (unused by both scripts) | `IA_report_new_video_IP.txt` | 618 x 143 | `534052ba99974182` |
| counterbalancing / demographics | `Non_verbal_communication_counterbalancing.csv` | 45 x 24 | `472d4e27b5acff5a` |
| pupil analysis code | `non_verbal_communication_pupil_size_analysis.Rmd` | R Markdown | `5edae5358442ff27` |
| dwell-time analysis code | `non_verbal_communication_dwell_time_analysis.rmd` | R Markdown | `dc40a485f6e26288` |

| Input | File | Size | sha256 (first 16) |
|-------|------|------|-------------------|
| published article | `1-s2.0-S0010027726001538-main.pdf` | 8 pages, text layer usable (True) | `d947171bf22cb83a` |

| | |
|--|--|
| Claims reviewed | 7 (C1–C7) |
| Extracted automatically from the PDF | 6 (statistic, p) pairs, 14 p-values, 9 sample-size statements, 2 intervals |
| Added by hand | 6 chi-square triples the extractor's regex misses (Greek glyph), transcribed from the sentences its p-value pass did capture; and the 7 claims themselves, since the paper states them in prose |
| Claims the data cannot address | C6, C7 — `not-assessable` |
| Seeds | `random_state = 20260101`, `n_perm = 10000`, `n_boot = 10000`, alpha = 0.05 |

**Coverage, stated honestly.** This reviews the seven claims the abstract makes, now checked
against the full text and the deposited data. 18 inferential tests are visible in the
authors' code and 14 p-values are reported in the main text; this review re-derives
five of them from raw data and arithmetic-checks 6 of the reported triples. What
remains outside it: **Tables S1–S4 and the supplementary figures were not supplied**, so the GAMM
and GLMM coefficient tables, their degrees of freedom, and the 1000-sample bootstrap CIs the
authors computed are unseen — that is why 6 of the extracted test statistics could
not be arithmetic-checked (section 5). The authors' three model classes could not be refitted
here, so the published fits themselves are not reproduced; this review audits the contrasts those
fits estimate.

**Which data correspond to which claim.**

| Claim | Supporting file |
|-------|-----------------|
| C1, C5 (pupil size) | `..._samplereport_ds.parquet` (5781 rows) and `..._ALL_DATA.parquet` (6232 rows) — 41 time bins per trial over 1000–5000 ms of the baseline+reaction interest period (video time 6000–10000 ms) |
| C2 (looking time) | `IA_report_IP_from_reaction.txt` plus `..._frame579.txt` for Zserbo session 1 — 304 agent interest-area rows from 152 sessions |
| C3 (gaze transitions) | the same IA reports, via `INTEREST_AREA_FIXATION_SEQUENCE` — 152 trial-level transition counts |
| C4 (N = 39) | both parquet files plus `Non_verbal_communication_counterbalancing.csv` |
| C6, C7 | no corresponding data |

`IA_report_new_video_IP.txt` covers the whole-video interest period and is used by neither of the
authors' scripts; it is not used here either.

**The models the authors actually specified**, read from the deposited `.Rmd` files. The task
brief described the analysis as "GLMMs with bootstrapped CIs"; that is accurate for looking time
and gaze transitions but not for the headline pupil analysis:

- **Pupil size (headline): a GAMM, not a GLMM** — `bam(pupil_base ~ trial_type + session_number +
  s(time, k=20) + s(time, by=trial_type, k=20) + s(Xgaze_RIGHT, Ygaze_RIGHT) +
  s(time, subject, bs='fs', m=1) + s(time, Event, bs='fs', m=1), method="ML")`, fitted twice: on
  complete cases (`dat3`, 28 dogs) and on `ALL_DATA`
  (39 dogs), with inference read off the `plot_diff` difference
  curve. The same script additionally runs a subject-aggregated paired t-test, a Shapiro-Wilk test
  on the difference scores, and an exact Wilcoxon signed-rank test on
  28 dogs.
- **Looking time:** `glmmTMB(prop.dwell.time.agents.scaled ~ trial_type + agent_type + z.trial +
  (1 + trial_type.c + agent_type.c + z.trial || subject), family = beta_family)`, evaluated by
  `drop1(test="Chisq")`, with 1000 parametric bootstrap CIs.
- **Gaze transitions:** `glmer(transition_count ~ trial_type + z.trial +
  offset(log(TRIAL_DWELL_TIME)) + (1 + z.trial | subject), family = poisson)`, `drop1(test="Chisq")`,
  1000 bootstrap CIs.

None of these three model classes is available in this environment, so **this review does not
reproduce the authors' models.** It re-derives the same contrasts with dog as the unit of analysis
and, for pupil size, adds a time-resolved cluster-based permutation test, which targets the same
estimand as the GAMM difference curve — a difference localised in time — without assuming a
parametric smooth.

---

## 2. Summary of findings

| Verdict | Claims |
|---------|--------|
| supported | 1 |
| supported-with-caveat | 4 |
| assumption-dependent | 0 |
| underpowered-resolution | 0 |
| not-reproducible | 0 |
| contradicted | 0 |
| not-assessable | 2 |

1. **C1 reproduces exactly, and the first pass's caveat against it was the reviewer's error.** The
   paper reports T+ = 313, N = 28,
   p = 0.011; recomputing from the raw parquet gives
   T+ = 313, N = 28, p = 0.0110, with all
   four routes agreeing (paired t 0.0204, Wilcoxon 0.0110, sign
   0.0357, permutation 0.0119) and a blind cluster permutation landing at
   7800–8700 ms inside the paper's reported 7697–9111 ms
   window. The 28-vs-39 split the first pass flagged is stated in the methods and the N = 39 null
   is reported by the authors themselves (t = 0.54,
   p = 0.592); verdict moves to `supported` (stages 4–5).
2. **One to two excluded trials are unaccounted for between the methods and the deposit.** The
   methods state 12 trials were excluded from the pupil analyses
   (4 for pupil detection below 70% of the interest period,
   8 for on-screen gaze below 70%); the deposited files differ by
   11 trials — `ALL_DATA` holds 152 and the analysis file
   141 — and only 10 of those carried pupil data to exclude,
   since Ace/CA is empty across all 41 bins. No verdict depends on it:
   all 11 belong to dogs already outside the 28-dog sample. Everything else in the
   accounting reconciles to the digit: 39 subjects,
   18 females against 18 reported, mean age
   59.33 against 59.3 months, and the
   4 missing sessions exactly as described (stage 1).
3. **Both nulls reproduce; neither interval establishes equivalence, and no correction is applied
   across 14 reported p-values.** Looking time differs by -3.8 ms
   (95% CI -459 to 419, on a congruent mean of
   4189) and transitions by 0.066
   (95% CI -0.38 to 0.46, mean 2.20) — roughly
   ±11% and
   ±21% of their means remain
   compatible with the data, and the paper states no smallest effect of interest to test against
   (stages 5–6).

---

## 3. Claim ledger — F4

![F4 — claim ledger]({{artifact:art_4a7f728b-43b5-4162-ace7-c34442e21d7b}})

*Recomputed effect with its 95% bootstrap CI (left) and reported against recomputed p-values on
log axes with alpha marked (right). Effects are signed incongruent − congruent, so a positive
pupil effect is the direction the paper claims. C1, C2 and C3 now carry a reported p-value; C5 does
not, because the paper reports no test of the quantity that claim asserts (section 5), and C4's
"effect" is a sample size, which is held out of the effect axis and verified in section 4 instead.*

| id | claim | verdict | effect rep. | effect recomp. | 95% CI | p rep. | p recomp. |
|----|-------|---------|-------------|----------------|--------|--------|-----------|
| C1 | pupils dilated more to incongruent responses | `supported` | - | 141.6 | 51.6 to 283 | 0.011 | 0.0119 |
| C2 | looking times did not differ cong vs incong | `supported-with-caveat` | - | -0.007462 | -0.0414 to 0.0264 | 0.731 | 0.6779 |
| C3 | gaze transitions did not differ cong vs incong | `supported-with-caveat` | - | 0.06579 | -0.382 to 0.461 | 0.951 | 0.771 |
| C4 | eye-tracking study with N = 39 dogs | `supported-with-caveat` | 39 | 39 | - | - | - |
| C5 | effect is context-specific (action-balanced design) | `supported-with-caveat` | - | -76.32 | -286 to 178 | - | 0.5237 |
| C6 | expectations are not self-centred (third party) | `not-assessable` | - | - | - | - | - |
| C7 | expectations form implicitly during ontogeny | `not-assessable` | - | - | - | - | - |

The full `reasons` text for every row is in `claim_ledger.csv`.

---

## 4. Design and sampling — stages 1–2

![F1 — design map]({{artifact:art_853987d5-7c56-4177-a900-faa17a6dd1e4}})

*Trials per dog and dogs per trial type, from `ALL_DATA`. The n available to a test is dogs — not
trials, and not the 6232 time-bin rows. The panel counts
151 trials rather than the 152 present in the file
because one trial (Ace/CA) has no pupil value in any of its 41 bins
and drops out of a figure keyed on that value; section 4(d) follows that trial through the
accounting.*

![F2 — raw data behind the claim]({{artifact:art_47ed4615-795a-4a14-80d9-a990f7f497e2}})

*Each dog's window-mean pupil size by trial type, with SD and SEM shown on one axis. The paper's
own error-bar definition could not be checked, because its text was not supplied.*

**Population.** The abstract claims a property of pet dogs in general ("pet dogs, during their
ontogeny"). The deposited data describe 39 dogs tested in a laboratory eye-tracker,
18 of them female, carrying 13 distinct breed labels. Whether
that sample represents pet dogs as a population is a methods question this review cannot answer:
recruitment frame, screening and the chin-rest training requirement are described only in the
paper. Chin-rest eye-tracking selects for trainable, cooperative dogs, and nothing in the data
bounds that selection.

**Unit of replication.** The pupil files hold 5781 and 6232 rows, but
these are 41 time bins × up to 4 trials × dog, so the row count is not an n. Two separate things
determine what n a test may use, and they point the same way here for different reasons.

The first is clustering, which turns out to be slight: across the 152 trial-level
pupil means from 39 dogs, the ICC within dog is 0.000
with a design effect of 1.00, giving an effective n of
151.0 — essentially the row count. A dog's four trial means are, on this
summary, about as informative as four independent observations, so a variance-inflation argument
against pooling trials is weak. (Within a trial, the 41 time bins are a different matter: they are
a smoothed time series, which is why the authors model them with `s(time, Event)` rather than
counting them.)

The second is the design, and it is the binding one: congruency was assigned **within** dog, so the
comparison each claim makes is a paired one and the number of independent contrasts is the number of
dogs — 39 at most, fewer wherever a dog is missing a condition. A near-zero
ICC does not license treating the trials as independent replicates of the contrast, because it is
the pairing, not the correlation, that defines the unit. Every test in this review therefore takes
dog as the unit by construction; the authors' GAMM reaches the same place through its
`s(time, subject)` and `s(time, Event)` random smooths.

**Who actually contributes.** The review's first substantive finding concerns a count, not a
p-value:

| Analysis | Dogs contributing |
|----------|-------------------|
| Abstract's stated N | 39 |
| Distinct subjects in either parquet file | 39 |
| Named dogs in the counterbalancing file | 39 (of 45 rows; 6 blank) |
| Authors' primary pupil GAMM and aggregated tests (`dat3`) | 28 |
| Pupil contrast available at subject level, `ds` file | 37 |
| Pupil contrast available at subject level, `ALL_DATA` | 38 |
| Looking time / transitions | 39 dogs, 152 trials; 38 with a paired contrast |

Three things follow.

**(a) N = 39 reproduces exactly** — 39 distinct subjects in both parquet files, 39 named dogs in
the counterbalancing file.

**(b) The counterbalancing file's 45 rows and 30 unique
`Dog_ID` values are not a subject-count discrepancy**, and this review checked rather than assumed
it. 6 rows are blank, and every `Dog_ID` appearing more than once maps to
exactly one counterbalancing `order` string (`order` unique within `Dog_ID`:
True, over 7 shared ids; `Dog_ID` 18, for
instance, is shared by three dogs that all ran `GG-GA-CA-CG`). `Dog_ID` indexes the
counterbalancing slot; `Name` is the subject key. The apparent 45/30 mismatch is a column-semantics
question, not a missing-dogs question, and it is resolved by that one sentence.

**(c) The analyses do not all run on 39.** Eleven dogs (Ace, Filou, Floki, George, Ignaz, Jimi, Lenny, Sahibu, Sheila, Timo, Zserbo) lack at least one of the
four conditions. The authors exclude exactly those eleven from the primary pupil analysis, leaving
28; this review confirms the hard-coded exclusion list in their script is
precisely the set of dogs with fewer than four conditions in the data, so the rule is principled
and reproducible. **The paper states this**: "To ensure a within-subject comparison, we only
retained dogs with four valid test trials, one per condition (N = 28)", and it reports the N = 39
re-analysis and its null result as well. The abstract's "N = 39" is the recruited and eye-tracked
sample and the pupil result rests on 28 dogs — a reader of the abstract alone could conflate
them, as the first pass of this review did, but the paper itself does not.

**(d) Trial accounting, which is where the one remaining gap is.** The two pupil files nest
exactly — every trial in the analysis file is also in `ALL_DATA` (`subset = True`)
— so the chain can be followed step by step:

| step | trials | source |
|------|--------|--------|
| intended: 39 dogs × 4 conditions | 156 | design |
| less sessions never recorded | −4 | methods: "one dog had 3 missing sessions – CC, CG, GC; the other dog had 1 missing session - GC" |
| present in `ALL_DATA` | 152 | 152 × 41 time bins = 6232 rows, its exact row count |
| of which carry no pupil data at all | 1 (Ace/CA: all 41 bins missing) | data |
| present in the analysis file | 141 | 141 × 41 = 5781 rows |
| removed between the two files | 11 | set difference |
| of those, carrying usable pupil data | 10 | data |
| excluded for quality per the methods | 12 | 4 pupil detected under 70% of the IP + 8 on-screen gaze under 70% |

Everything above the last line reconciles. The last line does not: the methods describe
12 quality exclusions, while 11 trials are absent
from the analysis file, of which only 10 carried any pupil data to
exclude. **So 1 trial is unaccounted for on the loosest reading and
2 on the strictest**, depending on whether the empty trial is counted
among the twelve. The removed trials are
Ace/CA, Filou/CA, Filou/CG, Filou/GG, George/CA, Ignaz/CG, Jimi/CG, Sahibu/GG, Sheila/CG, Timo/CA, Zserbo/CG.

This touches no verdict: all 11 belong to dogs already outside the 28-dog
complete-case sample, which is defined by missing conditions rather than by these exclusions. It
is the kind of discrepancy an author resolves in one sentence, and it is recorded here only
because a reader reconstructing the sample from the deposit will hit it.

A smaller point, now resolved rather than raised: the demographics join in the deposited
dwell-time script recovers the subject name by splitting a truncated EyeLink filename (`EDF3`) and
fails for 2 dogs (Floki, Knopfi). The
first pass flagged this as a possible contaminant of the reported age and sex summaries. It is
not: computed over all 39 named dogs the file gives 18 females and a
mean age of 59.33 months, against 18 and
59.3 months in the paper — both reproduce — so the published demographics
were evidently not taken through that join.

One naming note for anyone matching the paper to the deposit: the paper labels the four conditions
CC, CG, GC, GG while the data files use CA, CG, GA, GG (A for "approach" where the paper writes C
for "come closer"). The mapping is one-to-one and the congruency coding agrees — CA/CC and GG are
congruent, CG and GA/GC incongruent — but the strings differ.

**Uncertainty reporting — now checkable, and it checks out.** The paper prints
"5105 ± 1924 ms" and
"3218 ± 2065 ms" without naming what
the ± is. Recomputing the three candidates from the subject means settles it: for the reacting
agent the SD across dogs is 1924.3, the SEM is
308.1 and the 95% CI half-width is
623.8, so the printed value matches **SD**
(ratio 0.9999, unambiguous:
True); the communicating agent matches SD on the same test
(ratio 0.9999). SD is the right choice for a
descriptive statement about how dogs vary, and it is the conservative one — SEM would have been
6.24× narrower. The only reporting improvement left here
is to name it in the text. F2 draws SD and SEM on one axis so the distinction stays visible.

---

## 5. Method and assumptions — stages 3–4

![F3 — distribution shape]({{artifact:art_aa6b8775-dbdc-45ce-9bbb-941ff6a50c87}})

*Q-Q of the paired differences behind C1 with n, skewness, excess kurtosis and Shapiro-Wilk. This
is descriptive. Shape is not a gate, and nothing on this panel licenses calling any test invalid.*

**statcheck — the paper's own arithmetic.** 12 test statistics were recovered from
the text (6 by the automatic extractor, 6 chi-square triples transcribed by hand because
the extractor's regex does not match the Greek glyph). Of these, **6 carry the
degrees of freedom needed to recompute p and were checked; 6 could not be** — the
GAMM's t and F statistics are reported in the text without their df, which live in Tables S2–S4,
and the supplementary file was not supplied.

Of the 6 checkable triples, statcheck flags **3** and finds
**0 decision errors** (a decision error being a
reported p on the other side of alpha from the recomputed one). Every flagged row resolves once the
precision of the *statistic* is propagated, which the review computed rather than assumed: statcheck
recomputes p from the statistic as printed, so a chi-square given to two decimals pins p far more
tightly than the paper actually did —

- dwell proportion GLMM: trial type x agent type interaction: chi2 = 0.05, df = 1, reported p = 0.818; p at the printed value is 0.8231, but across the statistic's two-decimal rounding window the implied p runs 0.8146 to 0.8320, which contains the reported value (True)
- dwell proportion GLMM: trial type (reduced model): chi2 = 0.12, df = 1, reported p = 0.731; p at the printed value is 0.7290, but across the statistic's two-decimal rounding window the implied p runs 0.7237 to 0.7345, which contains the reported value (True)
- transition count GLMM: trial type (statistic reported as chi2 < 0.01): chi2 = 0.01, df = 1, reported p = 0.951; p at the printed value is 0.9203, but across the statistic's stated bound [0, 0.01] the implied p runs 0.9203 to 1.0000, which contains the reported value (True)

So **3 of the 3 flags are artefacts of the tool reading printed
precision as exact, and no genuine arithmetic inconsistency was found** in the statistics this
review could check. The 6 unchecked statistics are a gap in coverage, not a finding
about them: supplying Tables S2–S4 would close it.

One transcription decision is worth stating, because getting it wrong would have manufactured a
clean bill of health. For the gaze-transition trial-type test the paper prints "chi2 < 0.01,
df = 1, p = .951": the inequality bounds the *statistic*, while p is reported exactly. Coding that
row's `p_relation` as `<` would route it into statcheck's one-sided branch, which only asks whether
the recomputed p falls below the reported one — a near-automatic pass that would have reported the
row as consistent without ever comparing 0.9203
to 0.951. It is coded `=` here, which flags it, and the flag is then resolved on the statistic's
stated bound rather than waved through: any chi-square in [0, 0.01] implies p in
[0.9203,
1.0000], and .951 corresponds to chi2 ≈ 0.0038,
comfortably inside the bound the paper gave.

**Test selection.** On every axis checkable without the methods section the tests match the design:
the comparison is within-dog, every test used is paired or carries a per-subject random effect, and
all are two-sided. The substantive point concerns what a null means rather than which test was
chosen. The authors' pupil script runs a Shapiro-Wilk test on the difference scores, annotates the
result "Assumptions violated? --> Yes", and switches to an exact Wilcoxon signed-rank test. This
review reports the same diagnostic descriptively — Shapiro-Wilk p = 0.0027
on n = 28 paired differences, skewness 1.621, excess kurtosis
5.441, median difference 144.4 a.u.
against a mean of 141.6 — and does not treat it as a pass/fail gate.
The shape is right-skewed with heavy tails relative to a normal; that is a description, not a
verdict on any test. At this n the normality test has limited power, so a significant result need
not mean a consequential departure and a non-significant one would not have certified normality.
The reason to report all routes is not that a diagnostic licensed one of them, but that agreement
across them shows the conclusion does not depend on the choice. Note also that the Wilcoxon
signed-rank null is symmetry of the differences about zero rather than equality of means, so the
switch changes the hypothesis slightly; here it changes neither the direction nor the decision.
The skill's own note on this n is carried verbatim rather than paraphrased:
*"n=28 < 30: a normality test here has little power, so a non-significant result is not evidence of normality. Read the Q-Q plot for shape, and prefer a method that does not lean on the shape assumption."*

**Three routes — C1, pupil size, 28 dogs** (incongruent − congruent = 141.6 a.u.):

| route | test | p | assumes | smallest achievable p | reachable at 0.05 |
|-------|------|---|---------|-----------------------|-------------------|
| parametric | paired t-test | 0.0204 | differences are normal; pairs independent | 0.00e+00 | yes |
| non-parametric (ranks) | Wilcoxon signed-rank | 0.0110 | differences are symmetric about the null median | 0.00e+00 | yes |
| non-parametric (signs) | sign test | 0.0357 | pairs independent; nothing about shape | 7.45e-09 | yes |
| permutation | sign-flip permutation (mean difference) | 0.0119 | within-pair exchangeability under the null | 7.45e-09 | yes |

All routes agree at alpha = 0.05. The same four routes across every contrast the data support:

| contrast | n dogs | effect (incong − cong) | 95% bootstrap CI | paired t | Wilcoxon | sign | permutation | routes agree |
|----------|--------|------------------------|------------------|----------|----------|------|-------------|--------------|
| pupil, primary complete cases (dat3) | 28 | 141.6 | 51.64 to 282.6 | 0.0204 | 0.0110 | 0.0357 | 0.0119 | yes |
| pupil, ds file, no exclusion | 37 | 81.3 | -32.6 to 199.8 | 0.1807 | 0.1211 | 0.1877 | 0.1866 | yes |
| pupil, ALL_DATA (abstract's N = 39) | 38 | 55.58 | -107.6 to 176.2 | 0.4415 | 0.1365 | 0.1433 | 0.4588 | yes |
| looking time to agents (ms) | 38 | -3.796 | -458.7 to 418.8 | 0.9868 | 0.8974 | 1.0000 | 0.9913 | yes |
| looking time, proportion of trial | 38 | -0.007462 | -0.04142 to 0.02638 | 0.6760 | 0.7308 | 0.6271 | 0.6779 | yes |
| gaze transitions per trial | 38 | 0.06579 | -0.3816 to 0.4605 | 0.7632 | 0.6025 | 1.0000 | 0.7710 | yes |
| pupil, approach response (GA vs CA) | 28 | 103.4 | -46.49 to 264.5 | 0.2093 | 0.3620 | 1.0000 | 0.2100 | yes |
| pupil, leaving response (CG vs GG) | 28 | 179.7 | 34.53 to 366.3 | 0.0448 | 0.0946 | 0.5716 | 0.0403 | NO |

Two readings of that table matter.

First, the pupil p-value moves from 0.0119 (28 dogs, complete cases) to
0.187 (37 dogs, same file without the exclusion) to 0.459
(38 dogs, `ALL_DATA`), while the four routes agree with one another *within* every sample.
This is not a modelling-route dependence — the four routes agree with one another *within* every
sample — it is a sample-definition dependence, and **the paper owns it explicitly**: the methods
state the rule and its rationale ("to ensure a within-subject comparison, we only retained dogs
with four valid test trials, one per condition (N = 28)"), and the results report the all-data
re-analysis and its null (t = 0.54, p = 0.592). What
this review adds is the magnitude of the dependence at the subject level, which the paper's own
GAMM t-statistics do not convey: the same contrast moves from p = 0.0119 to
p = 0.459, and its interval widens from [52, 283] to
[-108, 176] a.u. Whether the rule was fixed before the outcomes were
examined is still not readable from either the paper or the code — it is the kind of thing a
preregistration answers, and none is cited — but the rule itself is principled, stated, and
reproduces exactly against the data (section 4).

Second, the leaving-response contrast (CG vs GG) is the one place in this review where the routes
**do** disagree at alpha: paired t p = 0.0448 and permutation p = 0.0403
fall below 0.05 while Wilcoxon p = 0.0946 and the sign test p = 0.5716 do
not. On the first pass, when this contrast was wrongly taken to carry C5, that disagreement drove
an `assumption-dependent` verdict. **The paper does not report this contrast**, so the
disagreement is a property of a contrast only this review computed, and it bears on no claim in the
ledger; it is reported because the route table above is reported in full, including the rows that
disagree, and because it is the reason the per-gesture contrasts should not be read individually.
No route was selected after the fact in any contrast here.

**What the review could not test.** Whether the GAMM's smoothing parameters, its `k`, or the
`s(Xgaze_RIGHT, Ygaze_RIGHT)` term change the `trial_type` effect; whether the beta model's
zero/one rescaling affects the looking-time null; whether `log(TRIAL_DWELL_TIME)` is the right
exposure for the transition counts. All three require refitting the authors' model classes, which
this environment cannot do.

---

## 6. Resolution and evidence — stage 5

![F5 — null distribution]({{artifact:art_4c85548c-49be-46bb-953a-f5c7f12296c0}})

*Sign-flip permutation null for C1 with the observed mean difference marked, and the bootstrap
distribution of the estimate with its 95% interval.*

![F7 — time-resolved pupil contrast]({{artifact:art_8c47fa3b-3efb-4848-a35b-17fdfc3c56a3}})

*The same contrast resolved over the 41 bins of the analysis window, in both samples. Left: the
authors' primary sample, bins exceeding |t| > 2.05 shaded, largest cluster
7800–8700 ms of video time, cluster-mass permutation p = 0.0196. Right:
`ALL_DATA`, where no bin reaches the threshold at all.*

**Resolution — every p-value next to the smallest its design could have produced.** No floor binds
anywhere in this review:

| contrast | n pairs | permutation p | smallest achievable p (sign-flip, 2^-n) | sign-test floor | resolution at n_perm = 10000 | floor binds |
|----------|---------|---------------|------------------------------------------|-----------------|---------------------------------|-------------|
| pupil, primary complete cases (dat3) | 28 | 0.0119 | 7.45e-09 | 7.45e-09 | 1e-04 | no |
| pupil, ds file, no exclusion | 37 | 0.1866 | 1.46e-11 | 1.46e-11 | 1e-04 | no |
| pupil, ALL_DATA (abstract's N = 39) | 38 | 0.4588 | 7.28e-12 | 7.28e-12 | 1e-04 | no |
| looking time to agents (ms) | 38 | 0.9913 | 7.28e-12 | 7.28e-12 | 1e-04 | no |
| looking time, proportion of trial | 38 | 0.6779 | 7.28e-12 | 7.28e-12 | 1e-04 | no |
| gaze transitions per trial | 38 | 0.7710 | 7.28e-12 | 7.28e-12 | 1e-04 | no |
| pupil, approach response (GA vs CA) | 28 | 0.2100 | 7.45e-09 | 7.45e-09 | 1e-04 | no |
| pupil, leaving response (CG vs GG) | 28 | 0.0403 | 7.45e-09 | 7.45e-09 | 1e-04 | no |
| C5 interaction | 28 | 0.5237 | 7.45e-09 | 7.45e-09 | 1e-04 | no |
| C1 cluster permutation | 28 | 0.0196 | 3.73e-09 | — | 1e-04 | no |

With 28 within-dog pairs the sign-flip design admits p-values down to 7.45e-09, so
every value reported here is one the design could have produced and none sits against a floor. The
binding limit is the **resolution** of the permutation itself: 10000 permutations resolve no
finer than p = 1e-04, and no number in this report should be read below that.
The place a floor would have mattered is the sign test, whose p = 0.0357 on C1 reflects
that 8 of 28 dogs moved against the claimed direction and
20 with it — a count-based statement that is weaker than the mean-based one and agrees
with it.

**How the permutation was done here.** Signs were flipped **within dog**, which encodes the null
that congruent and incongruent are exchangeable within a dog while preserving between-dog
variation. Shuffling trial labels across the whole dataset would encode a different and weaker
null. How the paper's own 1000 bootstraps were resampled cannot be stated: `boot.glmmTMB` and
`boot.glmm.pred` resample at the model level, and whether the published intervals are those
parametric-bootstrap intervals could not be verified without the paper.

**Intervals.** C1, 28 dogs: 141.6 a.u., 95% BCa bootstrap CI 51.6 to
282.6, n_boot = 10000, bootstrap SE 56.6, bias
0.74. For the two nulls the interval *is* the finding: looking time
-3.8 ms, CI -459 to 419 ms; transitions
0.066, CI -0.38 to 0.46. Both intervals contain zero, and
both also contain differences large enough to matter behaviourally — roughly a
11% change in looking time
and a 21% change in transition
rate. "Did not differ significantly" is a correct description of these data; "did not differ" is
not supported by them, and no equivalence bound was pre-specified in anything this review could
see. This is the single change that would most strengthen the two null claims: state a smallest
effect of interest and test against it, rather than reading non-significance as absence.

---

## 7. Multiplicity — stage 6

![F6 — multiplicity]({{artifact:art_6ddcc869-0b54-4d1d-bb0c-122df51def20}})

*Assessable on this pass. The 14 p-values reported in the paper's main text, against
the uniform expectation, with hits surviving each correction on the right. Bounds reported as
p < .001 are plotted at .001, which is the reading most favourable to the paper.*

**Counting what was run.** Two counts, which differ and should: 14 p-values appear
in the paper's main text, while the authors' deposited code contains 18 inferential tests.
The code is the larger number because it includes diagnostics and model-comparison steps the paper
does not print:

| analysis | inferential tests |
|----------|-------------------|
| pupil GAMM, complete cases (dat3): trial_type term + s(time, by=trial_type) | 2 |
| pupil GAMM, ALL_DATA (dat2): trial_type term + s(time, by=trial_type) | 2 |
| pupil, subject-aggregated paired t-test | 1 |
| pupil, Shapiro-Wilk on the difference scores | 1 |
| pupil, exact Wilcoxon signed-rank | 1 |
| dwell proportion beta GLMM, trial_type x agent_type interaction (drop1) | 1 |
| dwell proportion beta GLMM, main effects (drop1: trial_type, agent_type, z.trial) | 3 |
| dwell time, reacting vs communicating agent, paired t-test | 1 |
| dwell proportion beta GLMM, condition x agent_type vs reduced model (LR test) | 1 |
| transition count Poisson GLMM (drop1: trial_type, z.trial) | 2 |
| transition count, zero-inflation simulation | 1 |
| overdispersion tests (two models) | 2 |
| **total visible in deposited code** | **18** |

At alpha = 0.05 with m = 18, 0.90 significant results are expected by
chance alone; with m = 14 reported, 0.70. **The words
"Bonferroni", "false discovery" and "correct for multiple" do not appear in the paper**, and no
correction is described in its analysis section. Recounting the reported p-values under two
standard corrections:

| correction | controls | n significant at 0.05 | threshold on raw p |
|------------|----------|--------------------------|--------------------|
| none | nothing (per-test error rate) | 7 | 0.0500 |
| Bonferroni | FWER: probability of any false positive | 5 | 0.0036 |
| Benjamini-Hochberg | FDR: expected fraction of false calls <= 0.05 | 6 | 0.0110 |

Two things qualify that table before it is read as a criticism. First, the 14
reported p-values are not one family: they span three outcomes (pupil size, dwell proportion, gaze
transitions) with largely pre-specified hypotheses, and several are control predictors (session
number, trial number) or model-comparison steps rather than tests of the paper's claim. A blanket
correction across all 14 is the most severe reading available, not the correct one.
Second, the claim-bearing results are not marginal: the pupil Wilcoxon at
p = 0.011 and the agent-type effect at p < .001 are the two that carry
conclusions, and the two null results are unaffected by any correction. What the review can say
plainly is that the family over which the paper's alpha is controlled is not stated, and stating
it would let a reader do this arithmetic themselves.

**Correction applied to this review's own p-values.** 9 permutation p-values were
computed here; correcting *them* as a family is a statement about this review, not about the paper:

| correction | controls | n significant at 0.05 | threshold on raw p |
|------------|----------|--------------------------|--------------------|
| none | nothing (per-test error rate) | 3 | 0.0500 |
| Bonferroni | FWER: probability of any false positive | 0 | 0.0056 |
| Benjamini-Hochberg | FDR: expected fraction of false calls <= 0.05 | 0 | 0.0000 |

Expected by chance across these 9:
0.45. Uncorrected, 3 are significant
(C1 pupil, 28 dogs (permutation); C1 pupil, 28 dogs (cluster permutation); C5 leaving contrast (permutation)). Under Benjamini-Hochberg
0 survive (none); under Bonferroni
0 survive (none). The leaving-response
contrast of C5 (p = 0.0403) is
not among the survivors of either correction,
which — independently of the route disagreement — is a second reason not to read it as establishing
context-specificity.

**Distribution of p-values.** The diagnostic now has inputs but still has no resolution, and that
limit is the point rather than a preamble to it. A just-below-threshold diagnostic is a property of
a *collection* of p-values: it describes a set, says nothing about any individual result in it, and
says nothing about why any result looks the way it does. **14 reported p-values is
below the roughly 20 needed for the diagnostic to distinguish anything from chance**
(`interpretable = False`), so while the counts can be printed —
0 in the window (0.045, 0.05) against
1 in (0.04, 0.045), binomial p for an
excess = 1.000, 0 exactly at
threshold — they are not interpretable, and this review does not interpret them. The skill's own
note is carried rather than paraphrased: *"fewer than 20 p-values: this diagnostic is not interpretable at this scale; report the counts without drawing a conclusion."* No inference about analytic flexibility,
selective reporting or intent is drawn here; none could be drawn from 14 values, and
none would be drawn about an individual result from any number of them.

---

## 8. Per-claim supplement

![C1]({{artifact:art_308f4ff7-bf59-4808-a4fd-508e493db80a}})

*C1 — pupil size change, 28 dogs. Each dog's two values, the sign-flip permutation null, and
the four routes side by side. All agree at alpha.*

![C2]({{artifact:art_9156a7d5-b796-4f3f-8d9d-6ae8bdea65e0}})

*C2 — mean dwell time to the agents, 38 dogs. The interval, not the p-value, is the result.*

![C3]({{artifact:art_533eeb93-0df0-4c8d-8fb1-3f4c15204488}})

*C3 — gaze transitions per trial, 38 dogs.
34 of 152 trials contain no transition at all, so
the per-dog mean summarises a count with a floor.*

![C5]({{artifact:art_98d8b83e-6632-4102-af15-d71882cc8940}})

*C5 — the gesture-by-congruency interaction, 28 dogs. Left: each dog's congruency
contrast within each gesture context. Centre and right: the permutation null and the routes for
their difference, -76.3 a.u. (95% CI -285.7 to
177.8, permutation p = 0.524). The paper runs no
such test, and does not need to: its context-relativity argument is the action-balanced design,
not a contrast between these two columns. The null here is consistent with that argument.*

C4, C6 and C7 carry no per-claim panel: C4 is a count, verified in section 4, and C6 and C7 are
`not-assessable` with no contrast to draw.

---

## 9. What was not assessed

Diagnostics in this report describe distribution shape and spread only. They cannot evaluate:
  - independence of observations
  - the sample is representative of the stated population
  - absence of unmeasured confounding
  - measurements were not selected after seeing the outcome
These are design questions. Assess them from the methods section, and say so explicitly when the methods do not contain enough information to assess them.

Additionally, specific to this review:

- **The supplementary material.** Tables S1–S4 and the supplementary figures were not supplied.
  They hold the GAMM and GLMM coefficient tables with their degrees of freedom and the authors'
  1000-sample bootstrap CIs, which is why 6 of the 12 extracted test
  statistics could not be arithmetic-checked and why no reported interval appears anywhere in this
  report. Supplying them would close both gaps.
- **Reviewer error on the first pass, recorded rather than silently corrected.** With only the
  abstract, this review returned `supported-with-caveat` on C1 against an undisclosed sample split
  and `assumption-dependent` on C5 against a difference of significance. The paper discloses the
  split and makes no per-gesture claim; both verdicts were withdrawn on reading it. The arithmetic
  in both cases was correct — the reading of the prose was not, which is the failure mode of
  reviewing an abstract as though it were a paper.
- **The authors' models.** Neither the `mgcv::bam` factor-smooth GAMM, the `glmmTMB` beta model,
  nor the `glmer` Poisson model with offset was refitted — none is available here, and refitting
  them would be redoing the analysis rather than auditing it. The contrasts were re-derived at the
  subject level, and a cluster-based permutation test was added for pupil size because the window
  mean and the GAMM difference curve are different estimands. **A subject-level null on `ALL_DATA`
  does not demonstrate that the authors' GAMM is wrong;** it shows the effect is not visible in the
  dogs-as-units summary of the same data, and the discrepancy is one the authors can address in a
  sentence by stating which sample the published pupil result comes from.
- **Claims with no supporting data:** C6 (not self-centred) and C7 (ontogeny), both
  `not-assessable`. C6 would be changed by a first-person control condition; C7 by any
  developmental, rearing-history or longitudinal measure.
- **Inclusion and exclusion as written.** The eleven-dog exclusion was verified against the data
  (it is exactly the set lacking a condition) and the four excluded sessions (`Nala_2`, `Sahibu_2`,
  `Melody_1`, `Zserbo_1`) were reproduced from the code's own comments, but whether these criteria
  were fixed before the outcomes were examined is a methods question, readable only from the paper,
  and is not assessed.
- **Nothing was subsampled, filtered or approximated** beyond the authors' own documented filters.
  Every number here comes from the full deposited files. Where a dog drops out of a contrast it is
  because one of its two trial types is absent, and the n is stated on every panel and in every
  table.
- **Not computed, by choice:** residual diagnostics for the authors' models, power curves, and
  equivalence tests against a bound the paper did not state.

---

## 10. Reproduction

```bash
python lonardo_2026/review_lonardo_2026.py
```

| | |
|--|--|
| Seed | `random_state = 20260101` — sign-flip permutations, bootstraps and the cluster permutation |
| Permutations / resamples | `n_perm = 10000`, `n_boot = 10000`, cluster permutations = 10000 |
| Cluster-forming threshold | two-sided t at alpha = 0.05: 2.052 (28 dogs), 2.026 (38 dogs) |
| Bootstrap method | BCa percentile for paired differences (`resampling.bootstrap_ci`) |
| Package versions | numpy 2.4.6, scipy 1.17.1, pandas 2.3.3, pyarrow 25.0.0 |
| Data checksums | section 1 (sha256, first 16 hex) |
| Skill version | statistical-review v1 |
| Pipeline check | the authors' own commented values for mean dwell time to the reacting (5105.282 ms) and communicating (3217.737 ms) agent were reproduced to three decimals from the raw IA reports, confirming this review's dwell-time pipeline matches theirs |

---

*Findings in this report are arithmetic statements about recomputed values and about the inputs
named above. No reported value was available to compare them against, and where they differ from
the paper the cause may be an input this review could not see; the inputs used are named so that
any discrepancy can be resolved. Nothing here attributes intent.*
