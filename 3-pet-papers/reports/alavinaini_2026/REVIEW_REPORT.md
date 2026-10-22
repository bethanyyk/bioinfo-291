# Statistical review (v2 -- data-based)

**Alavi Naini L., Linde N., Malik S., Syed Y. & Cheek A.O. (2026)** "Diel Activity
Patterns of Bobcats and Domestic Cats in the Houston Metropolitan Area."
*microPublication Biology*. DOI: [10.17912/micropub.biology.001911](https://doi.org/10.17912/micropub.biology.001911)

Reviewed with the `statistical-review` skill, six stages in order.

> **This version supersedes v1.** v1 was a text-only review: the Dryad deposit the paper
> names could not be retrieved from the review sandbox. The deposit has since been
> supplied and verified against the Dryad API manifest, and every claim below is now
> re-derived from the detection-level data and the authors' own analysis scripts.

---

## 1. Header -- what was reviewed

| input | condition |
|---|---|
| paper | full text from the PMC XML, 20,192 characters; text layer usable |
| raw data | **obtained.** `Bobcat_Data.csv` (81 rows), `Domestic_Cat_Data.csv` (795 rows), both 22 columns, all rows Status "confirmed correct" |
| analysis code | `Watsons U2 and Polar Plot.m` (the published U2 values) and `Dhat R Script MicroPub.R` (the published Dhat values), both deposited |
| data availability | the paper **does** carry a statement: "Animal detection data, MatLab script and R script are available at datadryad.org, DOI: 10.5061/dryad.s4mw6m9jg" |
| provenance | v1 recorded the deposit as existing but unreachable (HTTP 401 on the file endpoint, then loss of network egress). The user downloaded it manually; both files match the Dryad API manifest byte-for-byte. **The v1 data gap was a retrieval failure in the review environment, not a reporting failure by the authors** -- the paper's data availability statement is accurate and the deposit is complete enough to reproduce the headline numbers. |

**Reproduction is exact.** Every headline statistic re-derives from the deposit:

| statistic | paper | this review |
|---|---|---|
| Watson's U2, single-species sites | 0.0753 | **0.075326** |
| Watson's U2, co-occurrence sites | 0.189 | **0.189346** |
| Dhat, single-species sites | 0.852 | **0.8523** |
| Dhat, co-occurrence sites | 0.647 (deposit plot: 0.646) | **0.6467** |
| nocturnal, bobcat at bobcat-only | 67% | **33/49 = 67.3%** |
| nocturnal, cat at cat-only | 85% | **594/696 = 85.3%** |
| nocturnal, cat at co-occurrence | 51% | **51/99 = 51.5%** |
| detections / sites | 81, 795 / 7, 19, 7 | 81, 795 / 7, 19, 7 |

**Coverage.** 11 claims, all with a verdict. Three were `not-assessable` in v1 and are now
assessable (C3, C8, C9); two remain `not-assessable` (C5, C10) because the deposit does not
contain the trap-days, camera coordinates or land-cover measurements they need.

---

## 2. Summary of findings

| verdict | n |
|---|---|
| `supported-with-caveat` | 5 |
| `assumption-dependent` | 3 |
| `not-assessable` | 2 |
| `supported` | 1 |

**Finding 1 -- the two analyses run on two different clocks (stage 3).** The deposited
MATLAB reads column 9, `Local Time Value (ABS)`, a day fraction built on a **fixed UTC-6
offset that does not follow daylight saving**: 703
of 876 detections are an hour behind true local time
there. The deposited R recomputes local time DST-aware from the UTC timestamp. Worse, the
Methods say timestamps were "converted to solar time using the sunTime function" and that
"Watson's U2 test was used to compare diel activity" -- but the solar conversion exists
only in the R script, which computes Dhat. **The published U2 values were computed on raw,
non-DST clock time.** Recomputing U2 on the basis the Methods describe gives
0.0922 (p = 0.324) and 0.2254
(p = 0.0231) -- different numbers from the published 0.0753 and 0.189.

**Finding 2 -- C2's significance does not survive the design (stages 2 and 5).** U2 = 0.189
reproduces exactly, but its p depends on what is assumed exchangeable: 0.0465
shuffling every detection freely (the published analysis), **0.0712
shuffling only within site**, and **0.0974** after a standard 30-minute
independent-event filter. v1 predicted this from the reported numbers alone -- it calculated
that a design effect of only 1.011 would cross alpha -- and the data confirm it. Two of the
seven shared sites supply 91 of the 99 co-occurrence cat detections.

**Finding 3 -- the central claim's interval includes zero, but a better test does not
(stage 5).** The paper's temporal-partitioning claim is the gap between Dhat 0.852 and
0.647. The deposited R script computes bootstrap CIs and the paper reports none; they are
[0.729, 0.897] and
[0.533, 0.805], and they overlap. The difference
is 0.206 with a 95% interval of
[-0.004, 0.308] -- **including zero**. The biology
is nonetheless on firm ground by a route the paper never takes: cat activity differs
between the two site classes (U2 = 1.515,
p <= 5e-05) while bobcat activity does not
(U2 = 0.079, p = 0.429).
**The paper should make that claim instead.**

---

## 3. Claim ledger -- F4

![F4]({{artifact:art_855682ae-cd1c-4a5d-8ec9-431ba502166c}})

| claim | verdict | v1 | effect (paper's units) | recomputed | p rep. | p recomp. |
|---|---|---|---|---|---|---|
| **C1** diel activity similar at single-species sites | `assumption-dependent` | confirmed | 0.0753 Watson's U2 | 0.07533 | 0.2000 | 0.4511 |
| **C2** cat vs bobcat activity differs at co-occurrence sites | `assumption-dependent` | downgraded | 0.189 Watson's U2 | 0.1893 | 0.0500 | 0.0465 |
| **C3** activity overlap much lower where species co-occur | `assumption-dependent` | now assessable | 0.205 difference in overlap coefficient Dhat1 | 0.206 [-0.00361, 0.308] | -- | 0.0305 |
| **C4** land cover differs between bobcat-only and cat-only sites | `supported-with-caveat` | unchanged | -32 percentage points of developed cover | -32 [-49.9, -14.1] | -- | 0.0031 |
| **C5** each species detected more often either side of 65% developed | `not-assessable` | confirmed not assessable | -- detection frequency | -- | -- | -- |
| **C6** bobcats primarily nocturnal wherever present | `supported-with-caveat` | confirmed, missing half now supplied | 17 pp above a 50% night share | 17.3 [3.38, 28.8] | -- | 0.0213 |
| **C7** cats shift to more daytime activity where bobcats occur | `supported-with-caveat` | confirmed and strengthened | 34 pp of nocturnal detections | 33.8 [-3.57, 41.1] | -- | 3.80e-13 |
| **C8** cat activity shifts into 03:00-10:00 and 16:00-21:00 | `supported-with-caveat` | now assessable | -- pp of detections in the named windows | 18.81 | -- | 5.40e-04 |
| **C9** few cats were collared, so they may survive by hunting | `supported-with-caveat` | partly assessable | -- % of scored cats collared | 7.576 | -- | -- |
| **C10** spatial and temporal partitioning facilitate coexistence | `not-assessable` | confirmed not assessable | -- n/a | -- | -- | -- |
| **C11** 81 bobcat and 795 cat detections across 33 sites | `supported` | confirmed against the rows | -- counts | -- | -- | -- |

Full `reasons`, both p-values, the v1 verdict and the design floor for every claim are in
`claim_ledger.csv`; every reported statistic is checked against the data in `statcheck.csv`.

### What the data confirmed, and what they overturned

**Confirmed from v1:**

- C1's printed "p < 0.20" points the wrong way. v1 derived p = 0.447 asymptotically; the
  permutation p from the authors' own data is 0.451.
- C2 is fragile to clustering. v1 computed that a design effect of 1.011 would cross alpha;
  the within-site permutation gives 0.0712.
- v1's GRIM reconstruction of the nocturnal counts: 33/49 recovered uniquely and correctly;
  594 falls inside the reconstructed range 589-595. Its 50/99 was one out -- the value is 51/99.
- C4 and C5 are unchanged: the deposit contains no land cover, no trap-days and no camera
  coordinates.

**Overturned or newly resolved:**

- **C2 downgraded** from `supported-with-caveat` to `assumption-dependent`: with the data,
  three of five defensible analyses put it above alpha.
- **C3 now assessable** and downgraded in substance: the bootstrap interval on the Dhat
  difference includes zero.
- **C8 and C9 now assessable.** C8's first half holds and its second half does not. C9's
  collar annotation exists in the deposit: 5 collared against
  61 not, on 99 annotated of 795.
- **C6's missing half is now supplied and does not support the claim.** Bobcats at
  co-occurrence sites are 21/32 =
  65.6% nocturnal, binomial p = 0.110,
  with a Wilson interval that includes 50%.

---

## 4. Design and sampling -- F1, F2

![F1]({{artifact:art_cf353464-87eb-4b7b-82d3-9b9ba9133903}})

**Stage 1.** The population is free-ranging bobcats and unowned domestic cats across the
22,891 km2 Houston metropolitan area; the sample is 33 camera sites in parks, cemeteries,
golf courses and natural areas along two transects. Whether that convenience sample
represents the metropolitan area is a design question the methods do not address and the
data cannot settle.

**Stage 2 -- pseudoreplication, now measured rather than bounded.** v1 could only show a
sensitivity curve over an unknown ICC. The site labels in the deposit make it concrete:
cat detections are distributed across 26 sites with mean 30.6
and s.d. 46.8, and **one site contributes
234 of the 795 cat detections (29%)**.
Bobcat detections are sparser and flatter (mean 5.8,
max 14 across 14 sites). Within the co-occurrence
comparison, two of seven sites supply 91 of 99 cat detections and 22 of 32 bobcat detections.

No independence filter is applied anywhere in the paper or the deposited code: consecutive
photographs of the same animal enter as separate detections. The largest single burst is 14
cat photographs spanning 6.8 minutes at one site. A 30-minute quiet interval reduces 876
detections to 777 events -- an 11% reduction overall, but
enough to move C2 (section 6).

The totals all check out against the rows (C11, `supported`): 81 and 795 detections, site
classes of 7/19/7 derived by the deposited script's own rule, arms of 49/696/32/99.

![F2]({{artifact:art_e251cef2-1ba6-49c0-8cce-4e97e75039d4}})

---

## 5. Method and assumptions -- F3

![F3]({{artifact:art_2f00a41c-f7b2-4c4a-b17a-73dca6f13115}})

**Stage 3 -- test selection.** Watson's U2 is an appropriate omnibus two-sample test for
circular data and the deposited implementation is the standard two-sample formula, correctly
rotation-invariant. Three observations:

1. **It is an omnibus test**, so it is least efficient against an alternative concentrated in
   one contrast. At single-species sites the night/day contrast is exactly such an
   alternative, and the two routes disagree (C1).
2. **The time base is wrong relative to the Methods** (Finding 1). The statistic is invariant
   to the minutes-versus-radians scaling -- that part is harmless, since the transform is
   monotone and U2 depends only on sort order -- but it is not invariant to the DST offset
   or to the solar-time transformation.
3. **Ties.** The implementation uses MATLAB's stable sort, so cross-group ties resolve in
   favour of the first-listed group. There is exactly one cross-group tie in the
   single-species comparison and none in the co-occurrence comparison, so tie handling does
   not affect either published value.

**Stage 4 -- assumptions.** Watson's U2 is distribution-free: no normality assumption is at
stake, no shape diagnostic applies, and nothing here is scored pass or fail. The assumption
that carries the inference is **exchangeability of detections**, and F3(a) prices it
directly -- the same observed statistic maps to 0.0465,
0.0712 or 0.0974 depending on which
exchangeability the null encodes.

The three routes in the skill's usual sense (parametric, rank, permutation) do not transfer
to a circular omnibus statistic with no parametric analogue in use here. The analogous set
is reported instead and in full: the asymptotic table the authors used, the free-shuffle
permutation, the within-site permutation, the event-filtered version, and both time bases.
They are listed together in panel C2 and none was selected after the fact.

On error bars: the land-cover "+/- " values are labelled 1 s.d. for two of three site
classes. Read as s.e.m. they would imply site-level s.d. of 50 and 61 percentage points on
a 0-100 bounded scale, which the bound rules out, so the s.d. label is the consistent
reading. The site values themselves are not in the deposit, so this remains an argument
from consistency, not a verification.

---

## 6. Resolution and evidence -- F5

![F5]({{artifact:art_92f88b1f-e483-46f1-a96f-50cd58acf90c}})

The permutation null computed from the authors' own detection times tracks the asymptotic
table closely (empirical 95th percentile 0.1845 against the tabulated
0.1869), which confirms both that reading p off that table was the right procedure
and that v1's asymptotic re-derivations were correct.

| | C1 single-species | C2 co-occurrence |
|---|---|---|
| U2 (reproduced) | 0.075326 | 0.189346 |
| printed | p < 0.20 | p < 0.05 |
| permutation p, published basis | **0.4511** | **0.0465** |
| within-site permutation | not defined (see below) | **0.0712** |
| 30-min independent events | 0.3719 | **0.0974** |
| permutation p, solar basis | 0.3237 | 0.0231 |
| solar basis, within-site | not defined | 0.0378 |
| smallest p this design could give | ~1e-77 (combinatorial) | ~1e-30 (combinatorial) |
| smallest p this review could give | 5.0e-05 (B = 20,000) | 5.0e-05 |
| smallest p the paper's table could print | p < 0.001 | p < 0.001 |

**No combinatorial floor binds** at these sample sizes. The binding limits are the Monte
Carlo resolution of this review and, for the paper, the granularity of the printed bracket.

**A within-site null does not exist for C1, and that is a structural finding.** At
single-species sites each site has, by construction, only one species, so species and site
are perfectly confounded: no permutation that respects site can test the species contrast.
Any comparison between single-species classes is a comparison between different places. The
same applies to C7. This is a limit of the design, not of the analysis.

---

## 7. Multiplicity -- F6

![F6]({{artifact:art_c82bf004-b096-4971-b7e4-447ef3fc6c20}})

The paper reports two inferential tests. The single uncorrected hit (C2,
p = 0.0465) survives neither Bonferroni nor Benjamini-Hochberg (both
0.0931). Whether two
pre-specified comparisons on disjoint site sets form a family that should be corrected at
all is a defensible argument for the authors; both numbers are given.

Untested comparisons described in the text remain untested: six land-cover contrasts, three
nocturnal percentages, the Dhat difference, and the hourly-window claim. This review supplies
tests for several of them, but those are the reviewer's and are not counted in the paper's
family.

**On the p-value distribution:** with m = 2 it cannot be read. The just-below-threshold
diagnostic is a property of a collection of p-values; it says nothing about any individual
result, and nothing about why any result looks the way it does. No inference about the
authors' conduct is drawn here or anywhere in this report.

---

## 8. Per-claim supplement

### C1 -- "the diel activity patterns were similar where only one species occurred" -- `assumption-dependent`
![C1]({{artifact:art_89ca2efa-585a-443b-bf27-5de9bcb3e10b}})

### C2 -- "domestic cat diel activity differed significantly from bobcat activity" -- `assumption-dependent`
![C2]({{artifact:art_dbce4892-83cc-4a34-a8ec-7bcbbe4f2afd}})

### C3 -- "activity patterns overlapped much less where they co-occurred" -- `assumption-dependent`
![C3]({{artifact:art_fc8b8332-3260-4dd4-9e64-6df1e5d38fa9}})

### C6 / C7 -- nocturnal fractions -- `supported-with-caveat`
![C7]({{artifact:art_0a6b93df-169b-4ae5-9fa3-1f6d4cb9eca3}})

### C8 -- the named activity windows -- `supported-with-caveat`
![C8]({{artifact:art_b3cf6e6d-430f-4b8b-a313-e955345c1a60}})

---

## 9. What was not assessed

Diagnostics in this report describe distribution shape and spread only. They cannot evaluate:
  - independence of observations
  - the sample is representative of the stated population
  - absence of unmeasured confounding
  - measurements were not selected after seeing the outcome
These are design questions. Assess them from the methods section, and say so explicitly when the methods do not contain enough information to assess them.

**Specific to this review:**

- **Independence is now measured but not resolved.** The site labels quantify clustering and
  the timestamps quantify photographic bursts, but the deposit contains no individual
  identities. 876 detections come from an unknown number of animals; the 1 km camera spacing
  separates cameras, not successive visits by one animal. Every p-value here treats
  detections, or 30-minute events, as the unit; none treats individuals as the unit, because
  nothing in the deposit can.
- **Confounding of species with site cannot be removed** for C1, C7 or C8: single-species and
  co-occurrence sites are different places, so "bobcat present" is confounded with everything
  else that differs between them. Only C2 compares species within the same sites.
- **Claims the deposit still cannot reach:** C5 (no per-site trap-days; no camera
  coordinates -- every row carries the same Latitude/Longitude, the downtown reference point,
  so Fig 1B is not reproducible either) and C10 (coexistence is not measured and no site
  changes predator composition within the study).
- **C4 still rests on published summaries**, not site-level values; its Welch tests apply
  normal theory to percentages bounded at 0-100 at n = 7 and support a direction, not a
  precise p.
- **The night boundary is not defined in the Methods** and no definition tested here
  reproduces all three printed percentages exactly. Sunrise/sunset on the DST-aware solar
  basis reproduces 67% and 85% and gives 51.5% where the paper prints 51%. The nocturnal
  figures in this review use that definition.
- **`overlap`'s compiled kernel could not be loaded in this environment**, so `Dhat1`, its
  bandwidth rule and its smoothed bootstrap resampler were transcribed from the package's own
  R source and verified by reproducing the deposit's published values (0.852 and 0.646/0.647).
  `suntools::sunriset` was likewise reimplemented (NOAA) and validated against published
  Houston sunrise and sunset times to within about two minutes.
- **Nothing was subsampled.** All 876 detections were used; the 30-minute event filter is
  reported as a sensitivity alongside the full-data result, never in place of it.

---

## 10. Reproduction

| item | value |
|---|---|
| script | `review_alavinaini_2026.py` -- regenerates every figure, CSV and this report |
| data | Dryad doi:10.5061/dryad.s4mw6m9jg, manifest-verified |
| seed | `SEED = 20261015`; per-analysis seeds fixed in the source |
| permutations | B = 20,000 per test (resolution 5.0e-05) |
| bootstrap | 2,000 smoothed replicates per Dhat; 600 site-level cluster replicates |
| alpha | 0.05 |
| Dhat estimator | Dhat1, selected by overlap's own rule (min n < 50) -- the paper does not name it |
| python / numpy / pandas / scipy | 3.11.16 / 2.4.6 / 2.3.3 / 1.17.1 |

Figures: F1_design_map.png, F2_raw_behind_claim.png, F3_assumption_shape.png, F4_claim_ledger.png, F5_null_distribution.png, F6_multiplicity.png, panel_C1.png, panel_C2.png, panel_C3.png, panel_C6_C7_nocturnal.png, panel_C8.png
