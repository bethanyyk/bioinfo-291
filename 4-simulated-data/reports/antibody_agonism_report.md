# Project Fulcrum — what drives agonism?

Analysis of the 500-antibody screen received from my partner
(`data/raw/part2-partner-antibody/antibody_agonism.csv`, SHA-256 `7858ccb2…`).
Reproduce with `python workflows/analyze_antibody_agonism.py`.

![Design rules for antibody agonism]({{artifact:d4be79c9-8f03-44ca-8819-936eab2dabb1}})

**Figure.** Each point is one antibody (n = 500). (a) Agonism against epitope distance;
black bars are medians. (b, c) Grey points are all antibodies; coloured curves are the
fitted model at each site, holding the other kinetic parameter at its median, receptor
density at 1.0 and the plate effect at the panel average. (d) Fitted plate multipliers
with 95% CIs. (e) Held-out predictions from the recommended model fitted on the 375
training antibodies only. "Agonism" is percent of the reference agonist.

## Headline

A plain linear model fails here, and it fails for three separate reasons that all have the
same fix. Agonism is **multiplicative, not additive**: within each site the SD tracks the mean
almost 1:1, so effects compound rather than add. Once that is handled with a log link, the
site and kinetic effects become simple and additive on the log scale, and the apparent
"interactions" that a linear model needs largely disappear.

**Recommended model** — Tweedie GLM, log link, variance power p = 1.3 selected by profile
log-likelihood on the training set:

```
agonism_pct ~ C(site) + log10_kon + log10_koff + log10_koff² + receptor_density + C(plate)
```

The Tweedie family is doing real work: 113 of 500 antibodies read **exactly 0**, which a
log-normal or Gamma model cannot represent and which an ordinary linear model handles by
predicting negative activity.

## Coefficients (fitted on all 500; 95% CIs)

| Effect | Estimate | Reading |
|---|---|---|
| Site 2 / 6 nm | ×5.3 [3.9, 7.3] | vs site 1 (3 nm) |
| Site 3 / 9 nm | **×9.7 [7.2, 13.1]** | best site |
| Site 4 / 12 nm | ×7.6 [5.6, 10.3] | |
| Site 5 / 15 nm | ×2.0 [1.4, 2.8] | |
| log₁₀ k_on | ×1.60 per log₁₀ [1.47, 1.74] | monotone, modest |
| k_off optimum | log₁₀ k_off = **−3.33 [−3.39, −3.27]** | ≈ 4.7×10⁻⁴ s⁻¹ |
| k_off curvature | −0.70 [−0.78, −0.62] | strong inverted U |
| Receptor density | ×1.22 per +1.0 unit, p = 0.17 | not resolved |
| Plate 2–5 vs 1 | ×1.38 to ×1.79 | plates 1 and 6 read low |

## Answers to the six questions

1. **Across sites — not monotone.** Agonism peaks at the middle epitope (9 nm, ×9.7 over the
   membrane-proximal site) and falls away at both ends. Any model treating `site` or
   `epitope_dist_nm` as a linear term will conclude distance barely matters.
2. **Kinetics — curved, and only in the log scale.** k_on acts monotonically but weakly
   (×1.6 per decade). k_off has a sharp interior optimum at 4.7×10⁻⁴ s⁻¹: both faster and
   slower dissociation lose activity. A straight line in either scale misses this entirely,
   which is why the naive `agonism ~ log10_KD` fit explains almost nothing
   (held-out R² = 0.13). Note `log10_KD` is exactly `log10_koff − log10_kon`, so only two
   of the three can enter a model; k_on and k_off enter separately here because their
   effects differ in both shape and size, which a single KD term cannot express.
3. **Interactions — apparent, but explained by the floor.** A linear model needs large
   site×kinetics interactions; the log link absorbs them. Adding them back to the GLM
   still improves fit (dispersion-scaled F = 3.70, p = 3×10⁻⁵), and the gain comes entirely
   from sites 1 and 5 — the two where the zero floor binds hardest (41% and 29% zeros).
   **Restricted to sites 2–4, where the floor barely binds, the k_on slopes are statistically
   identical (F = 0.03, p = 0.97).** I therefore read the interaction as a censoring artifact
   rather than a site-specific design rule, and do not include it.
4. **Plate matters and is not ignorable.** Plates 1 and 6 read ~1.4–1.8× low relative to
   plates 2–5 (dispersion-scaled F = 7.21, p = 2×10⁻⁶). It does not change any conclusion
   about the design rules, but it does reshuffle the candidate ranking, so the shortlist
   below is plate-corrected. With only 6 plates I fitted them as fixed effects; a random
   intercept (GLMM) would shrink these slightly and is the better choice if more plates arrive.
5. **Held-out performance** (125 antibodies, never used for model choice): Spearman 0.88,
   R² 0.85, MAE 5.1 percentage points, and **18 of the true top 20 antibodies recovered in the
   model's top 20**. Calibration is good (observed/predicted 0.75–1.14 across predicted
   quintiles). Note the two loss measures disagree about plate: including it improves held-out
   Tweedie deviance (357→331) but worsens RMSE (6.8→7.5), because RMSE is dominated by the
   handful of high-activity antibodies. Deviance is the matched loss and I weighted it accordingly.
6. **Recommendation for the next round.** Engineer at 9–12 nm, tune k_off to ~5×10⁻⁴ s⁻¹
   rather than maximising affinity, and treat k_on as a secondary knob. Expected gain from
   moving an antibody to the optimum is large: the site and k_off terms together span more
   than two orders of magnitude in predicted activity, whereas the full k_on range contributes
   about ×5.

## What this model cannot tell the team

- **Single dose, single replicate.** Every antibody was tested once at one dose, so none of
  this separates potency from efficacy. An antibody that looks inactive may simply be
  right-shifted. Dose-response on a subset would resolve it.
- **The zero floor censors 23% of the panel.** Below the assay floor the model cannot
  distinguish "weak" from "inactive", and this is precisely where the site×kinetics question
  lives. The site-1 and site-5 slopes should be treated as unresolved, not as zero.
- **Receptor density is not resolved** (×1.22 per unit, p = 0.17) over the narrow range
  sampled (0.49–1.97). This is a power limitation, not evidence of no effect.
- **Site and distance are perfectly confounded** — 5 sites, 5 fixed distances. Nothing here
  says the driver is geometric distance rather than some other property of those 5 epitopes.
  Testing a new epitope at an already-sampled distance would break the confound.
- **Plates were not randomised against anything measured**, so the plate effect is a
  correction, not a mechanism. A plate-1 antibody and a plate-3 antibody are not directly
  comparable without it.
- **Extrapolation.** The k_off optimum sits inside the sampled range, but no antibody has
  k_off below 10⁻⁵ s⁻¹; the downturn beyond the sampled range is an assumption of the
  quadratic, not an observation.

## Shortlist

`results/top20_candidates.csv` ranks the panel by plate-corrected observed agonism.
The top 20 are 15 site-3 and 5 site-4 antibodies with k_off between 10⁻³·⁸ and 10⁻²·⁵ —
i.e. the rules above, read off the data directly.
