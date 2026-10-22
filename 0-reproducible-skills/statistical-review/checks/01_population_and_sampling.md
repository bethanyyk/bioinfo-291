# Stage 1 — Population and Sampling

> Uncertainty comes from sampling. Everything you calculate from a sample is a
> random variable with its own distribution. Before testing anything, define the
> population you are making a claim about.

This stage has **no code**. Nothing in `utilities/` can check it. It is read from
the methods section, and where the methods are silent, the report says they are
silent. Run it first anyway: a sampling problem changes how every number later in
the review should be read.

## Checklist

- [ ] **The paper names the population.** Not "mice" but "vivarium-grown mice of
      this genetic background, this age range, this facility". A claim with no
      stated population cannot be evaluated, only described.
- [ ] **The sample could plausibly have come from that population.** Where were
      subjects recruited, which animals were used, which cells passed QC, which
      patients consented? Each filter narrows the population the claim actually
      covers.
- [ ] **Exclusions are stated and counted.** How many units were removed, at
      which stage, and on what criterion. An exclusion rule applied after seeing
      the outcome changes the population silently.
- [ ] **The comparison groups differ in the intended variable and not obviously
      in others.** Batch, cage, sequencing run, operator, collection date,
      scanner, clinical site. Confounders you can name are checkable; the
      threat is the ones nobody recorded.
- [ ] **Randomization and blinding are described,** or their absence is noted.
- [ ] **The conclusion does not exceed the population.** A claim stated about
      "patients" from a single-centre cohort of one ancestry group is a claim
      about that cohort.

## What to compute

Nothing inferential. Do cross-tabulate the grouping variable against every other
variable in the supplied data — batch, sex, plate, date, site — and report any
that are strongly associated with the group label. The
`data-visualization/abundance_and_proportion/` module covers how to show this;
`sample_composition()` from its utilities prints the cross-tabulation.

```python
import pandas as pd
pd.crosstab(data["arm"], data["batch"], margins=True)
```

An association between the treatment label and a nuisance variable is a
confounding finding, reportable without any test.

## Findings template

> **Population.** The paper claims *[X]* about *[population P as the paper
> defines it, or "a population the methods do not define"]*. The sample was
> *[S]*.
>
> **Representativeness.** *[Specific reason this does or does not plausibly
> represent P.]*
>
> **Confounding.** *[Variable V is associated with the group label: cross-tab
> here. / No recorded variable in the supplied data is associated with the group
> label, which does not rule out unrecorded ones.]*

## Verdict implications

| Observation | Verdict effect |
|-------------|----------------|
| Population not defined | note in report; does not by itself change a verdict |
| Recorded confounder associated with the group label | `supported-with-caveat` at best; name the variable |
| Exclusions applied after seeing the outcome | `supported-with-caveat`; the stated population no longer matches the sample |
| No data corresponding to the claim | `not-assessable`; state what is missing |

## Phrasing

Write what the sample was and what the claim covers. Do not write that the
authors "cannot conclude" something — write that the data support the claim for
the sampled population and that extending it to *[wider population]* needs
*[what]*.
