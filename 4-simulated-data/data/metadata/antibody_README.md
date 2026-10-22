# Project Fulcrum — Which Antibodies Switch the Receptor On?

## The story

A biotech team is developing **agonist antibodies**: antibodies that don't block a cell-surface receptor
but *activate* it, by grabbing receptors and pulling them into the right arrangement to signal.
Their screening campaign produced **500 candidate antibodies**. Each one binds the same receptor at one
of **5 epitope sites** (numbered 1–5), which sit at increasing distances from the cell membrane. Each
antibody was also characterised by surface plasmon resonance, giving its association rate (`kon`) and
dissociation rate (`koff`). Every antibody was then tested at the same dose in a cell-based
signalling assay, and its activity was recorded relative to a reference agonist.

The team has to pick which antibodies to take forward, and, more usefully, to learn **what design
rules drive agonism**, so the next round can be engineered rather than screened.

> **Your task: build linear models of `agonism_pct`. What drives agonism, how, and how well can
> you predict it? Which antibodies would you take forward, and what would you tell the team to
> look for in the next round?**

Expect a plain linear fit to disappoint. Part of the task is working out why, and what to do about it.

## The data

`antibody_agonism.csv` — 500 rows, one row per antibody.

| Column | Type | Description |
|---|---|---|
| `ab_id` | string | Antibody identifier (prefix is its site, e.g. `S3_17`). |
| `site` | integer | Epitope site, 1–5. |
| `epitope_dist_nm` | numeric | Distance of that epitope from the cell membrane, nm. Fixed per site. |
| `log10_kon` | numeric | log10 of association rate, M⁻¹ s⁻¹. |
| `log10_koff` | numeric | log10 of dissociation rate, s⁻¹. |
| `log10_KD` | numeric | log10 of equilibrium dissociation constant, M (`log10_koff − log10_kon`). |
| `receptor_density` | numeric | Relative receptor expression on the assay cells (1 ≈ typical). |
| `plate` | integer | Assay plate, 1–6. Antibodies were spread across plates. |
| `agonism_pct` | numeric | **Outcome.** Signalling activity as a percentage of the reference agonist, floored at 0. |

## Notes from the screening lab

- The panel is complete; no missing values.
- Every antibody was tested **once, at a single fixed dose**.
- Plates were run on different days. Plate-to-plate drift was not calibrated out.
- `kon` and `koff` were measured independently for each antibody, so they are not tied to each other
  (`log10_KD` is just their difference — don't put all three in a model without thinking about why).
- Per-site distances are the same for all antibodies at a site, so `site` and `epitope_dist_nm` carry the same information.

## Questions worth answering

1. How does agonism vary across the 5 sites? Is it monotone?
2. How does it vary with `kon` and `koff`? Is a straight line in the raw or log scale reasonable?
3. Do the kinetic effects depend on the site (interactions)?
4. How much does the plate matter, and does it change your conclusions?
5. What is your best predictive model, and what is its held-out performance?
6. What would you recommend for the next round of designs, and where is your model unable to say?

## Deliverable

A short write-up with one figure, the model you recommend, its coefficients with uncertainty, and an
honest account of what it can and can't tell the team. Hold some data out before you start tuning.
