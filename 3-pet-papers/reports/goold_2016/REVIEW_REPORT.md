# Statistical review - Goold, Vas, Olsen & Newberry (2016)

*Using network analysis to study behavioural phenotypes: an example using domestic dogs.* R. Soc. Open Sci. 3:160268. doi:10.1098/rsos.160268

## 1. Header - what was reviewed, and what could not be

| input | content | sha256 (first 16) |
|---|---|---|
| rawData.csv | 171 rows x 47 cols, 43 descriptors | `4713c991b05722f0` |
| finalMIdata_averaged_rounded.csv | 170 rows x 31 cols, the imputed dataset the networks were built from | `d4d76f3822830936` |
| Rscript.R | 1003 lines, the authors' own analysis script | `9a0af44356647278` |
| paper_text | full text (PMC), 63195 characters | `b78d068af9c3f100` |

The text layer is usable (`usable_text_layer = True`). Automatic extraction found 24 named-test mentions, 7 sample sizes and 7 reported intervals, and **zero p-values and zero test statistics**. `statcheck.csv` is therefore empty by construction, not by omission: this paper reports no null-hypothesis test anywhere in its results, so there is no reported-vs-recomputed p-value to check and no p-value distribution to describe. With zero reported p-values the just-below-threshold diagnostic is far below the 20-value floor at which it means anything, and it is not run.

Claims were therefore read by hand from the abstract, results and methods: **21 claims** are in the ledger, of which **2** are `not-assessable` from the supplied inputs. Having the authors' own `Rscript.R` means the pipeline could be re-derived step by step rather than guessed: the preprocessing, the polychoric correlations, the EBICglasso fit, the centrality measures and the resampling were all reimplemented independently in Python and checked against the paper's printed numbers.

What this review could **not** reach: the 15 individual imputed datasets (only their rounded average was supplied), the authors' bootstrap output tables (`bootTable*.csv`, regenerated here at a smaller replicate count), electronic supplementary tables S1-S5, and the per-descriptor likelihood-ratio tests used to screen handler non-independence (neither the paper nor the script records their output).

## 2. Summary of findings

| verdict | claims |
|---|---|
| `supported` | 5 |
| `supported-with-caveat` | 12 |
| `assumption-dependent` | 0 |
| `underpowered-resolution` | 0 |
| `not-reproducible` | 1 |
| `contradicted` | 1 |
| `not-assessable` | 2 |

Three findings change how the paper reads.

1. **The detection network is one rounding of the tuning parameter away from no network at all.** Its 70 edges beat the empty graph by 19.6 EBIC units, which at 11.78 units per edge is 1.7 edges' worth of penalty; the patrol network's margin is 471.4 units (37.6 edges). The paper says as much in one clause - above gamma = 0.65 the detection graph empties - and the sweep here confirms it exactly (C2, C15, stage 3).
2. **'The most central descriptors denoted desirable characteristics in both networks' holds in one network, not two.** In patrol the 14 desirable nodes average 9.77 against 2.08 for the 6 undesirable ones and every route agrees; in detection the same contrast is 7.86 against 5.02 with no route below 0.05 (C3, C4, stage 4).
3. **The patrol-versus-detection difference is described, never tested.** A label-permutation test of network-structure invariance run here gives max |edge difference| = 0.439 against a null 95th percentile of 0.485, p = 0.257 over 100 permutations (resolution 0.0099); 0 of 190 edges differ at FWER 0.05 (C14, stage 5-6).

## 3. Claim ledger - F4

![F4]({{artifact:art_88fe6778-fb9d-498a-acad-b1476450160b}})

`p_reported` is empty for every row because the paper reports no p-values; `p_recomputed` is filled wherever this review could construct a permutation test of the claim.

| claim_id | claim_short | verdict | effect_reported | effect_recomputed | p_reported | p_recomputed | stage |
|---|---|---|---|---|---|---|---|
| C1 | patrol network has 55 edges (28.95%) | supported | 55 | 55 |  |  | 3 |
| C2 | detection network has 70 edges (36.84%) | supported-with-caveat | 70 | 70 |  |  | 3 |
| C3 | central descriptors are desirable (patrol) | supported-with-caveat |  | 7.695 |  | 0.03129 | 3-4 |
| C4 | central descriptors are desirable (detection) | contradicted |  | 2.835 |  | 0.3564 | 3-4 |
| C5 | PLA most central; CUR top strength in patrol | supported-with-caveat |  |  |  |  | 5 |
| C6 | named descriptors differ in centrality between networks | supported-with-caveat |  |  |  |  | 3 |
| C7 | Cliff's delta quantifies centrality differences | supported-with-caveat |  |  |  |  | 5 |
| C8 | bootstrap analyses revealed network stability | supported-with-caveat |  |  |  |  | 5 |
| C9 | s.d. of edge count: 12.82 patrol, 29.75 detection | supported-with-caveat | 12.82 | 21.76 |  |  | 5 |
| C9b | s.d. of edge count, detection network | supported-with-caveat | 29.75 | 30.52 |  |  | 5 |
| C10 | subject-wise: patrol fails below ~60 dogs, detection ~40 | supported-with-caveat |  |  |  |  | 5 |
| C11 | patrol network more stable than detection | supported |  |  |  |  | 5 |
| C12 | missingness: 1.23+-1.38 missing, 4.16+-7.64 zeros per descriptor | supported | 1.23 | 1.233 |  |  | 2 |
| C13 | 171 dogs from 117 handlers, 1.46+-0.65 per handler | supported | 1.46 | 1.462 |  |  | 1-2 |
| C14 | network structure varied between dog types | not-reproducible |  |  |  | 0.2574 | 5-6 |
| C15 | gamma=0.65 is the most conservative non-empty value | supported | 0.65 | 0.65 |  |  | 3 |
| C16 | 8 descriptors removed for handler non-independence | supported-with-caveat |  |  |  |  | 2 |
| C17 | 15 imputed datasets averaged and rounded to one | supported-with-caveat |  |  |  |  | 2 |
| C18 | edges are 'credibly non-zero partial correlations' | supported-with-caveat |  |  |  |  | 6 |
| C19 | 8 descriptors failed a handler-independence LRT | not-assessable |  |  |  |  | 2-3 |
| C20 | 95% CIs on Cliff's delta (figure 3) | not-assessable |  |  |  |  | 5 |

Full `reasons` text for every row is in `claim_ledger.csv`.

## 4. Design and sampling - F1

![F1]({{artifact:art_3486fd73-5b24-44fa-9f67-9f244654f99b}})

227 handlers were invited, 174 surveys came back, and after removing three that covered more than one dog the analysis rests on 171 dogs rated by 117 handlers (1.46 +- 0.65 surveys each) - all three figures reproduce exactly (C13). One dog was dropped for zero-coded responses, leaving 116 patrol and 54 detection dogs, and these are the n values passed to `EBICglasso`.

The unit question matters here because 44 handlers rated more than one dog (97 of 170 ratings, up to 3 dogs from one handler). The authors screened for this and removed eight descriptors. Among the repeat handlers, 8 of the 20 **retained** descriptors still have an ICC above 0.10 (largest 0.28, FDA; median 0.09; `tables/handler_icc.csv`). A likelihood-ratio test that fails to reject is not a certificate of independence, and the screen was run on the same averaged dataset the networks use. Against that, the inflation is bounded by how few handlers repeat: with 1.45 ratings per handler overall the design effect is 1.12 and the effective n is about 152 rather than 170 (C16).

Three design questions are not answerable from any file supplied: whether the 117 responding handlers represent the 227 invited (no comparison of responders and non-responders is reported), whether ratings of different dogs by the same handler are conditionally independent given the descriptors retained, and whether the patrol/detection contrast is confounded with breed - the patrol dogs are overwhelmingly German shepherds and the detection dogs are not, so any network difference between dog types is also a difference between breed mixes.

## 5. Method and assumptions - F2, F3

![F2]({{artifact:art_f55a0677-ddc3-453f-88b1-469d4a706d51}})

![F3]({{artifact:art_aa6b8775-dbdc-45ce-9bbb-941ff6a50c87}})

The estimation pipeline reproduces exactly. Polychoric correlation matrices, projected to the nearest positive-definite correlation matrix, then `EBICglasso` at gamma = 0.65 over a 100-value lambda path, give **55 edges (28.95% of 190 pairs) for patrol and 70 (36.84%) for detection** - the paper's 55 (28.95%) and 70 (36.84%) (C1, C2). The centrality orderings reproduce too: PLA highest betweenness and CUR highest strength in patrol, PLA highest on both in detection (C5).

Two points about method follow from the data rather than from the text.

The composite the paper compares between networks - 'the overall mean of descriptors' mean betweenness and strength values' - averages two quantities on scales that differ by about two orders of magnitude. Strength here runs roughly 0.2-1.3 and betweenness 0-53, so the composite correlates with betweenness at r = 0.9999 (patrol) and 0.9999 (detection), against 0.77 and 0.74 for strength. Every statement the paper makes about 'mean centrality', including the Cliff's delta figure, is in practice a statement about betweenness alone. Standardising each measure before averaging would resolve this in one line.

The multiple imputation is disclosed and was reproduced from the script: 15 completions, averaged, rounded to integers, analysed as one dataset. That is the step at which imputation uncertainty leaves the analysis - every interval downstream is conditional on one synthetic dataset. The 15 completions were not supplied, so this review can name the consequence but cannot size it (C17).

On distributional shape: F3 describes the centrality distributions used in the C3/C4 comparison and nothing more. The undesirable group has 6 nodes; a normality test at that n has almost no power, so neither its Shapiro-Wilk p-value nor the Q-Q curvature is treated here as passing or failing anything. The reason the three routes are reported side by side below is precisely that the shape cannot settle the question.

## 6. Resolution and evidence - F5

![F5]({{artifact:art_4c85548c-49be-46bb-953a-f5c7f12296c0}})

The headline centrality claim was run through all three routes in both networks (`tables/three_routes_*.csv`); none was selected after the fact.

| network | measure | desirable mean | undesirable mean | Welch's t-test | Student's t-test | Mann-Whitney U | label permutation (mean difference) | smallest p this design could give |
|---|---|---|---|---|---|---|---|---|
| patrol | strength | 0.907 | 0.492 | 0.0003562 | 0.0007874 | 0.004644 | 0.00052 | 5.16e-05 |
| patrol | betweenness | 18.64 | 3.667 | 0.006533 | 0.03667 | 0.01808 | 0.03424 | 5.16e-05 |
| patrol | mean_cent | 9.775 | 2.079 | 0.005812 | 0.03396 | 0.01538 | 0.03129 | 5.16e-05 |
| detect | strength | 1.003 | 0.713 | 0.04383 | 0.02052 | 0.04076 | 0.02063 | 5.16e-05 |
| detect | betweenness | 14.71 | 9.333 | 0.2882 | 0.3596 | 0.3859 | 0.3802 | 5.16e-05 |
| detect | mean_cent | 7.859 | 5.023 | 0.2741 | 0.3418 | 0.397 | 0.3564 | 5.16e-05 |

The permutation floor for 14 nodes against 6 is p >= 5.16e-05, far below 0.05, so the design can in principle produce a significant result; the floor does not bind. The permutation resolution at 100,000 draws is 1.00e-05.

In patrol every route agrees that desirable descriptors are more central. In detection every route agrees that they are not, on betweenness and on the composite, while strength alone separates them. The routes do not disagree with each other anywhere, so no claim here is `assumption-dependent`; what the claim depends on is the centrality measure, not the test.

One caveat applies to all six tests and is not fixable by any choice of route: the units are the 20 nodes of a single estimated network, which are not independent of one another - they are functions of the same 190 correlations. The p-values describe how surprising the desirable/undesirable split is under random relabelling of nodes; they are not statements about dogs.

![C5 panel]({{artifact:art_12b749e6-ad53-4431-8b7a-ac5401313465}})

The paper withholds centrality confidence intervals, citing known bias in bootstrap intervals for network parameters - a defensible and explicitly argued choice. The consequence is that the centrality ordering in its figure 2 carries no stated uncertainty. Recomputing the bootstrap spread (100 and 100 replicate networks) shows the top node holds first place in 30% of the 100 non-empty patrol replicates on betweenness and 4% of the 93 non-empty detection replicates, while neighbouring nodes' intervals overlap heavily. 'Playful is the most central descriptor' is supported as a point estimate and is not separated from its neighbours by the resampling the paper itself ran (C5).

Two properties of the bootstrap itself qualify that number and apply equally to the paper's own 2000 replicates. Bootstrap networks are systematically denser than the network they resample - mean 75 edges against 55 observed for patrol and 96 against 70 for detection - because sampling dogs with replacement duplicates rows and inflates the correlations; and 7 of 100 detection replicates produced no network at all. A bootstrap replicate is therefore not a re-draw of the published network, and the rank of a node across replicates is a statement about a denser family of graphs than the one in figure 1.

![C8 panel]({{artifact:art_9c0721fb-44ee-4bcc-aa7d-ba0e25425b23}})

Edge counts across bootstrap networks have s.d. 21.76 (patrol) and 30.52 (detection) here against 12.82 and 29.75 reported; the detection network's instability is reproduced, though the replicate count differs (C9). Node-dropping is the weaker of the two stability analyses - it resamples nodes from the same dogs, so a high correlation at three nodes says the centrality ordering is internally consistent, not that it would survive a new sample of dogs. The subject-wise curves are the ones that address sampling, and they fall away well before the full sample in both networks, detection first (C8, C10, C11).

The gaps in the left panel are themselves a result. A level is plotted only where the resampled network has at least three nodes with varying centrality; below that the correlation with the original network is undefined because the resampled graph is empty or near-empty. No level below 0.65 is plottable for patrol and none below 0.95 for detection (levels attempted but unusable - patrol: 0.25, 0.33, 0.41, 0.49, 0.57; detection: 0.25, 0.33, 0.41, 0.49, 0.57, 0.65, 0.73, 0.81, 0.89). The detection network stops producing a network at all once a modest fraction of its 54 dogs is removed, which is a stronger version of the same point as the gamma sweep in section 5.

![C7 panel]({{artifact:art_4fb32135-3b4d-44bd-8056-470a38a5f355}})

Cliff's delta in the paper is computed between the patrol and detection **bootstrap distributions** of centrality, 2000 values against 2000. Those are not 2000 independent observations - they are 2000 re-descriptions of 116 and 54 dogs. The consequence is visible directly: redrawing delta from subsets of the bootstrap replicates, its point estimate stays near +0.796 while its spread falls from 0.092 at 20 replicates to 0.020 at 80. A confidence interval computed this way narrows with computation rather than with evidence, and the 95% intervals drawn in the paper's figure 3 inherit that property (C7).

![C14 panel]({{artifact:art_94a93d7f-d4a6-47e6-acd4-67a90dd31443}})

The patrol-versus-detection comparison is the paper's main substantive contrast and it is never formally tested: two networks are estimated separately, their centralities are compared through bootstrap effect sizes, and the structures are compared in prose. A label-permutation test of structure invariance - permute dog type, re-estimate both networks, record the largest absolute edge-weight difference - gives an observed 0.439 against a null 95th percentile of 0.485, p = 0.257 (100 permutations; resolution 0.0099, so no p below that value is reportable from this budget). Under a max-statistic correction 0 of the 190 edges and 0 of the 20 node strengths differ at FWER 0.05 (`tables/nct_edge_differences.csv`). Note that this permutation treats dogs as exchangeable and so ignores the repeated-handler structure described in section 4; a handler-stratified permutation would be the stricter version (C14).

## 7. Multiplicity - F6

![F6]({{artifact:art_77af8acb-44d6-4aa8-bb50-071f50528b86}})

![C18 panel]({{artifact:art_ff8468a2-a757-478d-baac-19fa4e37ec34}})

**The paper reports no p-values, so F6 is not a picture of the authors' results.** It is built from this review's own recomputation of the pairwise tests the network estimation implies, and the p-hacking diagnostic that normally accompanies it is deliberately absent: with zero reported p-values there is nothing to describe, and the just-below-threshold count says nothing about any individual result even when there are enough of them.

The multiplicity that does exist is structural. Estimating a network over 20 nodes decides 190 node pairs; two networks decide 380. Re-running those pairs as zero-order polychoric correlation tests:

| set | uncorrected | Bonferroni | BH | expected by chance at 0.05 | edges the lasso selected |
|---|---|---|---|---|---|
| patrol (190 pairs) | 128 | 73 | 121 | 9.5 | 55 |
| detection (190 pairs) | 113 | 74 | 109 | 9.5 | 70 |
| both networks (380) | 241 | 136 | 228 | 19 | 125 |

Because the 190 pairs are heavily dependent, Bonferroni is the wrong correction for them; the max-statistic permutation null is the one that respects the dependence. Permuting each descriptor independently (300 permutations per network) puts the FWER-0.05 threshold at |r| = 0.555 for patrol and 0.819 for detection (`tables/maxstat_*.csv`).

None of this is a charge of uncorrected testing, because the paper never claims an error rate: the lasso is a model-selection procedure and its output is a selected graph, not a set of rejected nulls. The point is the one the paper's own wording invites - edges are described as 'credibly non-zero partial correlations', and nothing in the procedure attaches a credibility level to any of them. Reporting the selection as what it is, and giving an edge-level stability proportion across bootstrap networks, would close the gap (C18).

A second kind of multiplicity sits in the comparison between networks: 20 nodes x 2 centrality measures = 40 Cliff's deltas, read for which are 'large', with no correction and no formal test. The max-statistic test in section 6 is the corrected version of that comparison.

## 8. Per-claim supplement

![C3]({{artifact:art_49eff456-f347-41c0-bced-7b65dd3145ed}})

![C4]({{artifact:art_12b92879-08c4-4527-8df7-a177231f8bf5}})

![C2]({{artifact:art_2d8aabb2-a121-4a13-9a4a-50984ab83c1f}})

![C17]({{artifact:art_d322db2c-9d2e-44a6-96ee-7d226576392d}})

![statcheck]({{artifact:art_664cdaea-453d-456a-a627-173eeea859b4}})

Panels for C5, C7, C8 and C14 appear in section 6 and are not repeated. The two `not-assessable` panels are kept rather than dropped: an absent diagnostic and an impossible one look identical in a report that simply omits the slot.

## 9. What was not assessed

```
Diagnostics in this report describe distribution shape and spread only. They cannot evaluate:
  - independence of observations
  - the sample is representative of the stated population
  - absence of unmeasured confounding
  - measurements were not selected after seeing the outcome
These are design questions. Assess them from the methods section, and say so explicitly when the methods do not contain enough information to assess them.
```

Specific to this paper:

- **Representativeness.** 117 of 227 invited handlers responded. Nothing in the supplied files allows responders to be compared with non-responders, so whether these 170 dogs represent Norwegian police dogs is unknown. Every claim in the ledger is conditional on this sample.
- **Breed confounding with dog type.** Patrol dogs are almost all German shepherds and Belgian malinois; detection dogs are mostly labradors, retrievers and spaniels. The network difference tested in C14 is a difference between work roles and between breed compositions at the same time, and no supplied file can separate them.
- **Conditional independence of repeated ratings.** The ICCs in section 4 describe between-handler variation; they cannot establish that the retained descriptors are independent, only that eight of the twenty still show appreciable clustering.
- **The authors' own likelihood-ratio screen.** Neither the paper nor the script records the per-descriptor LRT or delta-AIC values, and `clmm` was not re-run here, so the screen itself is reproduced in description only.
- **Imputation uncertainty (C17).** Only the rounded average of the 15 completions was supplied.
- **Supplementary tables S1-S5 and the authors' bootstrap tables** were not supplied; the reported raw centrality values (table S4) and Cliff's delta values (table S5) could be checked only against the nine deltas quoted in the running text.
- **Absolute betweenness scale.** The ordering of betweenness reproduces exactly, but `qgraph` and the implementation here may differ by whether shortest paths are counted once or twice per node pair (a factor of two). Table S4 would settle it; the comparison of nodes within a network is unaffected.

## 10. Reproduction

- `review_goold_2016.py` regenerates every figure, table and this report; `core.py` holds the numerics (polychoric ML, Higham nearest-correlation projection, block-coordinate-descent graphical lasso, EBIC selection, Brandes betweenness); `prep.py` re-derives the authors' preprocessing from `rawData.csv`; `heavy.py` produces the resampling.
- Seeds: three-routes and permutation tests `random_state = 1`; resampling seeds are stored inside each `results/*.npz` (`seed0` field).
- **Reduced resampling budget.** The paper used 2000 bootstrap samples per analysis. Shared-machine limits forced a smaller budget here: 100 non-parametric bootstrap networks per group, 100 label permutations for the network-difference test, 300 permutations for the edge max-statistic null, and 6 replicates per level for the case-dropping and node-dropping curves. The lambda path in resampling uses 25 values with `lambda.min.ratio = 0.05` against 100 and 0.01 for the two reported networks; on the observed data that coarser grid returns the identical patrol network (55 edges) and 71 rather than 70 edges for detection. Consequences, stated rather than hidden: the permutation p-values in C14 cannot resolve below 0.0099 and the max-statistic thresholds below 0.0033, so any claim needing a smaller p-value is reported as unresolvable at this budget rather than given a number; bootstrap s.d.s (C9) and stability thresholds (C10) are quantitative at low precision and are labelled accordingly; the Cliff's delta demonstration (C7) is about how the spread scales with replicate count and is unaffected.
- Three-route and permutation tests use n_perm = 100,000 (resolution 1.00e-05).
- Environment: Python, numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, matplotlib 3.11.2. `scikit-learn`, `networkx` and R were unavailable in this sandbox, so the graphical lasso, the nearest-correlation projection and weighted betweenness were implemented directly; the lasso solution was checked against the KKT conditions of its own objective (maximum stationarity residual ~1e-5) and the whole pipeline against the paper's printed edge counts (55 and 70) and centrality ordering.
- Input checksums are in section 1. Tables behind every number are in `tables/`.
