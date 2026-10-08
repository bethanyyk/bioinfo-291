# Immune response signatures in melanoma checkpoint blockade

**Single-cell RNA-seq analysis of melanoma tumour biopsies before and during immune checkpoint blockade (anti-PD-1), asking which immune populations expand or contract with treatment and whether any baseline signature stratifies responders.**

**Data:** GEO [GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575) — Sade-Feldman et al., *Cell* 2018 (PMID 30388456). 16,291 CD45+ cells, Smart-seq2, from 48 biopsies of 32 patients: 19 baseline and 29 on-treatment, each annotated responder or non-responder. Raw reads are controlled-access (dbGaP phs001680.v1.p1); this analysis uses only the public processed TPM matrix and its annotation.

## Quick start

```bash
# Set up the environment (one-time)
conda env create -f environment/environment.yml
conda activate icb-scrna

# Run the full pipeline
python workflows/run_pipeline.py --config configs/analysis.yaml

# Validate with tests
pytest tests/

# See SETUP.md for detailed instructions, troubleshooting, and how to extend the analysis
```

**Runtime:** ~20 minutes on 8 cores, ~1.5 GB peak RAM  
**For detailed setup and workflows:** see [SETUP.md](SETUP.md)  
**For design decisions and reproducibility framework:** see [DESIGN.md](DESIGN.md)

## Project layout

```
├── configs/analysis.yaml         # All parameters
├── src/icb_scrna/                # Analysis modules (ingest, preprocess, annotate, etc.)
├── workflows/run_pipeline.py     # Main entry point
├── tests/                        # 32 regression tests (all passing)
├── data/
│   ├── raw/                      # Downloaded GEO files
│   ├── processed/                # Intermediate .h5ad files
│   └── metadata/                 # Manifests, annotations, checksums
├── results/                      # Analysis tables
├── figures/                      # Publication-grade figures
├── reports/analysis_report.md    # Full methods, findings, limitations
├── environment/                  # Conda spec + lockfile (106 packages)
├── CLAUDE.md                     # Design decisions, data contracts, constraints
└── README.md                     # This file
```

## Main findings

### No cell-type abundance changes with treatment
Across 48 biopsies (11 paired patients, 29 on-treatment samples), no cell type differs between baseline and on-treatment at FDR < 0.05 (smallest adjusted p = 0.44). This suggests gross immune composition is not where the treatment effect lives.

### Baseline composition separates responders
7 of 11 populations differ between responders and non-responders at baseline (FDR < 0.05):

| Population | Responder mean | Non-responder mean | OR | FDR |
|---|---|---|---|---|
| B cells | 0.137 | 0.020 | 10.14 | 0.0016 |
| CD4 T | 0.332 | 0.145 | 3.18 | 0.0001 |
| CD8 T | 0.327 | 0.473 | 0.49 | 0.0105 |
| Cycling T | 0.021 | 0.080 | 0.20 | 0.0001 |
| pDC | 0.006 | 0.025 | 0.20 | 0.0016 |
| Monocyte/Macrophage | 0.018 | 0.078 | 0.16 | 0.0040 |
| cDC | 0.006 | 0.020 | 0.25 | 0.0078 |

Responders had *lower* CD8 at baseline and higher B + CD4, suggesting the immune composition in responders is shifted toward lymphocytes, not uniformly more lymphoid.

### On-treatment CD8 T cells carry a chronic interferon program
109 genes at FDR < 0.1 are differentially expressed in CD8 T cells comparing non-responders vs. responders on-treatment. All but one (108) are elevated in non-responders, dominated by interferon-stimulated transcripts: *IFI6*, *HLA-DPA1*, *CD38*, *IFITM1*, *IFITM2*, *UBE2L6*, *RARRES3*, *PSME1*.

This interferon signature in non-responders is **independently** found in the baseline population analysis: on-treatment CD8 T cells (which are enriched in non-responders) form a distinct cluster that is **absent** in baseline samples, and a synthetic "CD8 interferon state" distinguishes responders at baseline (OR 0.46, FDR 0.011).

### No baseline signature predicts response beyond chance
Five candidate signatures were evaluated by leave-one-patient-out cross-validation:
1. **Whole-compartment composition** — AUC 0.689
2. **CD8 state composition** — AUC 0.644
3. **TCF7+ CD8 fraction** (published comparator) — AUC 0.689
4. **CD8 baseline gene program** — AUC 0.711
5. **CD8 interferon response program** — AUC 0.733 (best)

The best signature reaches AUC **0.733** with 95% CI **0.48–0.96** and permutation p = **0.055**. Every confidence interval includes 0.5 (random performance). At n = 19 baseline patients, power is insufficient to distinguish any signature from chance.

## Key decisions documented in CLAUDE.md

1. **Biopsy as the unit of inference** — Response is annotated per biopsy. Three patients show discordant labels across lesions, so patient-level averaging would be wrong.

2. **Why DESeq2 was not used** — GEO publishes only TPM; raw reads are dbGaP-controlled. The negative-binomial model is misspecified for TPM. Expression comparisons use patient-level pseudobulk with Welch's t-test.

3. **Why binomial GLM with robust errors** — A standard approach for testing proportions. Propeller is R-only, but the GLM with patient-clustered errors achieves the same goal.

4. **Windows sandbox constraints** — harmonypy pinned to 0.0.9 (later releases import torch, whose DLLs fail to load); pytest rootdir discovery fails on `C:\` (workaround: run tests from outside the drive-root ancestry).

## Methods in brief

**Ingest & QC:** Genes detected in <10 cells dropped. Cells filtered at ≥500 genes and ≤20% mitochondrial TPM. `log2(TPM+1)` rescaled to natural log. 14,503 of 16,291 cells retained (89%); all loss on the 20% mitochondrial threshold.

**Integration & clustering:** 2,000 HVGs, 30 PCA components, Harmony integration over patient, Leiden clustering (resolution 1.0). One NK/T mixed cluster was sub-clustered. 11 populations assigned from canonical markers (CELLxGENE CellGuide).

**Cell filtering:** 991 cells from lineage-enriched sorts (recorded only in the expression header's second row) are excluded from all composition analyses to avoid sort-driven proportions. Analysis uses 14,503 unbiased CD45+ cells.

**Abundance:** Binomial GLM per biopsy with patient-clustered robust standard errors (not propeller; R-only).

**Expression:** Pseudobulk aggregated to patient level, compared with Welch's t-test (raw reads unavailable; count models inapplicable to TPM).

**Signatures:** Leave-one-patient-out CV with gene selection inside each fold, bootstrap CIs (2,000 replicates), and permutation null (200 shuffles).

## Reproducibility

- **Environment:** All dependencies pinned in `environment.lock.txt` (106 packages)
- **Configuration:** Parameters live in `configs/analysis.yaml`, not hardcoded
- **Tests:** 32 regression tests covering parsing, design matrices, and signature validation (all passing)
- **Reproducibility:** Full re-run from clean data reproduces all outputs exactly

Run the test suite to verify:
```bash
pytest tests/ -v
```

## Limitations

1. **Small sample size** (n = 19 baseline). Signature evaluation is underpowered; no signature clears the permutation null.

2. **Biopsy-level labels** with patient-level inconsistency. Three patients show discordant response labels across lesions; all tests use biopsy-level grouping with patient-clustered errors.

3. **No raw-read data** in public GEO (dbGaP-controlled). Analysis uses published TPM; count-based methods (DESeq2) are misspecified and not used.

4. **No validation cohort.** Findings describe the training set; no external validation.

5. **Response labels only** — binary outcome (responder/non-responder). No survival data, progression-free survival, or continuous measures of response.



## Output files

**In `results/` (committed to repo):**
- `differential_abundance.csv` — Abundance contrasts (baseline vs. on-treatment, responder vs. non-responder)
- `cluster_annotation_evidence.csv` — Marker gene evidence for each cell type assignment
- `pseudobulk_de_significant.csv` — Differential expression at FDR < 0.1
- `markers_cell_type_top.csv` — Top 25 markers per cell type
- `signature_performance.csv` — AUC, CI, and permutation p per signature
- `signature_scores_loo.csv` — Leave-one-patient-out signature scores

**In `data/processed/` (gitignored, regenerable):**
- `gse120575_raw.h5ad` — 16,291 cells, 36,602 genes
- `gse120575_annotated.h5ad` — Annotated with cell types and metadata

**In `figures/` (committed):**
- `fig0_qc.png` — Quality control thresholds (genes, mitochondrial fraction)
- `fig1_atlas.png` — UMAP of all 11 immune populations
- `fig2_composition.png` — Baseline and on-treatment proportions; responder vs. non-responder
- `fig3_cd8_volcano.png` — Differential expression in on-treatment CD8 T (responders vs. non-responders)
- `fig4_signatures.png` — Signature performance (ROC curves, composition-based predictions)

---

**Last updated:** October 2026

For installation and development instructions, see `../0-reproducible-skills/README.md`.
