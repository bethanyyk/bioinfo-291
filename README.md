# scRNA-seq Immunotherapy Tumor Response Analysis

Single-cell analysis of melanoma tumour biopsies taken **before and during immune checkpoint
blockade**, asking which immune populations expand or contract with treatment, what marks them,
and whether any baseline signature stratifies responders.

**Data:** GEO [GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575) —
Sade-Feldman et al., *Cell* 2018 (PMID 30388456). 16,291 CD45+ cells, Smart-seq2, from 48
biopsies of 32 patients: 19 baseline and 29 on-treatment, each annotated responder or
non-responder. Raw reads are controlled-access (dbGaP phs001680.v1.p1); this project uses only
the public processed TPM matrix and its annotation.

## What the analysis found

- **Treatment does not detectably change gross immune composition** — no cell type differs
  between baseline and on-treatment at FDR < 0.05 across 48 biopsies.
- **Baseline composition does separate responders** — 7 of 11 populations differ before
  treatment, led by B cells (odds ratio 10.1) and CD4 T cells (3.18) in responders.
- **Non-responder CD8 T cells carry a chronic interferon program on treatment** — 109 genes at
  FDR < 0.1, 108 of them elevated in non-responders.
- **No baseline signature predicts response beyond chance** at n = 19; the best reaches
  AUC 0.73 (95% CI 0.48–0.96, permutation p = 0.055).

Full write-up with figures, tables, and limitations: [`reports/analysis_report.md`](reports/analysis_report.md).

## Running it

```bash
conda env create -f environment/environment.yml
conda activate icb-scrna
python workflows/run_pipeline.py --config configs/analysis.yaml
pytest
```

The driver downloads the GEO files (127 MB), builds the atlas, and writes every table and figure.
Stages skip when their outputs exist; `--force` re-runs, `--stage NAME` runs one. About 20 minutes
from scratch on 8 cores, peaking under 2 GB of RAM — the expression matrix is streamed into a
sparse representation rather than loaded dense.

## Layout

```
configs/analysis.yaml      every parameter; nothing hardcoded in src/
src/icb_scrna/             ingest, preprocess, annotate, abundance, de, features, signature, plots
workflows/run_pipeline.py  stage driver
tests/                     32 tests; no network, no data/raw needed
results/                   analysis tables (committed)
figures/                   fig0 QC, fig1 atlas, fig2 composition, fig3 CD8 volcano, fig4 signatures
reports/analysis_report.md findings, methods, limitations
data/                      raw downloads and derived .h5ad (gitignored)
```

Conventions, the directory contract, and the memory-budget constraint on the loader are documented
in [`CLAUDE.md`](CLAUDE.md).

## Methods in one paragraph

Genes detected in <10 cells dropped at ingest; cells filtered at ≥500 genes and ≤20% mitochondrial
TPM (14,503 of 16,291 retained). `log2(TPM+1)` rescaled to natural log. 2,000 HVGs, 30 PCs, Harmony
over patient, Leiden resolution 1.0, with one NK/T mixed cluster sub-clustered; 11 populations
assigned from canonical markers cross-checked against CELLxGENE CellGuide. 991 cells from
lineage-enriched sorts are excluded from all composition analyses. Abundance tested per biopsy with
a binomial GLM and patient-clustered robust standard errors; expression tested on patient-level
pseudobulk with Welch's t-test (DESeq2 is inapplicable — there are no counts). Signatures evaluated
by leave-one-patient-out CV with gene selection inside each fold, bootstrap CIs, and a permutation
null.