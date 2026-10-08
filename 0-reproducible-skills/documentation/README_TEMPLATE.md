# [Project Title]

**One-sentence summary:** [What does this analysis do? What question does it answer?]

**Data:** [Source (GEO, dbGaP, etc.), accession, link, sample size, key attributes]

**Output:** [One figure / table summarizing the main findings]

## Running the analysis

```bash
# 1. Set up the environment (one-time)
conda env create -f environment/environment.yml
conda activate <env-name>

# 2. Run the full pipeline (default: skips completed stages)
python workflows/run_pipeline.py --config configs/analysis.yaml

# 3. Or run one stage
python workflows/run_pipeline.py --config configs/analysis.yaml --stage ingest

# 4. Or re-run everything from scratch
python workflows/run_pipeline.py --config configs/analysis.yaml --force

# 5. Verify with tests
pytest tests/ -v
```

**Runtime:** [Expected wall-clock time for a full run]

**Memory:** [Peak RAM usage, e.g., "1.5 GB"]

## Project layout

```
├── 0-reproducible-skills/        # Reusable tools (separate from this project)
├── 1-scRNA/                       # This project
│   ├── configs/
│   │   └── analysis.yaml          # All parameters, thresholds, paths
│   ├── src/
│   │   └── <package>/             # Analysis package (one module per stage)
│   │       ├── ingest.py
│   │       ├── preprocess.py
│   │       ├── annotate.py
│   │       └── ...
│   ├── workflows/
│   │   └── run_pipeline.py        # Stage orchestration driver
│   ├── tests/
│   │   ├── conftest.py            # Shared fixtures
│   │   ├── test_ingest.py
│   │   └── ...
│   ├── data/
│   │   ├── raw/                   # (gitignored) external downloads
│   │   ├── processed/             # (gitignored) intermediate .h5ad, .parquet
│   │   └── metadata/              # (committed) small metadata files
│   ├── results/
│   │   ├── *.csv                  # (committed) analysis tables
│   │   └── *.csv.gz               # (committed, or LFS) large result matrices
│   ├── figures/
│   │   └── *.png                  # (committed) publication-grade figures
│   ├── reports/
│   │   └── analysis_report.md     # Findings, methods, limitations, caveats
│   ├── environment/
│   │   ├── environment.yml        # Conda spec
│   │   └── environment.lock.txt   # Exact versions from the successful run
│   ├── CLAUDE.md                  # Implementation notes & design decisions
│   ├── README.md                  # This file
│   └── pytest.ini                 # Test configuration
```

## What was analyzed

### Stage: [Stage Name]
[Brief description of what this stage does, why it matters]

**Input:** [File path, format, description]

**Output:** [File path(s), what each contains]

**Key parameters:** [Names and values of parameters that drove this stage]

**Validations:** [What was checked to ensure correctness?]

### Example from bioinfo-291:

### Stage: ingest
Parses the public expression matrix (TPM) and cell metadata from GEO. The matrix has a
two-row header and special encoding; the metadata file includes a 19-line preamble. Both
are validated against expected shape and presence of required columns.

**Input:**
- `GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz` (127 MB)
- `GSE120575_patient_ID_single_cells.txt.gz` (83 KB)

**Output:**
- `data/processed/gse120575_raw.h5ad` — 16,291 cells × 36,602 genes

**Key parameters:**
- Min genes per cell: 500
- Min cells per gene: 10

**Validations:**
- Cell count matches GEO metadata (16,291 cells)
- No duplicate cell barcodes
- No missing genes or cells in expression matrix

## Methods in brief

[1-2 paragraph methods section with exact parameter values, explaining how each stage works]

### Example:
Genes detected in <10 cells dropped. Cells filtered at ≥500 genes, ≤20% mitochondrial TPM.
`log2(TPM+1)` rescaled to natural log. 2,000 HVGs selected by variance; 30 PCA components;
Harmony integration over patient covariate; Leiden clustering (resolution 1.0) on 30 PCs.
Cell types assigned by marker expression (canonical markers from CellGuide) and cross-checked
against flow-sorted reference signatures.

Abundance tested per biopsy (unit of inference) with binomial GLM and patient-clustered robust
errors. Expression tested on patient-level pseudobulk (no raw counts available) with Welch's
t-test. Signatures evaluated by leave-one-patient-out CV with gene selection inside each fold
and 200 permutation nulls.

## Main findings

[Bullet list of key results]

### Example:
- **No treatment-induced cell-type abundance change** (FDR < 0.05); smallest adjusted p = 0.44
- **Baseline composition separates responders** — 7 of 11 populations differ (B cells 10× higher in responders)
- **On-treatment CD8 T cells show an interferon program in non-responders** (109 genes, FDR < 0.1)
- **No baseline signature predicts response beyond chance** (best AUC 0.73, 95% CI 0.48–0.96, permutation p = 0.055)

## Limitations

[Known caveats, sample size constraints, alternative interpretations, generalizability questions]

### Example:
- Small sample size (n=19 baseline patients). Signature evaluation is underpowered.
- Biopsy-level response labels. Three patients show discordant response across lesions; inference is biopsy-level, not patient-level.
- No raw sequencing data in public GEO; analysis uses published TPM. (Raw reads are dbGaP-controlled; count-based tests would require access.)
- No validation cohort. Findings are described, not predictive.

## For more details

- **Full report:** `reports/analysis_report.md` — extended methods, all results, supplementary analyses
- **Implementation notes:** `CLAUDE.md` — design decisions, data format contracts, platform constraints
- **Source code:** `src/` — each module is self-contained and can be explored independently
- **Test suite:** `tests/` — regression tests for data parsing, statistical correctness, and pipeline reproducibility

## Extending or modifying this project

1. **Running a different analysis:** Edit `configs/analysis.yaml` (parameters stay in one place)
2. **Adding a stage:** Create a new module in `src/`, add it to the stage list in `workflows/run_pipeline.py`, and add tests
3. **Updating the environment:** Edit `environment/environment.yml`, re-run `conda env create`, commit the new `environment.lock.txt`
4. **Regenerating outputs:** Run with `--force` or delete the output file to re-trigger a stage

Every change should be committed. If a result differs, `git diff` will show what changed.

## Reproducibility

- All parameters are in `configs/`, not hardcoded.
- Environment is pinned (`environment.lock.txt`); exact replay is possible.
- Test suite validates data parsing and statistical correctness.
- `CLAUDE.md` documents design decisions and constraints.
- A full re-run from clean data reproduces all outputs (verified by spot-checking key results).

---

**Last updated:** [Date]

**By:** [Your name]
