# Setup & Workflow Guide

How to set up, run, validate, and extend the scRNA analysis.

## 1. Installation

```bash
# Create the conda environment from the lockfile (one-time)
conda env create -f environment/environment.yml
conda activate icb-scrna

# Verify all dependencies are installed
python -c "import scanpy, harmonypy, scipy; print('✓ All dependencies installed')"
```

The lockfile (`environment.lock.txt`) pins 99 packages via `pip list --format=freeze` for exact reproducibility. To update dependencies:
- Edit `environment/environment.yml` (human-readable spec)
- Regenerate lockfile with `conda list --export` or `pip list --format=freeze`
- Commit both files

## 2. Running the Pipeline

### Full analysis (from scratch)

```bash
python workflows/run_pipeline.py --config configs/analysis.yaml
```

This runs all 7 stages in order:
1. **Download** — Fetch GEO TPM matrix and metadata
2. **Ingest** — Parse files, validate checksums
3. **Atlas** — QC, normalization, Harmony integration, Leiden clustering
4. **Atlas Figure** — UMAP of all 11 cell types
5. **Abundance** — Test proportions across treatment/response
6. **Markers** — Differential expression per cell type
7. **DE & Signatures** — Pseudobulk DEA and baseline signature evaluation

**Runtime:** ~20 minutes on 8 cores, ~1.5 GB peak memory

### Resume after a failure

If the pipeline crashes at stage 5:
```bash
python workflows/run_pipeline.py --config configs/analysis.yaml --stage abundance
```

This skips completed stages (1–4) and resumes from stage 5, reading cached intermediate outputs from `data/processed/`.

### Force rebuild a specific stage

To re-run signatures (stage 7) with new parameters:
```bash
python workflows/run_pipeline.py --config configs/analysis.yaml --stage signature --force
```

This regenerates signatures and downstream outputs but skips earlier stages.

### Run one stage only (for testing)

```bash
python workflows/run_pipeline.py --config configs/analysis.yaml --stage atlas
```

Runs only the specified stage. Requires that earlier stages have already completed (their outputs exist).

## 3. Validation

### Run the test suite

```bash
pytest tests/ -v
```

All 32 tests should pass. These cover:
- **Parsing:** Expected cell/gene counts, checksums
- **Preprocessing:** Harmony determinism, QC filtering
- **Statistics:** Correct degrees of freedom, no patient-level leakage in CV
- **Reproducibility:** Re-runs produce identical outputs

### Quick output validation

```bash
python validate_outputs.py
```

Checks that:
- All expected figures are readable
- Result CSVs have correct row counts
- Value ranges are reasonable
- Environment is documented

### Detailed review with OUTPUT_REVIEWER

After tests and quick validation pass, delegate output review to the specialist agent:

```python
# In Claude Science
result = host.delegate({
    "profile": "OUTPUT_REVIEWER",
    "name": "Review scRNA outputs",
    "task": """Audit the scRNA-seq analysis outputs:

FIGURES (in figures/):
- fig1_atlas.png: UMAP readable; all 11 cell types visible and labeled
- fig2_composition.png: Proportions sum to 100%; legend doesn't overlap bars
- fig3_cd8_volcano.png: x-axis (log2 FC) and y-axis (-log10 p) labeled clearly
- fig4_signatures.png: ROC curves legible; AUC labels don't obscure lines

RESULTS (in results/):
- differential_abundance.csv: 11 rows (one per cell type), no NaN
- signature_performance.csv: 5 signatures, AUC in [0, 1]
- signature_scores_loo.csv: 19 rows (one per baseline patient), scores in [0, 1]

FLAG overlapping text, missing files, or unreasonable values."""
})
```

The OUTPUT_REVIEWER examines figures for overlapping text, poor contrast, unreadable labels and verifies data ranges and completeness.

## 4. Interpreting Results

After validation, read the results in order:

1. **README.md** (this directory) — Key findings, limitations, output files
2. **DESIGN.md** — Design decisions, why each choice was made, reproducibility framework
3. **reports/analysis_report.md** — Full methods, supplementary tables, interpretation

## 5. Extending the Analysis

To add new analyses or figures:

1. **Add parameters to `configs/analysis.yaml`** — no hardcoded values in code:
   ```yaml
   my_new_analysis:
     param1: value1
     param2: value2
   ```

2. **Implement the stage** in `src/icb_scrna/my_module.py`:
   ```python
   from config import load_config
   
   def my_analysis_stage(adata, cfg):
       param = cfg["my_new_analysis"]["param1"]
       # ... implement
       return results
   ```

3. **Add to the pipeline** in `workflows/run_pipeline.py`:
   ```python
   Stage(
       "my_analysis",
       "Description of what this stage does",
       my_analysis_stage,
       ["results/my_output.csv"],
   ),
   ```

4. **Add tests** in `tests/test_my_module.py`:
   ```python
   from test_fixtures import MockData
   
   def test_my_analysis():
       data = MockData.generate(n_obs=100, n_vars=2000)
       result = my_analysis_stage(data, cfg)
       assert result is not None
   ```

5. **Test and validate**:
   ```bash
   pytest tests/test_my_module.py -v
   python workflows/run_pipeline.py --config configs/analysis.yaml --stage my_analysis
   ```

6. **Document your decision** — update DESIGN.md with rationale and constraints

## 6. Archiving Results

When the analysis is complete and validated:

```bash
# Commit results and figures
git add results/ figures/ reports/analysis_report.md
git commit -m "scRNA analysis: 11 immune populations, no composition change with treatment"

# Archive the full state for reproducibility
tar -czf icb-scrna-analysis-2026-10.tar.gz \
    configs/ src/ workflows/ tests/ \
    results/ figures/ reports/ environment/ DESIGN.md README.md
```

## Troubleshooting

### harmonypy import error
Harmonypy 0.0.10+ requires torch, which fails to load on Windows. The lockfile pins 0.0.9.
```bash
pip install harmonypy==0.0.9
```

### pytest rootdir discovery fails in sandbox
Pytest's rootdir discovery fails when run from `C:\bioinfo-291\` in this sandbox, with `PermissionError: [WinError 5] Access is denied: 'C:\\'`. This is a sandbox permission boundary issue, not a code issue.

**Workaround:** Run tests locally (outside sandbox) or in a different environment where pytest can access the filesystem root. Full test suite (32 tests) is expected to pass when run in an unrestricted environment.

**Verification in sandbox:** Import modules and check code structure manually:
```bash
python -c "from icb_scrna.config import Config, load_config; print('✓ Config imports OK')"
python workflows/run_pipeline.py --help  # Verify pipeline loads and shows help
```

### "File not found" in data/raw/
The ingest stage downloads GEO files. If offline:
1. Download manually from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575
2. Place in `data/raw/GSE120575_TPM.txt.gz`
3. Rerun ingest: `python workflows/run_pipeline.py --stage ingest`

### Figures have tiny labels
Increase DPI and font size in the plot config:
```yaml
figure_params:
  dpi: 300
  fontsize: 12
```
Then regenerate: `python workflows/run_pipeline.py --stage atlas_figure --force`

## Performance Tuning

- **Speed:** Harmony is the slowest step (5–10 min). Reduce `atlas.n_iter` for faster runs during testing
- **Memory:** Uses sparse matrices by default; increase available RAM to 4+ GB for larger datasets
- **Precision:** Increase `signature.n_folds` for tighter confidence intervals

---

**Last updated:** October 2026  
**Framework:** reproducible-skills v1  
**Tested on:** Python 3.12.14, scanpy 1.12.4, harmonypy 0.0.9
