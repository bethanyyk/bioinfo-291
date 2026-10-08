# bioinfo-291 Repository Structure & Guide

This document describes the overall organization and how to use each part.

## Directory Structure

```
bioinfo-291/
│
├── 0-reproducible-skills/           # REUSABLE TOOLS (use in any project)
│   ├── README.md                    # Master guide: what's available and how to use it
│   ├── config/
│   │   ├── __init__.py
│   │   └── config_base.py          # YAML config loader with path resolution
│   ├── pipeline/
│   │   ├── __init__.py
│   │   └── pipeline_driver.py      # Stage orchestration with skip logic
│   ├── testing/
│   │   ├── __init__.py
│   │   └── test_fixtures.py        # Mock data, validators, comparators
│   ├── utilities/
│   │   ├── __init__.py
│   │   └── data_validation.py      # Checksums, CSV validation, parse logs
│   ├── environment/
│   │   ├── environment_template.yml     # Template for conda environments
│   │   └── environment_notes.md         # Guide to pinning strategy
│   └── documentation/
│       ├── CLAUDE_TEMPLATE.md      # Template for CLAUDE.md
│       └── README_TEMPLATE.md      # Template for README.md
│
├── 1-scRNA/                         # WORKED EXAMPLE (complete analysis)
│   ├── README.md                    # What the analysis does, findings, how to run it
│   ├── CLAUDE.md                    # Design decisions, data contracts, constraints
│   ├── configs/
│   │   └── analysis.yaml           # All parameters (not hardcoded in code)
│   ├── src/icb_scrna/              # Analysis modules (one per stage)
│   │   ├── __init__.py
│   │   ├── ingest.py               # Load and parse raw data
│   │   ├── preprocess.py           # QC and normalization
│   │   ├── annotate.py             # Cell type assignment
│   │   ├── abundance.py            # Test cell-type proportions
│   │   ├── de.py                   # Pseudobulk differential expression
│   │   ├── features.py             # Extract signature features
│   │   ├── signature.py            # Cross-validated signature evaluation
│   │   ├── plots.py                # Figure builders
│   │   └── config.py               # Config loading helper
│   ├── workflows/
│   │   └── run_pipeline.py         # Main entry point; orchestrates all stages
│   ├── tests/
│   │   ├── conftest.py             # Shared test fixtures
│   │   ├── test_ingest.py          # Tests for ingest module
│   │   ├── test_abundance.py       # Tests for abundance module
│   │   ├── test_config_and_de.py   # Tests for config and de modules
│   │   └── test_signature.py       # Tests for signature module
│   ├── data/
│   │   ├── raw/                    # (gitignored) Downloaded GEO files
│   │   ├── processed/              # (gitignored) Intermediate .h5ad files
│   │   └── metadata/               # (committed) Manifests, annotations, checksums
│   ├── results/
│   │   ├── *.csv                   # (committed) Analysis tables
│   │   └── *.csv.gz                # (committed) Large result matrices
│   ├── figures/
│   │   └── fig*.png                # (committed) Publication-grade figures
│   ├── reports/
│   │   └── analysis_report.md      # Full findings, methods, limitations
│   ├── environment/
│   │   ├── environment.yml         # Conda spec with load-bearing harmonypy pin
│   │   ├── environment.lock.txt    # Exact versions from successful run
│   │   └── environment_notes.md    # Guide to environment management
│   └── pytest.ini                  # Pytest configuration
│
├── README.md                        # Repo-level README (you are here)
├── STRUCTURE.md                     # This file
├── CLAUDE.md                        # Repo-level design notes
├── .gitignore                       # Excludes data/, large results
└── pytest.ini                       # Root-level pytest config
```

## When to use what

### Using `0-reproducible-skills/` in a new project

1. **Copy the base modules:**
   ```bash
   cp 0-reproducible-skills/config/config_base.py       <your-project>/src/
   cp 0-reproducible-skills/pipeline/pipeline_driver.py <your-project>/src/
   ```

2. **Use the templates:**
   ```bash
   cp 0-reproducible-skills/environment/environment_template.yml  <your-project>/environment/environment.yml
   cp 0-reproducible-skills/documentation/CLAUDE_TEMPLATE.md      <your-project>/CLAUDE.md
   cp 0-reproducible-skills/documentation/README_TEMPLATE.md      <your-project>/README.md
   ```

3. **Import the utility modules:**
   ```python
   # Adjust sys.path as needed if reproducible-skills is a sibling
   from config_base import load_config
   from pipeline_driver import PipelineDriver, Stage
   from test_fixtures import MockData, Validator
   from data_validation import DataManifest, ParseLog
   ```

### Understanding the scRNA example (`1-scRNA/`)

Read in this order:

1. **`1-scRNA/README.md`** — What the analysis does, how to run it, main findings
2. **`1-scRNA/CLAUDE.md`** — Why key design decisions were made, what constraints exist
3. **`1-scRNA/reports/analysis_report.md`** — Full methods, results, limitations
4. **`1-scRNA/workflows/run_pipeline.py`** — How stages are orchestrated
5. **`1-scRNA/configs/analysis.yaml`** — What parameters exist
6. **`1-scRNA/src/icb_scrna/`** — Individual stage implementations
7. **`1-scRNA/tests/`** — Regression tests showing expected behavior

## Key patterns demonstrated

### 1. Configuration as code
```python
from config_base import load_config

cfg = load_config('configs/analysis.yaml')
cfg['signature']['n_genes']        # Access parameters
cfg.path('results')                # Resolved path object
```
Everything lives in YAML; no hardcoded values in source.

### 2. Stage-based pipeline
```python
from pipeline_driver import PipelineDriver, Stage

driver = PipelineDriver(
    stages=[Stage('ingest', '...', func=ingest.main), ...],
    outputs={'ingest': ['data.h5ad'], ...}
)
driver.run(config, force=args.force, stage=args.stage)
```
Each stage is independent; skip logic handles re-runs.

### 3. Test fixtures without external data
```python
from test_fixtures import MockData, Validator

cfg = MockData.config_yaml(tmp_path / 'test.yml')
X = MockData.sparse_matrix(nrows=100, ncols=50)
Validator.assert_valid_matrix(X, min_rows=50)
```
Tests build mock data; they never depend on `data/raw/`.

### 4. Data validation and logging
```python
from data_validation import DataManifest, ParseLog

manifest = DataManifest('data/metadata/manifest.json')
manifest.add('data/raw/gse.gz', 'TPM matrix', compute=True)
manifest.validate('data/raw/gse.gz')

logger = ParseLog('data/metadata/parse_log.json')
logger.add('ingest', 'genes_dropped', '342', 'detected in < 10 cells')
```
Record what happened during parsing so the next person knows.

### 5. Documentation of decisions
See `1-scRNA/CLAUDE.md`:
- Why biopsy is the unit of inference (not patient)
- Why DESeq2 wasn't used (no raw counts)
- What Windows sandbox constraints exist (and how to work around them)

And `1-scRNA/README.md`:
- What the analysis does and why
- How to run it
- Main findings and limitations

## Design principles

**Separation of concerns:**
- Analysis logic → `src/`
- Orchestration → `workflows/`
- Parameters → `configs/` (YAML)
- Tests → `tests/` (no external data needed)
- Tools → `0-reproducible-skills/`

**Single source of truth:**
- Parameters: YAML config files (not scattered through code)
- Environment: `environment.yml` spec + `environment.lock.txt` record
- Decisions: CLAUDE.md (why, not how)
- Data: Manifests with checksums (what was downloaded and verified)

**Reproducibility:**
- Fixed random seed in config
- Pinned environment (both spec and lockfile)
- All parameters tracked
- Data checksums recorded
- Full re-runs from clean data reproduce prior outputs (verified)
- Test suite catches regressions

**Fail early, document well:**
- Tests catch errors during development
- DataManifest catches corrupted downloads
- ParseLog documents what was excluded
- CLAUDE.md captures constraints for the next person

## Common tasks

### Running the scRNA analysis
```bash
cd 1-scRNA
conda env create -f environment/environment.yml
conda activate icb-scrna
python workflows/run_pipeline.py --config configs/analysis.yaml
pytest tests/
```

### Re-running one stage
```bash
python workflows/run_pipeline.py --config configs/analysis.yaml --stage abundance
```

### Re-running everything from scratch
```bash
python workflows/run_pipeline.py --config configs/analysis.yaml --force
```

### Starting a new analysis project
1. Create a new directory: `mkdir 2-newproject`
2. Copy the structure:
   ```bash
   mkdir -p 2-newproject/{configs,src,workflows,tests,data,results,figures,reports,environment}
   cp 0-reproducible-skills/config/config_base.py 2-newproject/src/
   cp 0-reproducible-skills/environment/environment_template.yml 2-newproject/environment/environment.yml
   cp 0-reproducible-skills/documentation/{CLAUDE,README}_TEMPLATE.md 2-newproject/
   mv 2-newproject/CLAUDE_TEMPLATE.md 2-newproject/CLAUDE.md
   mv 2-newproject/README_TEMPLATE.md 2-newproject/README.md
   ```
3. See `0-reproducible-skills/README.md` for what modules to import and how to use them.

### Adding a new dependency
1. Edit `environment/environment.yml` (the spec)
2. Recreate the environment:
   ```bash
   conda env create --force -f environment/environment.yml
   ```
3. Validate with tests:
   ```bash
   pytest
   ```
4. Commit the new `environment.lock.txt`:
   ```bash
   conda list --export > environment/environment.lock.txt
   git add environment/environment.lock.txt
   ```

---
