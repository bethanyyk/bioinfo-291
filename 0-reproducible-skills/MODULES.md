# Reproducible Skills Modules

This document describes what each subfolder contains and when to use it.

## `config/`

**What it does:** YAML-based configuration system with automatic path resolution.

**Files:**
- `config_base.py` — Config loader and Config class

**Use it when:** You need to externalize parameters from code. Loads a YAML file, gives you dict-style access and resolved Path objects for data directories.

**Example:**
```python
from config_base import load_config

cfg = load_config('configs/analysis.yaml')
n_genes = cfg['signature']['n_genes']              # dict-style access
results_dir = cfg.path('results')                   # auto-creates, returns Path object
stage_config = cfg.stage('abundance')               # get stage-specific params
```

**Key design:** Parameters live in YAML, not scattered through code. Makes experiments reproducible and easy to compare (just diff the config files).

---

## `pipeline/`

**What it does:** Orchestrates analysis stages with skip logic and CLI flags.

**Files:**
- `pipeline_driver.py` — PipelineDriver and Stage classes

**Use it when:** You have a multi-stage analysis (ingest → preprocess → annotate → analyze). The driver handles:
- **Skip logic** — If stage outputs exist, skip the stage (re-run only if needed)
- **`--force` flag** — Re-run everything even if outputs exist
- **`--stage NAME` flag** — Run only one stage
- **Logging** — Reports what ran and why

**Example:**
```python
from pipeline_driver import PipelineDriver, Stage

stages = [
    Stage('ingest', 'Load raw data', func=ingest.main),
    Stage('preprocess', 'QC and normalization', func=preprocess.main),
    Stage('annotate', 'Cell type assignment', func=annotate.main),
]

driver = PipelineDriver(
    stages=stages,
    outputs={
        'ingest': ['data.h5ad'],
        'preprocess': ['data.h5ad'],
        'annotate': ['data.h5ad'],
    }
)

driver.run(config, force=args.force, stage=args.stage)
```

**Key design:** Separates orchestration (which stages to run, in what order) from implementation (what each stage does). Makes pipelines debuggable and re-runnable.

---

## `testing/`

**What it does:** Test fixtures and validators for reproducible tests.

**Files:**
- `test_fixtures.py` — MockData, Validator, and Compare classes

**Use it when:** Writing tests that shouldn't depend on external data. Provides:
- **MockData** — Generate random matrices, DataFrames, and config files on the fly
- **Validator** — Check shape, dtype, missing values, numerical closeness
- **Compare** — Detailed failure messages when sets or DataFrames don't match

**Example:**
```python
from test_fixtures import MockData, Validator, Compare

def test_parsing():
    cfg = MockData.config_yaml(tmp_path / 'test.yml', seed=42)
    X = MockData.sparse_matrix(nrows=100, ncols=50, density=0.1)
    
    result = parse_function(X)
    
    Validator.assert_valid_matrix(result, min_rows=50, max_cols=5000)
    Validator.assert_no_missing(result)
    Compare.dataframe_equal(result, expected_df)
```

**Key design:** Tests are fast, deterministic, and run without downloading anything. Use these fixtures instead of committing large test data files.

---

## `utilities/`

**What it does:** Data validation, checksums, and parsing transparency.

**Files:**
- `data_validation.py` — DataManifest, CSVValidator, and ParseLog classes

**Use it when:**
- You download external data and want to verify it wasn't corrupted
- You parse messy files and want to document what was excluded
- You need to detect encoding issues early

**Three classes:**

### `DataManifest`
Record and verify file integrity via SHA-256.
```python
from data_validation import DataManifest

manifest = DataManifest('data/metadata/manifest.json')
manifest.add('data/raw/gse_expression.gz', description='TPM matrix', compute=True)
manifest.save()

# Later, verify the file is intact
manifest.validate('data/raw/gse_expression.gz')  # Raises if checksum mismatch
```

### `CSVValidator`
Validate encoding, shape, and column names before parsing.
```python
from data_validation import CSVValidator

# Detect encoding
enc = CSVValidator.check_encoding('file.csv', expected=['utf-8', 'latin-1'])

# Check shape
info = CSVValidator.check_shape(
    'file.csv',
    expected_columns=['gene', 'count'],
    expected_rows=50000,
)
```

### `ParseLog`
Document parsing decisions (what was excluded, why).
```python
from data_validation import ParseLog

logger = ParseLog('data/metadata/parse_log.json')
logger.add('ingest', 'genes_dropped', '342', 'detected in < 10 cells')
logger.add('preprocess', 'cells_filtered', '1204 of 16291 (7.4%)', 'QC thresholds')
logger.save()
```

**Key design:** Validate early, log decisions, make reproducibility auditable. The next person reading your code knows exactly what happened and why.

---

## `environment/`

**What it does:** Templates and guidance for conda environment management.

**Files:**
- `environment_template.yml` — Starter conda environment spec
- `environment_notes.md` — Pinning strategy and reproducibility guide

**Use it when:** Setting up a new analysis project.

**Workflow:**
1. Copy `environment_template.yml` to your project as `environment/environment.yml`
2. Add your domain-specific packages (scanpy, torch, whatever you need)
3. Run `conda env create -f environment/environment.yml`
4. Work and validate
5. After success, commit `environment.lock.txt` (exact resolved versions)

**Key design:** The `.yml` is what you write (human-readable, minimal). The `.lock.txt` is what conda resolves (exact versions, for replay). You commit both.

**Pinning strategy** (from `environment_notes.md`):
- **Never pin minor versions** unless load-bearing (incompatible, broken, or fails to load)
- **Pin load-bearing dependencies** with justification (e.g., `harmonypy==0.0.9` because 0.0.10+ imports torch, which fails on Windows sandboxes)
- **Record exact versions** after successful runs in `environment.lock.txt`

---

## `documentation/`

**What it does:** Templates for documenting your project.

**Files:**
- `CLAUDE_TEMPLATE.md` — Template for design decisions, data contracts, constraints
- `README_TEMPLATE.md` — Template for findings, how to run it, limitations

**Use it when:** Starting a new project.

**CLAUDE_TEMPLATE.md** documents the "why":
- Design decisions (why algorithm A over B?)
- Data contracts (format, encoding, special cases)
- Discovered constraints and workarounds (platform limitations, gotchas)

**README_TEMPLATE.md** documents the "what":
- What the analysis does and why
- How to run it
- Main findings and limitations
- Methods in brief

**Key design:** Future you (and your team) reads CLAUDE.md to understand decisions and constraints, README.md to understand findings and reproduction.

---

## Quick Reference

| Folder | Purpose | Copy into your project? |
|--------|---------|--------|
| `config/` | Externalize parameters to YAML | Yes — copy `config_base.py` to `src/` |
| `pipeline/` | Orchestrate multi-stage analysis | Yes — copy `pipeline_driver.py` to `src/` or `workflows/` |
| `testing/` | Test fixtures without external data | Yes — copy `test_fixtures.py` to `tests/` |
| `utilities/` | Data validation and checksums | Yes — copy `data_validation.py` to `src/` |
| `environment/` | Environment spec and pinning guide | Mostly reference — copy `.yml` template to `environment/` |
| `documentation/` | Design and findings templates | Copy templates to project root, fill them in |

---

## All together: How a project uses these

```
your-project/
├── configs/analysis.yaml           # Created from config template
├── src/
│   ├── config_base.py              # Copied from config/
│   ├── pipeline_driver.py           # Copied from pipeline/
│   ├── data_validation.py           # Copied from utilities/
│   ├── ingest.py                   # Your code
│   ├── preprocess.py               # Your code
│   └── ...
├── workflows/run_pipeline.py       # Uses PipelineDriver
├── tests/
│   ├── test_fixtures.py            # Copied from testing/
│   ├── test_ingest.py              # Your tests, using MockData
│   └── ...
├── environment/
│   ├── environment.yml             # Copied from environment/, then customized
│   └── environment.lock.txt        # Created after conda env create
├── CLAUDE.md                       # Filled in from documentation template
└── README.md                       # Filled in from documentation template
```

Each module is self-contained and meant to be copied into your project. You're not importing from `0-reproducible-skills/`; you're copying patterns.
