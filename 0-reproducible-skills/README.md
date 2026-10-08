# Reproducible Bioanalysis Skills

This folder contains reusable tools, patterns, and templates for building reproducible analysis pipelines. Use these across multiple projects to avoid repeating boilerplate and maintain consistency.

## What's here

### `config/`
**Configuration management system.**

- `config_base.py` — YAML-based config loader with path resolution
  - Automatically creates output directories
  - Dict-style and method-style access (`config['key']` or `config.path('results')`)
  - Per-stage config access (`config.stage('abundance')`)

**Use this:** In every project, to keep parameters separate from code.

### `pipeline/`
**Stage-based pipeline orchestration.**

- `pipeline_driver.py` — Orchestrates analysis stages with skip logic
  - Automatically skips stages if outputs exist (unless `--force` is passed)
  - Supports `--stage NAME` to run a single stage
  - Logs what ran and why

**Use this:** To separate data processing logic from orchestration.

Example:
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

### `testing/`
**Test fixtures and validators for reproducibility.**

- `test_fixtures.py` — Mock data generators, validators, and comparison utilities
  - `MockData.sparse_matrix()`, `MockData.dataframe()`, `MockData.config_yaml()` — generate test data without external dependencies
  - `Validator.assert_valid_matrix()`, `assert_no_missing()`, `assert_close()` — numerical validation helpers
  - `Compare.sets_equal()`, `Compare.dataframe_equal()` — detailed assertion failures

**Use this:** To write tests that don't depend on external data, and to get readable error messages on failure.

Example:
```python
from test_fixtures import MockData, Validator

# Test fixture: a mock config file
cfg_path = MockData.config_yaml(tmpdir / 'test.yml', seed=42)

# Validation: check shape and dtype
X = some_function()
Validator.assert_valid_matrix(X, min_rows=100, max_cols=5000, dtype=np.float32)
Validator.assert_no_missing(X)

# Comparison: detailed failure message
Compare.dataframe_equal(result, expected)  # Shows exactly what differs
```

### `environment/`
**Environment specification and pinning guide.**

- `environment_template.yml` — Template for conda environment files
- `environment_notes.md` — Guide to pinning strategy and reproducibility

**Use this:**
1. Copy `environment_template.yml` to your project as `environment/environment.yml`
2. Add your domain-specific packages
3. Follow the pinning guide: only pin versions that are load-bearing (incompatible, broken, or fail to load)
4. After a successful run, commit `environment.lock.txt` alongside the `.yml`

### `documentation/`
**Documentation templates.**

- `CLAUDE_TEMPLATE.md` — Template for CLAUDE.md (implementation notes, design decisions, constraints)
- `README_TEMPLATE.md` — Template for README.md (what the analysis does, how to run it, findings)

**Use this:**
1. Copy these templates to your project
2. Fill in sections specific to your analysis
3. Commit alongside your code

### `utilities/`
**Common utility modules.**

- `data_validation.py` — File manifest (checksums), CSV validation, parse logging
  - `DataManifest` — Record and verify file integrity
  - `CSVValidator` — Check encoding, shape, and column names
  - `ParseLog` — Log what was excluded during parsing (e.g., "342 genes with <10 cells dropped")

**Use this:** To validate external data, record what you decided during parsing, and catch errors early.

Example:
```python
from data_validation import DataManifest, ParseLog

# Create a manifest for external downloads
manifest = DataManifest('data/metadata/manifest.json')
manifest.add('data/raw/gse_expression.gz', description='TPM matrix', compute=True)
manifest.save()

# Later, verify the file wasn't corrupted
manifest.validate('data/raw/gse_expression.gz')

# Log parsing decisions
logger = ParseLog('data/metadata/parse_log.json')
logger.add('ingest', 'genes_dropped', '342', 'detected in < 10 cells')
logger.add('preprocess', 'cells_filtered', '1204 of 16291 (7.4%)', 'QC thresholds')
logger.save()
```

## How to use these in a new project

### 1. Start with templates
```bash
# Copy config, environment, and documentation templates
cp 0-reproducible-skills/config/config_base.py          1-myproject/src/
cp 0-reproducible-skills/environment/environment_template.yml  1-myproject/environment/environment.yml
cp 0-reproducible-skills/documentation/CLAUDE_TEMPLATE.md      1-myproject/CLAUDE.md
cp 0-reproducible-skills/documentation/README_TEMPLATE.md      1-myproject/README.md
```

### 2. Import reusable modules
```python
# In your analysis code
from reproducible_skills.config import Config, load_config
from reproducible_skills.pipeline import PipelineDriver, Stage
from reproducible_skills.testing import MockData, Validator
from reproducible_skills.utilities import DataManifest, ParseLog
```

### 3. Set up your pipeline structure
```bash
mkdir -p 1-myproject/{src,tests,configs,workflows,data,results,figures,reports,environment}
python 1-myproject/workflows/run_pipeline.py --config 1-myproject/configs/analysis.yaml
```

### 4. Write tests
Use `test_fixtures.MockData` to generate fixtures without external dependencies:
```python
def test_parsing():
    cfg = MockData.config_yaml(tmp_path / 'test.yml', seed=42)
    result = parse_file(cfg)
    Validator.assert_valid_matrix(result, min_rows=100)
```

### 5. Document as you go
- Fill in CLAUDE.md with design decisions and constraints as you discover them
- Update README.md with findings and runtime expectations
- Keep environment.lock.txt in sync after successful runs

## Design principles

**Separation of concerns:**
- Code logic lives in `src/`
- Orchestration lives in `workflows/`
- Parameters live in `configs/`
- Tests live in `tests/` (no external data needed)

**Single source of truth:**
- Parameters: YAML config files (not hardcoded)
- Environment: `environment.yml` spec + `environment.lock.txt` record
- Decisions: CLAUDE.md (why, not how)

**Fail early, document well:**
- Tests catch errors during development
- DataManifest catches corrupted downloads
- ParseLog documents what was excluded and why
- CLAUDE.md captures constraints and workarounds for the next person

**Reproducibility:**
- Fixed random seed in config
- Environment is pinned (both spec and lock)
- All parameters are tracked
- Data checksums are recorded
- A full re-run from clean data reproduces prior outputs

## Examples

See `1-scRNA/` for a complete, working example:
- How to organize a multi-stage pipeline
- How to use ConfigBase for parameters
- How to write tests with MockData and Validator
- How to document decisions in CLAUDE.md and README.md
- How to handle environment constraints (Windows sandbox, harmonypy pin)

## Contributing

When you solve a problem that isn't specific to your project—a new config pattern, a testing utility, a documentation template—consider extracting it here so the next project can reuse it.

Examples of what to extract:
- A new fixture that generates mock data for your domain
- A new validation check that you'll use again
- A documented constraint or workaround
- A new stage pattern or orchestration approach

---

**Last updated:** October 2026

This folder is meant to grow. Add to it as you solve problems across projects.
