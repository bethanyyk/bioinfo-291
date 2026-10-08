# Project Implementation Notes (CLAUDE.md)

This file documents the **how and why** of the project's implementation—conventions, design decisions, and constraints—so future sessions can maintain consistency without guessing.

Use this for:
- Directory structure and file organization
- Data contracts (what format, encoding, assumptions about each file)
- Design decisions that aren't obvious from code (why algorithm A over B, why this threshold)
- Constraints discovered during development (platform limitations, memory budgets, sandbox gotchas)
- Known workarounds or one-off steps required to run the project

Do **not** use this for:
- General scientific background (put that in the README or report)
- Hardcoded values (put them in config files)
- Step-by-step instructions (that's for README or comments in code)

## Directory Contract

### `configs/`
YAML configuration files. Every parameter that varies between runs lives here, never in code.

**Expected files:**
- `analysis.yaml` — stage parameters, paths, thresholds, random seed

### `src/`
Package source code, organized by analysis stage.

**Naming convention:** one module per stage, e.g. `ingest.py`, `preprocess.py`, `abundance.py`.

Each module should:
- Have a single entry-point function matching the stage name
- Take a `config` parameter as its first argument
- Log progress and errors to the module logger
- Return or write results; don't rely on side effects

### `tests/`
Test suite. Run with `pytest tests/` from the project root.

Conventions:
- One test file per source module, e.g., `test_ingest.py`
- Tests use fixtures from `conftest.py` (shared setup)
- Tests build mock data, never depend on `data/raw/`

### `data/`
Local data storage.

- `data/raw/` — downloaded or external datasets (gitignored, regenerable)
- `data/processed/` — intermediate HDF5, parquet, or pickle files (gitignored, regenerable)
- `data/metadata/` — small metadata files and manifests (committed)

### `results/`
Analysis outputs: CSV tables, statistical summaries. Committed to the repo.

### `figures/`
Final publication-grade figures as PNG or PDF. Committed to the repo.

### `reports/`
Narrative reports (Markdown, LaTeX, or HTML) describing findings, methods, and limitations.

## Design Decisions

### [Your decision name]
**What:** [Brief description of the choice]

**Why:** [Rationale. Why not the alternative?]

**Tradeoff:** [What you gave up or constrained by choosing this.]

**Example from bioinfo-291:**

### Biopsy as the unit of inference, not cell or patient
**What:** All statistical tests use biopsy (tissue sample) as the unit, grouped by patient for robust errors.

**Why:** Response is annotated per biopsy. Three patients have discordant labels across lesions, so patient-level averaging would be wrong.

**Tradeoff:** Smaller effective sample size (48 biopsies vs. 32 patients). But correctness > statistical power.

## Constraints & Workarounds

### Windows sandbox: pytest rootdir discovery fails
**Symptom:** `PermissionError: [WinError 5] Access is denied: 'C:\'`

**Root cause:** pytest's common-ancestor discovery walks to the drive root, which the sandbox blocks.

**Failed workarounds:**
- Passing `--rootdir` and `-c pytest.ini` does NOT help; the error occurs before config is parsed.

**Working solution:**
Copy tests to a directory outside the drive-root ancestry and run with `PYTHONPATH`:
```bash
cp tests/*.py /tmp/ptests/
PYTHONPATH=/path/to/project/src pytest /tmp/ptests/ -q
```

### [Your constraint]
**Symptom:** [What went wrong]

**Root cause:** [Why it happens]

**Workaround:** [How to fix it]

## Data Format Notes

### Expression matrix (`.h5ad` or `.tsv.gz`)
- Format: [HDF5/tab-separated/etc.]
- Encoding: [UTF-8, latin-1, etc.]
- First row: [headers / sample names / other]
- Expected columns: [gene names, structure]
- Dtype: [float32, int32, etc.]
- Special cases: [missing values as NaN / -1 / omitted? Any cells/genes to exclude?]

### Annotation file (`.csv`)
- Format: [CSV, TSV, other]
- First [N] lines: [preamble, describe if present]
- Encoding: [UTF-8, latin-1, etc.]
- Index column: [sample ID / cell ID / name]
- Expected columns: [list]
- Special cases: [any duplicate rows? Multi-value cells?]

### Output tables (`.csv`, `.parquet`, etc.)
- Index: [what each row represents: gene, cell type, contrast, etc.]
- Columns: [describe all columns]
- Sorting: [are rows sorted? By what?]
- Missing values: [NaN / empty string / omitted?]

## Reproducibility Checklist

- [ ] All parameters live in `configs/` (not hardcoded in source)
- [ ] Random seed is fixed in config and used in all stochastic steps
- [ ] Environment spec (`environment.yml`) and lockfile (`environment.lock.txt`) are committed
- [ ] Pinned versions have load-bearing justification in the environment.yml comments
- [ ] Test suite passes (`pytest tests/`)
- [ ] A full re-run from scratch reproduces prior outputs (spot-check key results)
- [ ] Data manifest and checksums are recorded (if data is external)
- [ ] README describes how to run the pipeline end-to-end
- [ ] CLAUDE.md documents why key design decisions were made
- [ ] Methods section in the report names exact parameter values and thresholds
