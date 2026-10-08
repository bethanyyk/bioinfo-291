# bioinfo-291 — Reproducible Analysis Framework

This repository contains a **complete, reproducible bioanalysis pipeline** and a **reusable toolkit** for building similar pipelines in the future.

## Two things live here

### `0-reproducible-skills/` — Reusable across all projects
**Configuration, orchestration, testing, and documentation patterns** that apply to any analysis pipeline. Start here when building a new project; these tools eliminate boilerplate and ensure consistency.

See [`0-reproducible-skills/README.md`](0-reproducible-skills/README.md) for what's available:
- **Config system** — YAML → Python objects with path resolution
- **Pipeline driver** — Stage orchestration with skip logic and `--force`/`--stage` flags
- **Test fixtures** — Mock data, validators, comparison utilities
- **Documentation templates** — CLAUDE.md (decisions), README.md (findings)
- **Utility modules** — Data validation, checksums, parse logging

### `1-scRNA/` — A complete worked example
**A full analysis of melanoma tumour biopsies** asking which immune populations expand with checkpoint blockade and whether any baseline signature predicts response. This demonstrates:
- How to organize a multi-stage analysis pipeline
- How to use the reusable tools above
- How to document design decisions and constraints
- How to validate data and reproducibility

---

## Quick start: Using the reusable tools in your own project

1. **Copy the templates into your project:**
   ```bash
   cp 0-reproducible-skills/config/config_base.py       <your-project>/src/
   cp 0-reproducible-skills/environment/environment_template.yml  <your-project>/environment/
   cp 0-reproducible-skills/documentation/CLAUDE_TEMPLATE.md  <your-project>/CLAUDE.md
   cp 0-reproducible-skills/documentation/README_TEMPLATE.md  <your-project>/README.md
   ```

2. **Import the utilities you need:**
   ```python
   from config_base import load_config
   from pipeline_driver import PipelineDriver, Stage
   from test_fixtures import MockData, Validator
   from data_validation import DataManifest, ParseLog
   ```

3. **Organize your project like `1-scRNA/`:**
   ```
   your-project/
   ├── configs/analysis.yaml      # Parameters (YAML, not hardcoded)
   ├── src/                       # Analysis code
   ├── workflows/run_pipeline.py  # Stage orchestration
   ├── tests/                     # Tests with fixtures (no external data)
   ├── data/                      # Local data (raw/, processed/, metadata/)
   ├── results/                   # Analysis tables
   ├── figures/                   # Publication-grade figures
   ├── reports/                   # Narrative reports
   ├── environment/               # Conda spec + lockfile
   ├── CLAUDE.md                  # Design decisions & constraints
   └── README.md                  # How to run it, findings
   ```

See [`0-reproducible-skills/README.md`](0-reproducible-skills/README.md) for detailed examples.

---

## The scRNA example in detail (`1-scRNA/`)

A complete analysis of **melanoma tumour biopsies** before and during immune checkpoint blockade, asking:
- Which immune populations expand or contract with treatment?
- What marks them?
- Does any baseline signature predict response?

**Data:** GEO [GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575) (Sade-Feldman et al., *Cell* 2018). 16,291 CD45+ cells from 48 biopsies of 32 patients.

**Key findings:**
- No cell type changes abundance with treatment (FDR < 0.05)
- Baseline composition differs between responders and non-responders (7 of 11 populations, FDR < 0.05)
- On-treatment CD8 T cells carry an interferon program in non-responders (109 genes, FDR < 0.1)
- No baseline signature predicts response beyond chance (best AUC 0.73, 95% CI 0.48–0.96)

### Running the scRNA analysis
```bash
cd 1-scRNA
conda env create -f environment/environment.yml
conda activate icb-scrna
python workflows/run_pipeline.py --config configs/analysis.yaml
pytest tests/
```

**Runtime:** ~20 minutes on 8 cores, ~1.5 GB peak RAM

### Validating outputs with OUTPUT_REVIEWER

After the pipeline completes, use the **OUTPUT_REVIEWER** agent to audit figures and results:

```python
# In Claude Science
result = host.delegate({
    "profile": "OUTPUT_REVIEWER",
    "task": "Review 1-scRNA/figures/ and 1-scRNA/results/ for visual correctness, data reasonableness, and completeness."
})
```

The OUTPUT_REVIEWER checks for:
- **Overlapping text, poor contrast, unreadable labels** in figures
- **Missing files or incorrect row counts** in results
- **Values that fall outside expected ranges** (negative counts, AUC > 1, etc.)

This offloads visual QC and data sanity checks to a specialist, freeing you to focus on interpretation.

### Project layout
```
1-scRNA/
├── configs/analysis.yaml         # All parameters
├── src/icb_scrna/                # Analysis modules (ingest, preprocess, annotate, etc.)
├── workflows/run_pipeline.py     # Stage driver
├── tests/                        # 32 regression tests (all passing)
├── data/
│   ├── raw/                      # Downloaded GEO files
│   ├── processed/                # Intermediate HDF5 files
│   └── metadata/                 # Manifests and annotations
├── results/                      # Analysis tables (committed)
├── figures/                      # Final figures (committed)
├── reports/analysis_report.md    # Full methods, findings, limitations
├── environment/                  # Conda spec + lockfile (106 packages)
├── CLAUDE.md                     # Implementation notes & design decisions
└── README.md                     # How to run it, findings summary
```

**Key decisions documented in CLAUDE.md:**
- Biopsy (not patient) as the unit of inference
- Why DESeq2 wasn't used (no raw counts; only TPM available)
- Why binomial GLM with patient-clustered errors (propeller is R-only)
- Windows sandbox constraints (harmonypy pin, pytest rootdir workaround)

---

## Design principles

**Separation of concerns:**
- Code logic → `src/`
- Orchestration → `workflows/`
- Parameters → `configs/` (YAML, not hardcoded)
- Tests → `tests/` (no external data)

**Single source of truth:**
- Parameters: YAML config files
- Environment: `environment.yml` spec + `environment.lock.txt` record
- Decisions: CLAUDE.md

**Fail early, document well:**
- Tests catch errors during development
- DataManifest verifies downloads
- ParseLog documents what was excluded
- CLAUDE.md captures constraints for the next person

**Reproducibility:**
- Fixed random seed
- Pinned environment (spec + lock)
- All parameters tracked
- Data checksums recorded
- Full re-runs from clean data reproduce prior outputs

---

## For instructors

This repository demonstrates:

1. **Reproducible research practices** — configuration, testing, documentation, environment pinning
2. **Software engineering patterns** — modularity, separation of concerns, stage orchestration
3. **Defensive programming** — validation, checksums, regression tests
4. **Collaboration** — reusable tools so students don't repeat boilerplate across projects

Each of these is **generalizable and demonstrated** in both the abstract tools (`0-reproducible-skills/`) and a concrete, working example (`1-scRNA/`).

---

**Last updated:** October 2026

See individual READMEs in `0-reproducible-skills/` and `1-scRNA/` for more detail.
