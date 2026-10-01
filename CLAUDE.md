# bioinfo-291 — repository conventions

Single-cell analysis of melanoma tumour biopsies taken before and during immune
checkpoint blockade (GEO **GSE120575**, Sade-Feldman et al., *Cell* 2018).

## Directory contract

Every directory has one job. Nothing is written outside this contract.

| Path | Contents | In git? |
|---|---|---|
| `configs/` | YAML parameter files. No parameter is hardcoded in `src/`. | yes |
| `data/raw/` | Files downloaded verbatim from GEO. Never edited, never written to by analysis code. | no (gitignored) |
| `data/processed/` | Derived `.h5ad` objects. Regenerable from `data/raw/` + code. | no (gitignored) |
| `data/metadata/` | Small per-sample / per-cell annotation tables and download checksums. | yes |
| `src/icb_scrna/` | Importable Python package. One module per pipeline stage. | yes |
| `workflows/` | Stage driver that runs the pipeline end to end from a config. | yes |
| `results/` | Analysis outputs as CSV/TSV. One file per question asked. | yes |
| `figures/` | Publication-intent figures, PNG at 300 dpi. | yes |
| `reports/` | Human-readable write-ups referencing `results/` and `figures/`. | yes |
| `tests/` | pytest suite. Runs without network access and without `data/raw/`. | yes |
| `environment/` | Conda environment specification and a realized lockfile. | yes |

## Conventions

- **Config-driven.** Stages read `configs/analysis.yaml`. Changing a threshold means
  editing the YAML, not the code. `src/icb_scrna/config.py` is the only loader.
- **Stages are idempotent.** Each stage writes a declared output path and is skipped
  by the driver when that path already exists, unless `--force` is passed.
- **Raw data is immutable.** Code reads `data/raw/` and never writes there.
- **Random seeds are fixed** in the config (`seed: 0`) and passed explicitly to every
  stochastic step (PCA, Harmony, Leiden, UMAP, cross-validation splits).
- **Sample is the unit of inference**, not the cell. Any test comparing patient groups
  aggregates to the biopsy or patient level first; per-cell tests across patients
  inflate significance and are not used for clinical contrasts.
- **Figures** are built with an explicit `fig, ax = plt.subplots()` handle and saved
  with `fig.savefig(...)`. No global `plt.savefig`.
- **Tables** carry their statistics in full: effect size, raw p, and the multiple-testing
  adjusted value. No table reports a p-value without the correction method named in its
  column header.

## Memory budget

This analysis is developed on a machine with roughly 2–3 GB of usable RAM. The GEO
expression matrix is a ~55k gene x 16k cell dense text file and **must not** be read
with a plain `pandas.read_csv`. `src/icb_scrna/ingest.py` streams it in row chunks and
converts to sparse `float32` on the fly. Keep that property if you touch the loader.

## Running

```bash
conda env create -f environment/environment.yml
conda activate icb-scrna
python workflows/run_pipeline.py --config configs/analysis.yaml
```

Individual stages: `python workflows/run_pipeline.py --config configs/analysis.yaml --stage atlas`

## Tests

```bash
pytest
```

`pytest.ini` pins `testpaths` and `pythonpath`, so no editable install is needed.
The suite builds its own fixtures and runs without network access and without
`data/raw/`. 32 tests, all passing as of the last run.

On a restricted filesystem where `stat()` on the drive root is denied, pytest's
rootdir discovery fails before collection with
`PermissionError: [WinError 5] Access is denied: 'C:\'`. Passing `--rootdir`
and `-c pytest.ini` does **not** help — it was tried and fails identically,
because the error occurs in pytest's common-ancestor resolution before the
config is consulted. The workaround that does work is to run the tests from a
directory that does not share a drive-root ancestry with the interpreter, with
the package on `PYTHONPATH`:

```bash
cp tests/*.py /some/other/dir/ptests/
PYTHONPATH=/path/to/bioinfo-291/src pytest /some/other/dir/ptests -q
```

On an unrestricted filesystem a plain `pytest` from the repository root works.

## Provenance

Raw sequencing reads for GSE120575 are controlled-access (dbGaP phs001680.v1.p1).
This project uses only the public processed TPM matrix and its sample annotation,
both downloaded from the GEO series supplementary files.
