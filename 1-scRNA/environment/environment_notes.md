# Environment Management for Reproducible Workflows

## Pinning Strategy

### Never pin minor versions (X.Y) unless load-bearing.

A "load-bearing" pin is required when:
- A later release is incompatible with your code (e.g., an API break).
- A later release has a critical bug affecting your results.
- A dependency fails to load in your sandboxed environment (e.g., torch DLLs on Windows).

If you can run with `numpy` and `scipy` as-is, don't write `numpy==1.24.0`. Pin only `harmonypy==0.0.9` if 0.0.10+ breaks your code.

### Record exact versions after a successful run.

After your pipeline runs successfully and you've validated the results:
1. Export the realized environment:
   ```bash
   conda env export --name <env-name> > environment.lock.txt
   ```
   Or for pip:
   ```bash
   pip freeze > environment.lock.txt
   ```
2. Commit `environment.lock.txt` alongside `environment.yml`.
3. Document which dependencies are pinned and why (see the example in `bioinfo-291/1-scRNA/environment/`).

This enables:
- **Exact replay** — future runs use the same package versions.
- **Debugging** — if a result differs later, you can see what changed.
- **Transparency** — new team members see what was resolved.

## Load-Bearing Pins: Windows Sandbox Example

From bioinfo-291:

```yaml
# harmonypy 0.0.9 — releases from 0.0.10 import torch at module load
# torch's CUDA/CPU DLLs fail to load on Windows sandboxes
# The integration itself needs no GPU and no deep learning; 0.0.9 works.
- harmonypy==0.0.9
```

When documenting a load-bearing pin:
- Name the exact issue (incompatibility, import failure, API break).
- Confirm the code works with the pinned version.
- Add a regression test that would catch a future pin lift.

## Conda vs. pip

- Use conda for compiled C/Fortran dependencies (scipy, statsmodels, scikit-learn, torch).
- Use pip in the conda environment for pure-Python packages or when conda isn't available.
- For a conda env, declare both in `dependencies` and list `pip` as a dependency, then use a `pip:` subsection.

Example:
```yaml
dependencies:
  - python=3.12
  - numpy
  - scipy
  - scikit-learn
  - pip
  - pip:
      - harmonypy==0.0.9
      - some-pypi-only-package==1.2.3
```

## Reproducing Later

On a new machine or future session:

1. **Exact replay:**
   ```bash
   conda create --name analysis --file environment.lock.txt
   ```

2. **Latest compatible versions:**
   ```bash
   conda env create -f environment.yml
   conda list --export > environment.lock.txt  # then commit the new lock
   ```

3. **Adding a package:**
   Edit `environment.yml`, re-run `conda env create`, and commit the new lock.

## Testing the Environment

After creating the environment, run the test suite to confirm all dependencies resolve correctly:

```bash
conda activate <env-name>
pytest -v
```

If a test fails due to a missing or broken dependency, it caught the problem early.
