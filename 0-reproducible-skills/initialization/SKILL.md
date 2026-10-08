---
name: initialize-project-structure
description: "Initialize a new bioanalysis project with standard directory structure. Enforces exactly 9 top-level directories: data (with raw/, metadata/, processed/), src, workflows, configs, results, figures, reports, tests, environment."
---

# Initialize Project Structure

Initialize a reproducible bioanalysis project with the standard directory layout used across bioinfo-291.

## Standard Structure

Every project must have **exactly** these 9 top-level directories (no more, no fewer):

```
your-project/
├── data/              # Raw data, metadata, processed outputs
│   ├── raw/           # Downloaded or original data (do not modify)
│   ├── metadata/      # Manifests, annotations, checksums
│   └── processed/     # Intermediate & final output files
├── src/               # Python/R source code (modules, utilities)
├── workflows/         # Pipeline orchestration scripts
├── configs/           # YAML configuration files (parameters, thresholds)
├── results/           # Analysis tables, metrics, statistics
├── figures/           # Publication-quality figures
├── reports/           # Narrative reports, extended methods
├── tests/             # Unit tests, regression tests, fixtures
└── environment/       # Conda spec (environment.yml) and lock (environment.lock.txt)
```

**No other top-level directories are allowed.** If your project needs additional structure, create it *within* one of these directories (e.g., `data/downloads/`, `src/utils/`, `results/supplementary/`).

## Usage

### Initialize a new project

```python
from initialize_project_structure import init_project_structure

# Create the standard directory structure
init_project_structure("./my-analysis")
# → Prints: ✓ my-analysis initialized
#           ✓ Standard structure created
```

### Validate an existing project

```python
from initialize_project_structure import validate_project_structure

# Check that a project has the correct structure
issues = validate_project_structure("./my-analysis")

if issues:
    print("Structure issues found:")
    for issue in issues:
        print(f"  - {issue}")
else:
    print("✓ Project structure is valid")
```

## Why this structure?

1. **Separation of concerns** — Code, configuration, data, and results each have their place
2. **Clarity** — Anyone picking up the project knows where to find things
3. **Reproducibility** — Parameters live in YAML; intermediate outputs are gitignored; code is version-controlled
4. **Consistency** — All projects follow the same layout, making templates reusable
5. **Strictness** — Only 9 top-level directories; no project-specific folders at the root level

See the `1-scRNA/` example in `../` for a complete worked example.

## Examples

### From the command line

```bash
# Initialize a new spatial transcriptomics project
python -c "from initialize_project_structure import init_project_structure; init_project_structure('spatial-txn-analysis')"

# Validate an existing project
python -c "from initialize_project_structure import validate_project_structure; validate_project_structure('spatial-txn-analysis')"
```

### In a Python script

```python
from initialize_project_structure import init_project_structure, validate_project_structure

# Set up a new analysis
project_path = "my-analysis-2026"
init_project_structure(project_path)

# Later, verify the structure before committing
issues = validate_project_structure(project_path)
assert not issues, f"Structure problems: {issues}"
```

### In a notebook or REPL

```python
# After loading this skill, the functions are available in the kernel:
init_project_structure("./new-project")
issues = validate_project_structure("./new-project")
print("✓ Project is valid" if not issues else f"Issues: {issues}")
```

## What it does

### `init_project_structure(path)`

- Creates the specified directory if it doesn't exist
- Creates all 9 required top-level directories
- Creates the 3 required subdirectories within `data/` (raw, metadata, processed)
- Skips directories that already exist (idempotent — safe to run multiple times)
- Prints progress to stdout
- Raises `OSError` if creation fails

### `validate_project_structure(path)`

- Checks that the directory exists and is a directory
- Checks that all 9 required top-level directories exist
- Checks that `data/` contains exactly the 3 required subdirectories (no more, no fewer)
- Checks that no unexpected top-level directories exist
- Returns a list of issue strings (empty if valid)
- Example output:
  ```
  [
    "Missing directory: workflows/",
    "Unexpected directory: scratch/",
    "data/ missing subdirectory: processed/",
    "data/ unexpected subdirectory: downloads/"
  ]
  ```

## Integration with reproducible-skills

This module is part of the reproducible-skills toolkit. To use it in a new project:

1. Copy the initialization module to your project:
   ```bash
   cp -r 0-reproducible-skills/initialization/ your-project/
   ```

2. Call from your project setup script:
   ```python
   from initialization.initialize_project_structure import init_project_structure
   init_project_structure(".")
   ```

Or use it standalone to scaffold new projects:
```bash
cd ~/projects
python -c "import sys; sys.path.insert(0, '../bioinfo-291/0-reproducible-skills'); from initialization import init_project_structure; init_project_structure('new-spatial-analysis')"
```

## Testing

The skill includes validation, so use it to check any project:
```python
issues = validate_project_structure("../some-other-project")
if issues:
    print("Fix these issues:")
    for issue in issues:
        print(f"  {issue}")
```

---

See `../README.md` for the full reproducible-skills toolkit.
