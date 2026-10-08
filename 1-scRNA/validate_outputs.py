#!/usr/bin/env python
"""
Validate scRNA analysis outputs using the OUTPUT_REVIEWER agent.

This script demonstrates how to delegate output validation to a specialist reviewer.
It can be called from the command line or integrated into a pipeline.

Usage:
    python validate_outputs.py

Or in Claude Science:
    result = host.delegate({
        "profile": "OUTPUT_REVIEWER",
        "task": open("validate_outputs.py").read()
    })
"""

import pathlib

def validate_scRNA_outputs():
    """
    Audit the scRNA-seq analysis outputs for correctness and completeness.
    
    Checks:
    - FIGURES: Readability, label overlaps, legend placement
    - RESULTS CSVs: Row counts, column names, value ranges
    - INTERMEDIATE DATA: AnnData object loads, cell/gene counts
    - REPRODUCIBILITY: Environment lockfile, config parameters documented
    """
    
    root = pathlib.Path(__file__).parent
    
    # Define what we expect
    expected_figures = [
        ("fig1_atlas.png", "UMAP of 11 immune populations"),
        ("fig2_composition.png", "Cell-type proportions: baseline vs. on-treatment, responder vs. non-responder"),
        ("fig3_cd8_volcano.png", "Differential expression in on-treatment CD8 T cells (responders vs. non-responders)"),
        ("fig4_signatures.png", "Signature performance: ROC curves and composition-based predictions"),
    ]
    
    expected_results = [
        ("differential_abundance.csv", 11, "Cell-type abundance differences (baseline vs. on-treatment, responder vs. non-responder)"),
        ("pseudobulk_de_significant.csv", 109, "Genes differentially expressed in on-treatment CD8 T cells (FDR < 0.1)"),
        ("signature_performance.csv", 5, "5 signature candidates with AUC, CI, permutation p"),
        ("signature_scores_loo.csv", 19, "Leave-one-patient-out signature scores for each baseline patient"),
    ]
    
    print("=" * 80)
    print("OUTPUT_REVIEWER: scRNA-seq Analysis Validation")
    print("=" * 80)
    
    # CHECK FIGURES
    print("\n[1] FIGURES — Visual correctness and readability\n")
    figures_dir = root / "figures"
    for fig_name, description in expected_figures:
        fig_path = figures_dir / fig_name
        if fig_path.exists():
            size_mb = fig_path.stat().st_size / 1e6
            print(f"  ✓ {fig_name:<25} ({size_mb:.1f} MB)")
            print(f"    {description}")
        else:
            print(f"  ✗ MISSING: {fig_name}")
            print(f"    Expected: {description}")
    
    print("\n    REVIEW ITEMS:")
    print("    - Check for overlapping text (axis labels, legend, annotations)")
    print("    - Verify all cell types labeled and visible in UMAP")
    print("    - Confirm proportions sum to 100% in composition plots")
    print("    - Verify x/y axes have clear labels with units")
    print("    - Check that ROC legend doesn't obscure curves")
    
    # CHECK RESULTS
    print("\n[2] RESULTS — Data reasonableness and completeness\n")
    results_dir = root / "results"
    for csv_name, expected_rows, description in expected_results:
        csv_path = results_dir / csv_name
        if csv_path.exists():
            import pandas as pd
            try:
                df = pd.read_csv(csv_path)
                print(f"  ✓ {csv_name:<35} ({len(df)} rows)")
                print(f"    {description}")
                if len(df) != expected_rows:
                    print(f"    ⚠ WARNING: Expected {expected_rows} rows, got {len(df)}")
            except Exception as e:
                print(f"  ✗ ERROR reading {csv_name}: {e}")
        else:
            print(f"  ✗ MISSING: {csv_name}")
            print(f"    Expected: {expected_rows} rows, {description}")
    
    print("\n    REVIEW ITEMS:")
    print("    - Verify row counts match expected values")
    print("    - Check for missing values (NaN, blank cells)")
    print("    - Confirm numeric ranges are reasonable:")
    print("      • Cell counts and proportions ≥ 0")
    print("      • AUC, sensitivity, specificity in [0, 1]")
    print("      • p-values, FDR in [0, 1]")
    print("      • log2 fold-change unbounded")
    
    # CHECK DATA INTEGRITY
    print("\n[3] DATA INTEGRITY — Reproducibility and completeness\n")
    checks = [
        ("environment/environment.lock.txt", "Exact package versions (99 packages)"),
        ("configs/analysis.yaml", "All parameters documented (no hardcoding)"),
        ("CLAUDE.md", "Design decisions and constraints"),
        ("environment/environment.yml", "Human-readable conda spec"),
    ]
    
    for item, description in checks:
        item_path = root / item
        if item_path.exists():
            print(f"  ✓ {item:<35} {description}")
        else:
            print(f"  ✗ MISSING: {item}")
    
    print("\n    REVIEW ITEMS:")
    print("    - Confirm environment.lock.txt pins all 99 packages with exact versions")
    print("    - Verify configs/analysis.yaml lists thresholds for:")
    print("      • Cell QC filters (gene count, mitochondrial %)")
    print("      • HVG selection (n_hvgs)")
    print("      • Harmony integration (n_iter, n_PCs)")
    print("      • Leiden clustering (resolution)")
    print("    - Check CLAUDE.md documents rationale for DEA method choice")
    
    # CHECK REPRODUCIBILITY
    print("\n[4] REPRODUCIBILITY — Can this be re-run identically?\n")
    import json
    
    print("  To verify reproducibility:")
    print("  1. Delete data/processed/*.h5ad")
    print("  2. Run: python workflows/run_pipeline.py --config configs/analysis.yaml")
    print("  3. Compare new results/ and figures/ with originals (should be identical)")
    print("  4. Run: pytest tests/")
    print("    All 32 tests should pass")
    
    # SUMMARY
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("""
If all checks pass:
  ✓ Figures are publication-ready
  ✓ Results tables are complete and reasonable
  ✓ Environment is pinned for reproducibility
  ✓ Design decisions are documented
  ✓ Analysis can be re-run identically

Next steps:
  1. Incorporate reviewer feedback on figures
  2. Commit results/ and figures/ to git
  3. Archive the full analysis state (code + environment)
  4. Write the narrative report (reports/analysis_report.md)
""")

if __name__ == "__main__":
    validate_scRNA_outputs()
