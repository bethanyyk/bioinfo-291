"""
Data inspection utilities for visualization planning.

Use these before creating any visualization to understand your data:
- What are the units, range, distribution?
- What is the sample composition?
- Are there outliers, missingness, or other issues?

These utilities answer: "Can I visualize this variable responsibly?"
"""

import numpy as np
import pandas as pd
from typing import Optional, Union


def inspect_variable(
    data: Union[pd.Series, np.ndarray],
    name: Optional[str] = None,
    groups: Optional[Union[pd.Series, np.ndarray]] = None,
) -> dict:
    """
    Inspect a variable before visualization.
    
    Returns:
        dict with keys: name, dtype, missing, range, mean, median, sd, 
                       distribution_shape, unique_count, min, max, q25, q75
    """
    data = pd.Series(data) if isinstance(data, np.ndarray) else data
    
    stats = {
        'name': name or data.name or 'variable',
        'dtype': str(data.dtype),
        'n': len(data),
        'missing': data.isna().sum(),
        'pct_missing': (data.isna().sum() / len(data) * 100),
        'min': float(data.min()),
        'q25': float(data.quantile(0.25)),
        'median': float(data.median()),
        'mean': float(data.mean()),
        'q75': float(data.quantile(0.75)),
        'max': float(data.max()),
        'sd': float(data.std()),
        'range': float(data.max() - data.min()),
        'iqr': float(data.quantile(0.75) - data.quantile(0.25)),
        'skewness': float(data.skew()),
        'kurtosis': float(data.kurtosis()),
        'unique_count': data.nunique(),
    }
    
    print(f"\n{stats['name']} ({stats['dtype']})")
    print(f"  n={stats['n']}, missing={stats['missing']} ({stats['pct_missing']:.1f}%)")
    print(f"  Range: {stats['min']:.3g} to {stats['max']:.3g}")
    print(f"  Mean: {stats['mean']:.3g} ± {stats['sd']:.3g}")
    print(f"  Median (IQR): {stats['median']:.3g} ({stats['q25']:.3g}–{stats['q75']:.3g})")
    print(f"  Skewness: {stats['skewness']:.2f}, Kurtosis: {stats['kurtosis']:.2f}")
    print(f"  Unique values: {stats['unique_count']}")
    
    # By group (if provided)
    if groups is not None:
        groups = pd.Series(groups) if isinstance(groups, np.ndarray) else groups
        print(f"\n  By group:")
        for group in groups.unique():
            mask = groups == group
            subset = data[mask]
            print(f"    {group}: n={mask.sum()}, mean={subset.mean():.3g}, median={subset.median():.3g}")
    
    return stats


def sample_composition(
    data: pd.DataFrame,
    group_col: str,
    additional_cols: Optional[list] = None,
) -> dict:
    """
    Understand sample composition (group sizes, nesting, etc.).
    
    Args:
        data: DataFrame with samples as rows
        group_col: column name for grouping variable
        additional_cols: columns to show cross-tabulation (e.g., ['response', 'treatment'])
    
    Returns:
        dict with group_counts, imbalance_ratio, and cross-tabulation
    """
    groups = data[group_col].value_counts().sort_index()
    
    print(f"\nSample Composition (by {group_col}):")
    print(f"  Total n: {len(data)}")
    for group, count in groups.items():
        pct = count / len(data) * 100
        print(f"    {group}: {count} ({pct:.1f}%)")
    
    # Imbalance ratio
    imbalance = groups.max() / groups.min() if groups.min() > 0 else np.inf
    print(f"  Imbalance ratio (max/min): {imbalance:.2f}")
    
    # Nesting (if applicable)
    if 'subject_id' in data.columns:
        subjects_per_group = data.groupby(group_col)['subject_id'].nunique()
        print(f"\n  Nesting (subjects per group):")
        for group, n_subjects in subjects_per_group.items():
            cells_per_subject = data[data[group_col] == group].groupby('subject_id').size()
            print(f"    {group}: {n_subjects} subjects, " 
                  f"{cells_per_subject.min()}–{cells_per_subject.max()} cells per subject")
    
    # Cross-tabulation (if additional_cols provided)
    if additional_cols:
        print(f"\n  Cross-tabulation:")
        for col in additional_cols:
            if col in data.columns:
                ct = pd.crosstab(data[group_col], data[col], margins=True)
                print(f"\n    {group_col} × {col}:")
                print(ct)
    
    return {
        'group_counts': groups.to_dict(),
        'total_n': len(data),
        'imbalance_ratio': imbalance,
    }


def check_proportions(
    data: pd.DataFrame,
    numerator_col: str,
    numerator_val: Optional[Union[str, bool, int]] = True,
    denominator_col: Optional[str] = None,
    by_col: Optional[str] = None,
) -> pd.DataFrame:
    """
    Calculate proportions and confidence intervals.
    
    Args:
        data: DataFrame
        numerator_col: column name for the category of interest
        numerator_val: value that indicates "positive" (default: True, 1, or "yes")
        denominator_col: column name for denominator (default: all rows)
        by_col: grouping column (e.g., "treatment")
    
    Returns:
        DataFrame with proportions and 95% CIs
    """
    from scipy.stats import binom
    
    if denominator_col is None:
        denominator_col = '__all__'
        data[denominator_col] = True
    
    results = []
    
    if by_col is None:
        # Single proportion
        n_pos = (data[numerator_col] == numerator_val).sum()
        n_total = len(data)
        prop = n_pos / n_total
        ci_lo, ci_hi = binom.interval(0.95, n_total, prop)
        results.append({
            'group': 'overall',
            'n_positive': n_pos,
            'n_total': n_total,
            'proportion': prop,
            'ci_lower': ci_lo / n_total,
            'ci_upper': ci_hi / n_total,
        })
    else:
        # Proportions by group
        for group in data[by_col].unique():
            subset = data[data[by_col] == group]
            n_pos = (subset[numerator_col] == numerator_val).sum()
            n_total = len(subset)
            prop = n_pos / n_total
            ci_lo, ci_hi = binom.interval(0.95, n_total, prop)
            results.append({
                'group': group,
                'n_positive': n_pos,
                'n_total': n_total,
                'proportion': prop,
                'ci_lower': ci_lo / n_total,
                'ci_upper': ci_hi / n_total,
            })
    
    result_df = pd.DataFrame(results)
    print(result_df.to_string(index=False))
    return result_df
