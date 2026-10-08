"""
Visualization utilities for data inspection, color scales, and annotations.

Import common functions:
    from utilities import inspect_variable, categorical_palette, add_comparison_annotation
"""

from .data_inspection import inspect_variable, sample_composition, check_proportions
from .color_scales import categorical_palette, sequential_palette, diverging_palette, check_colorblind
from .annotations import add_comparison_annotation, add_sample_size_labels, add_statistical_summary, annotate_points

__all__ = [
    'inspect_variable',
    'sample_composition',
    'check_proportions',
    'categorical_palette',
    'sequential_palette',
    'diverging_palette',
    'check_colorblind',
    'add_comparison_annotation',
    'add_sample_size_labels',
    'add_statistical_summary',
    'annotate_points',
]
