"""
Utilities for annotations and labels.

Good annotations:
- Name what is being compared
- State the outcome or finding
- Include units if not obvious
- Retain context (n, effect size, p-value)
"""

import matplotlib.pyplot as plt
from typing import Optional, Union, List


def add_comparison_annotation(
    ax: plt.Axes,
    groups: List[str],
    p_value: Optional[float] = None,
    test_name: Optional[str] = None,
    effect_size: Optional[float] = None,
    effect_size_name: Optional[str] = None,
    y_position: Optional[float] = None,
    fontsize: int = 9,
) -> None:
    """
    Add annotation showing comparison outcome.
    
    Args:
        ax: matplotlib axis
        groups: list of group names being compared (e.g., ["responder", "non-responder"])
        p_value: p-value from statistical test
        test_name: name of statistical test (e.g., "Welch's t-test")
        effect_size: effect size (e.g., fold-change)
        effect_size_name: name of effect size (e.g., "fold-change")
        y_position: y-coordinate for annotation (default: top of axis)
        fontsize: font size for annotation
    
    Example:
        add_comparison_annotation(ax, groups=['responder', 'non-responder'],
                                 p_value=0.001, test_name="Welch's t-test")
        → Adds "p < 0.001 (Welch's t-test)" at top of plot
    """
    parts = [f"{groups[0]} vs. {groups[1]}"]
    
    if effect_size is not None:
        if effect_size_name is None:
            effect_size_name = "effect size"
        parts.append(f"{effect_size_name}: {effect_size:.2f}")
    
    if p_value is not None:
        if p_value < 0.001:
            p_str = "p < 0.001"
        elif p_value < 0.01:
            p_str = f"p = {p_value:.3f}"
        else:
            p_str = f"p = {p_value:.2f}"
        
        if test_name:
            parts.append(f"{p_str} ({test_name})")
        else:
            parts.append(p_str)
    
    text = ", ".join(parts)
    
    if y_position is None:
        y_position = ax.get_ylim()[1] * 0.95
    
    ax.text(0.5, y_position, text, 
           transform=ax.transAxes,
           ha='center', va='top',
           fontsize=fontsize,
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


def add_sample_size_labels(
    ax: plt.Axes,
    n_values: dict,
    position: str = "bottom",
    offset: float = 0.02,
    fontsize: int = 8,
) -> None:
    """
    Add sample size labels to groups (bars, boxes, etc.).
    
    Args:
        ax: matplotlib axis
        n_values: dict mapping group name (x-position) to n
        position: "bottom", "top", or "inside"
        offset: offset from axis edge (as fraction of y-range)
        fontsize: font size for labels
    
    Example:
        add_sample_size_labels(ax, {'responder': 24, 'non-responder': 24})
    """
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    
    for x_pos, (group_name, n) in enumerate(n_values.items()):
        if position == "bottom":
            y = ax.get_ylim()[0] + offset * y_range
            va = 'bottom'
        elif position == "top":
            y = ax.get_ylim()[1] - offset * y_range
            va = 'top'
        else:  # inside
            y = ax.get_ylim()[0] + 0.5 * y_range
            va = 'center'
        
        ax.text(x_pos, y, f"n={n}", ha='center', va=va, fontsize=fontsize)


def add_statistical_summary(
    ax: plt.Axes,
    summary_text: str,
    position: str = "upper right",
    fontsize: int = 8,
) -> None:
    """
    Add statistical summary box to figure.
    
    Args:
        ax: matplotlib axis
        summary_text: text to display (can include newlines)
        position: legend position (e.g., "upper right", "lower left")
        fontsize: font size
    
    Example:
        summary = "Median (IQR):\\nResponders: 5.2 (2.1–8.3)\\nNon-responders: 2.1 (0.5–4.2)"
        add_statistical_summary(ax, summary)
    """
    ax.text(0.98, 0.97, summary_text,
           transform=ax.transAxes,
           ha='right', va='top',
           fontsize=fontsize,
           family='monospace',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))


def annotate_points(
    ax: plt.Axes,
    x: list,
    y: list,
    labels: list,
    fontsize: int = 8,
    offset: float = 0.02,
) -> None:
    """
    Annotate points in a scatter plot with labels.
    
    Args:
        ax: matplotlib axis
        x, y: coordinates
        labels: label for each point
        fontsize: font size
        offset: offset from point (in data units)
    """
    for xi, yi, label in zip(x, y, labels):
        ax.annotate(label, xy=(xi, yi), xytext=(offset, offset),
                   textcoords='offset points',
                   fontsize=fontsize,
                   ha='left')
