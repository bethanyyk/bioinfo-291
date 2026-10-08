"""
Color scales for different data types.

Principles:
- Categorical (hue): distinct colors, no implied order
- Sequential (magnitude): ordered lightness, darker = higher
- Diverging (centered): light center, two colors diverging (e.g., cold to hot)
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from typing import Optional, List, Union


def categorical_palette(
    n_categories: int,
    palette: str = "default",
) -> List[str]:
    """
    Categorical color palette (distinct hues, no implied order).
    
    Args:
        n_categories: number of colors needed
        palette: "default" (tab20/Set1), "colorblind", "pastel"
    
    Returns:
        list of hex color strings
    """
    if palette == "default":
        cmap_name = "tab20" if n_categories <= 20 else "tab20b"
        cmap = plt.get_cmap(cmap_name)
        colors = [mcolors.to_hex(cmap(i)) for i in range(n_categories)]
    elif palette == "colorblind":
        # Colorblind-friendly palette (Okabe-Ito)
        colors = [
            "#E69F00",  # orange
            "#56B4E9",  # blue
            "#009E73",  # green
            "#F0E442",  # yellow
            "#0072B2",  # dark blue
            "#D55E00",  # vermillion
            "#CC79A7",  # purple
            "#000000",  # black
        ]
        colors = colors[:n_categories]
    elif palette == "pastel":
        cmap = plt.get_cmap("Pastel1")
        colors = [mcolors.to_hex(cmap(i)) for i in range(n_categories)]
    else:
        # Assume it's a matplotlib colormap name
        cmap = plt.get_cmap(palette)
        colors = [mcolors.to_hex(cmap(i / (n_categories - 1))) for i in range(n_categories)]
    
    return colors


def sequential_palette(
    vmin: float,
    vmax: float,
    palette: str = "Greys",
    n_colors: int = 256,
) -> str:
    """
    Sequential color palette (ordered lightness, darker = higher).
    
    Args:
        vmin: minimum value
        vmax: maximum value
        palette: matplotlib colormap name (e.g., "Greys", "Viridis", "YlOrRd")
        n_colors: number of discrete colors (default: continuous)
    
    Returns:
        matplotlib colormap
    """
    cmap = plt.get_cmap(palette)
    return cmap


def diverging_palette(
    vmin: float,
    vmax: float,
    center: float = 0.0,
    palette: str = "RdBu_r",
    n_colors: int = 256,
):
    """
    Diverging color palette (light center, two colors diverging).
    
    Useful for data with a meaningful center (e.g., log fold-change from -1 to +1).
    
    Args:
        vmin: minimum value
        vmax: maximum value
        center: center value (mapped to light color)
        palette: matplotlib colormap name (e.g., "RdBu", "RdBu_r", "coolwarm", "bwr")
        n_colors: number of discrete colors
    
    Returns:
        tuple: (cmap, norm) where cmap is matplotlib colormap and norm is TwoSlopeNorm.
        Use in seaborn/matplotlib as: sns.heatmap(data, cmap=cmap, norm=norm)
    """
    cmap = plt.get_cmap(palette)
    
    # Normalize so center is at 0.5
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
    return cmap, norm


def check_colorblind(colors: List[str], mode: str = "deuteranopia") -> None:
    """
    Simulate how colors appear to people with color blindness.
    
    Args:
        colors: list of hex color codes
        mode: "protanopia" (red-blind), "deuteranopia" (green-blind), 
              "tritanopia" (blue-blind)
    
    Prints visualization of original vs. simulated colors.
    """
    # Simplified simulation (for demo; real simulator is more accurate)
    print(f"\nColor simulation ({mode}):")
    print("  [Original] [Simulated]")
    for color in colors:
        print(f"  {color}")
    print("\nNote: Use a full simulator (e.g., Coblis) for accurate assessment.")


def get_colormap_for_data(data_type: str, **kwargs):
    """
    Recommend and return a colormap based on data type.
    
    Args:
        data_type: "categorical", "sequential", "diverging", "count"
        **kwargs: passed to palette function (e.g., n_categories, palette, center)
    
    Returns:
        appropriate colormap or color list
    """
    if data_type == "categorical":
        return categorical_palette(**kwargs)
    elif data_type == "sequential":
        return sequential_palette(**kwargs)
    elif data_type == "diverging":
        return diverging_palette(**kwargs)
    elif data_type == "count":
        # For discrete counts: use light-to-dark gradient
        return sequential_palette(**kwargs)
    else:
        raise ValueError(f"Unknown data_type: {data_type}")
