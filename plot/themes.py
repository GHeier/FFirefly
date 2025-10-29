"""
Theme management for FFirefly plotting package.
Provides predefined themes and custom theme support.
"""

from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List
import matplotlib.pyplot as plt
from cycler import cycler


@dataclass
class Theme:
    """Theme configuration for plots."""

    name: str
    colors: List[str]
    linewidth: float = 1.5
    markersize: float = 5.0
    figsize: tuple = (8, 6)
    dpi: int = 100
    font_size: int = 11
    font_family: str = 'sans-serif'
    grid_alpha: float = 0.2
    grid_style: str = '--'
    background_color: str = 'white'
    axes_background: str = 'white'
    text_color: str = 'black'
    axes_edge_color: str = 'black'
    tick_color: str = 'black'
    cmap: str = 'viridis'
    cmap_diverging: str = 'RdBu_r'
    legend_frameon: bool = True
    legend_framealpha: float = 0.8

    def apply(self) -> None:
        """Apply theme to matplotlib rcParams."""
        plt.rcParams.update({
            'axes.prop_cycle': cycler(color=self.colors),
            'lines.linewidth': self.linewidth,
            'lines.markersize': self.markersize,
            'figure.figsize': self.figsize,
            'figure.dpi': self.dpi,
            'font.size': self.font_size,
            'font.family': self.font_family,
            'grid.alpha': self.grid_alpha,
            'grid.linestyle': self.grid_style,
            'figure.facecolor': self.background_color,
            'axes.facecolor': self.axes_background,
            'text.color': self.text_color,
            'axes.labelcolor': self.text_color,
            'axes.edgecolor': self.axes_edge_color,
            'xtick.color': self.tick_color,
            'ytick.color': self.tick_color,
            'legend.frameon': self.legend_frameon,
            'legend.framealpha': self.legend_framealpha,
        })


# Predefined themes
THEMES = {}

THEMES['firefly'] = Theme(
    name='firefly',
    colors=[
        "#9a05fc",  # Purple
        "#f00524",  # Red
        "#f80af1",  # Pink
        "#0a68f8",  # Blue
        "#1c841f",  # Green
        "#fc9303",  # Orange
        "#865522",  # Brown
        "#00c7a9",  # Teal
        "#C9C22A",  # Yellow
        "#7f7f7f",  # Gray
        "black"
    ],
    linewidth=1.5,
    markersize=5.0,
    cmap='viridis',
    cmap_diverging='RdBu_r'
)

THEMES['dark'] = Theme(
    name='dark',
    colors=[
        "#bb86fc",  # Purple
        "#cf6679",  # Pink/Red
        "#03dac6",  # Teal
        "#80cbc4",  # Light cyan
        "#ffd54f",  # Yellow
        "#ff8a65",  # Orange
        "#90caf9",  # Blue
        "#a5d6a7",  # Green
        "#ce93d8",  # Light purple
        "#bcaaa4",  # Brown gray
    ],
    linewidth=2.0,
    markersize=6.0,
    background_color='#121212',
    axes_background='#1e1e1e',
    text_color='#e0e0e0',
    axes_edge_color='#424242',
    tick_color='#e0e0e0',
    grid_alpha=0.15,
    cmap='plasma',
    cmap_diverging='RdBu_r',
    legend_framealpha=0.9
)

THEMES['publication'] = Theme(
    name='publication',
    colors=[
        "#000000",  # Black
        "#e41a1c",  # Red
        "#377eb8",  # Blue
        "#4daf4a",  # Green
        "#984ea3",  # Purple
        "#ff7f00",  # Orange
        "#a65628",  # Brown
        "#f781bf",  # Pink
    ],
    linewidth=1.2,
    markersize=4.0,
    figsize=(6, 4),
    dpi=300,
    font_size=10,
    font_family='serif',
    grid_alpha=0.25,
    cmap='viridis',
    cmap_diverging='RdBu_r',
    legend_frameon=True,
    legend_framealpha=1.0
)

THEMES['presentation'] = Theme(
    name='presentation',
    colors=[
        "#e41a1c",  # Red
        "#377eb8",  # Blue
        "#4daf4a",  # Green
        "#984ea3",  # Purple
        "#ff7f00",  # Orange
        "#ffff33",  # Yellow
        "#a65628",  # Brown
        "#f781bf",  # Pink
    ],
    linewidth=3.0,
    markersize=8.0,
    figsize=(12, 8),
    dpi=150,
    font_size=16,
    grid_alpha=0.3,
    cmap='viridis',
    cmap_diverging='RdBu_r',
    legend_frameon=True,
    legend_framealpha=0.95
)

THEMES['minimal'] = Theme(
    name='minimal',
    colors=[
        "#2E3440",  # Dark gray
        "#5E81AC",  # Blue
        "#BF616A",  # Red
        "#A3BE8C",  # Green
        "#B48EAD",  # Purple
        "#D08770",  # Orange
    ],
    linewidth=1.5,
    markersize=4.0,
    figsize=(8, 5),
    font_size=11,
    grid_alpha=0.15,
    background_color='#ECEFF4',
    axes_background='white',
    cmap='cividis',
    cmap_diverging='RdBu_r'
)


def get_theme(name: str = 'default') -> Theme:
    """Get a theme by name.

    Args:
        name: Theme name (default, dark, publication, presentation, minimal)

    Returns:
        Theme object

    Raises:
        ValueError: If theme name not found
    """
    if name not in THEMES:
        available = ', '.join(THEMES.keys())
        raise ValueError(f"Theme '{name}' not found. Available: {available}")
    return THEMES[name]


def set_theme(name: str = 'default') -> None:
    """Set the active plotting theme.

    Args:
        name: Theme name to activate
    """
    theme = get_theme(name)
    theme.apply()


def register_theme(theme: Theme) -> None:
    """Register a custom theme.

    Args:
        theme: Custom Theme object
    """
    THEMES[theme.name] = theme
