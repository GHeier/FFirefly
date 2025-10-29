"""
FFirefly Plotting Package

A modern, object-oriented plotting package for visualizing FFirefly computational results.

Features:
- Multiple plot types: line, scatter, histogram, bar, heatmap, path, frequency
- Smart data loading for all field types (Field_R, Field_C, Field_RM, Field_CM)
- Theme system with predefined themes
- Plotter class for building complex multi-panel plots
- Returns (fig, ax, artists) for full control

Example:
    >>> import plot
    >>> fig, ax, artists = plot.sketch('data.h5', 'line')
    >>> plt.show()

    >>> # Multi-panel plot
    >>> p = plot.Plotter(theme='dark')
    >>> p.add_subplot(2, 1, 1, 'line', 'G_iw.h5')
    >>> p.add_subplot(2, 1, 2, 'frequency', 'Sigma_iw.h5')
    >>> fig, axes, artists = p.build()
    >>> plt.show()
"""

from typing import Optional, Tuple, Any
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from .plotter import Plotter
from .themes import get_theme, set_theme, register_theme, Theme, THEMES
from .loaders import load, FieldData, get_label_from_filename

# Import plot type modules for direct access
from .plot_types import line, scatter, histogram, bar, heatmap, path, frequency


def sketch(
    file: str,
    plot_type: str = 'line',
    ax: Optional[Axes] = None,
    fig: Optional[Figure] = None,
    theme: Optional[str] = None,
    **kwargs
) -> Tuple[Figure, Axes, list]:
    """
    Create a plot from a data file.

    This is the main entry point for quick plotting. For more complex multi-panel
    plots, use the Plotter class.

    Args:
        file: Path to data file (.h5, .hdf5, .dat, .txt, .csv)
        plot_type: Type of plot to create:
            - 'line': Line plot (1D data, frequency data)
            - 'scatter': Scatter plot
            - 'histogram': Histogram
            - 'bar': Bar chart
            - 'heatmap': 2D heatmap/colormap
            - 'path': Band structure along k-path
            - 'frequency': Frequency-dependent plot (G, Σ, etc.)
        ax: Existing axes to plot on (optional)
        fig: Existing figure (optional)
        theme: Theme to use ('firefly', 'dark', 'publication', 'presentation', 'minimal')
        **kwargs: Additional arguments passed to the specific plot type

    Returns:
        (fig, ax, artists) tuple where:
            - fig: matplotlib Figure object
            - ax: matplotlib Axes object (or list of Axes for multi-panel)
            - artists: list of plot artists (lines, patches, collections, etc.)

    Examples:
        >>> # Simple line plot
        >>> fig, ax, artists = sketch('data.dat', 'line')

        >>> # Frequency plot with custom styling
        >>> fig, ax, artists = sketch('G_iw.h5', 'frequency',
        ...                          plot_real=True, plot_imag=True,
        ...                          marker='o', markersize=4)

        >>> # Heatmap with custom range
        >>> fig, ax, artists = sketch('field.h5', 'heatmap',
        ...                          xlim=(-np.pi, np.pi), ylim=(-np.pi, np.pi),
        ...                          resolution=300, cmap='RdBu_r')

        >>> # Band structure plot
        >>> fig, ax, artists = sketch('bands.h5', 'path',
        ...                          path='gxmg', hline=True)
    """
    # Apply theme if specified
    if theme is not None:
        set_theme(theme)

    # Dispatch to appropriate plot type
    plot_type = plot_type.lower()

    if plot_type == 'line':
        return line.plot_line(file, ax=ax, fig=fig, **kwargs)
    elif plot_type == 'scatter':
        return scatter.plot_scatter(file, ax=ax, fig=fig, **kwargs)
    elif plot_type in ['histogram', 'hist']:
        return histogram.plot_histogram(file, ax=ax, fig=fig, **kwargs)
    elif plot_type == 'bar':
        return bar.plot_bar(file, ax=ax, fig=fig, **kwargs)
    elif plot_type in ['heatmap', 'colormap']:
        return heatmap.plot_heatmap(file, ax=ax, fig=fig, **kwargs)
    elif plot_type in ['path', 'band', 'bands']:
        return path.plot_path(file, ax=ax, fig=fig, **kwargs)
    elif plot_type in ['frequency', 'freq']:
        return frequency.plot_frequency(file, ax=ax, fig=fig, **kwargs)
    else:
        raise ValueError(
            f"Unknown plot type: {plot_type}. "
            f"Available: line, scatter, histogram, bar, heatmap, path, frequency"
        )


# Convenience aliases
plot = sketch


__all__ = [
    # Main functions
    'sketch',
    'plot',

    # Plotter class
    'Plotter',

    # Theme management
    'get_theme',
    'set_theme',
    'register_theme',
    'Theme',
    'THEMES',

    # Data loading
    'load',
    'FieldData',
    'get_label_from_filename',

    # Plot type modules (for advanced usage)
    'line',
    'scatter',
    'histogram',
    'bar',
    'heatmap',
    'path',
    'frequency',
]

__version__ = '1.0.0'
__author__ = 'FFirefly Team'
