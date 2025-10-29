"""
Scatter plot type for FFirefly plotting package.
"""

from typing import Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ..loaders import FieldData, load, get_label_from_filename


def plot_scatter(
    filename: str,
    ax: Optional[Axes] = None,
    fig: Optional[Figure] = None,
    label: Optional[str] = None,
    color_by: Optional[str] = None,  # 'magnitude', 'phase', 'real', 'imag'
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    **kwargs
) -> Tuple[Figure, Axes, list]:
    """Create a scatter plot from field data.

    Args:
        filename: Path to data file
        ax: Existing axes to plot on (optional)
        fig: Existing figure (optional)
        label: Label for the plot
        color_by: Color scatter points by this quantity (for complex data)
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        **kwargs: Additional arguments passed to ax.scatter()

    Returns:
        (fig, ax, artists) tuple
    """
    # Create figure/axes if not provided
    if ax is None:
        if fig is None:
            fig, ax = plt.subplots()
        else:
            ax = fig.add_subplot(111)
    elif fig is None:
        fig = ax.get_figure()

    # Load data
    field_data = load(filename)

    artists = []

    if label is None:
        label = get_label_from_filename(filename)

    # Handle different field types
    if field_data.field_type == 'text':
        # Text file: assume first column is x, second is y (or values for coloring)
        data = field_data.data

        if data.shape[1] == 2:
            x, y = data[:, 0], data[:, 1]
            scatter = ax.scatter(x, y, label=label, **kwargs)
            artists.append(scatter)
        elif data.shape[1] >= 3:
            x, y, c = data[:, 0], data[:, 1], data[:, 2]
            scatter = ax.scatter(x, y, c=c, label=label, **kwargs)
            artists.append(scatter)
            fig.colorbar(scatter, ax=ax)
        else:
            raise ValueError("Scatter plot requires at least 2 columns")

    else:
        raise NotImplementedError("Scatter plots from field objects not yet fully implemented")

    # Set labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)

    # Add legend
    if label:
        ax.legend()

    # Add grid
    ax.grid(True, alpha=0.2)

    return fig, ax, artists
