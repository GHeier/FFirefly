"""
Bar plot type for FFirefly plotting package.
"""

from typing import Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ..loaders import load, get_label_from_filename


def plot_bar(
    filename: str,
    ax: Optional[Axes] = None,
    fig: Optional[Figure] = None,
    label: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    **kwargs
) -> Tuple[Figure, Axes, list]:
    """Create a bar plot from field data.

    Args:
        filename: Path to data file
        ax: Existing axes to plot on (optional)
        fig: Existing figure (optional)
        label: Label for the plot
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        **kwargs: Additional arguments passed to ax.bar()

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

    if label is None:
        label = get_label_from_filename(filename)

    # Extract values
    if field_data.field_type == 'text':
        data = field_data.data
        x = data[:, 0]
        y = data[:, 1]
    else:
        raise NotImplementedError("Bar plots from field objects not yet implemented")

    # Create bar plot
    bars = ax.bar(x, y, label=label, **kwargs)

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
    ax.grid(True, alpha=0.2, axis='y')

    return fig, ax, [bars]
