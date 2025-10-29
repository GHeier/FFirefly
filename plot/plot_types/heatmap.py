"""
Heatmap/colormap plot type for FFirefly plotting package.
"""

from typing import Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ..loaders import load, get_label_from_filename


def plot_heatmap(
    filename: str,
    ax: Optional[Axes] = None,
    fig: Optional[Figure] = None,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    xlim: Tuple[float, float] = (-np.pi, np.pi),
    ylim: Tuple[float, float] = (-np.pi, np.pi),
    resolution: int = 200,
    cmap: Optional[str] = None,
    component: str = 'real',  # 'real', 'imag', 'magnitude', 'phase'
    colorbar: bool = True,
    w: float = 0.0,  # Frequency value
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    **kwargs
) -> Tuple[Figure, Axes, list]:
    """Create a 2D heatmap/colormap from field data.

    Args:
        filename: Path to data file
        ax: Existing axes to plot on (optional)
        fig: Existing figure (optional)
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        xlim: X-axis limits
        ylim: Y-axis limits
        resolution: Number of points in each direction
        cmap: Colormap name
        component: Which component to plot ('real', 'imag', 'magnitude', 'phase')
        colorbar: Whether to add a colorbar
        w: Frequency value to evaluate at
        vmin: Minimum value for colormap
        vmax: Maximum value for colormap
        **kwargs: Additional arguments passed to ax.pcolormesh()

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

    if field_data.dimension != 2:
        raise ValueError(f"Heatmap requires 2D field, got {field_data.dimension}D")

    # Create meshgrid
    x = np.linspace(*xlim, resolution)
    y = np.linspace(*ylim, resolution)
    X, Y = np.meshgrid(x, y)

    # Evaluate field on grid
    points = []
    for i in range(resolution):
        for j in range(resolution):
            points.append([X[i, j], Y[i, j]])

    values = field_data(points, w)
    Z = values.reshape(resolution, resolution)

    # Extract component
    if field_data.is_complex:
        if component == 'real':
            Z = Z.real
            comp_label = 'Re'
        elif component == 'imag':
            Z = Z.imag
            comp_label = 'Im'
        elif component == 'magnitude':
            Z = np.abs(Z)
            comp_label = '|·|'
        elif component == 'phase':
            Z = np.angle(Z)
            comp_label = 'arg'
        else:
            raise ValueError(f"Unknown component: {component}")
    else:
        comp_label = ''

    # Auto-select colormap based on data
    if cmap is None:
        if vmin is None and vmax is None:
            z_min, z_max = np.min(Z), np.max(Z)
        else:
            z_min = vmin if vmin is not None else np.min(Z)
            z_max = vmax if vmax is not None else np.max(Z)

        # If data crosses zero, use diverging colormap
        if z_min * z_max < 0:
            cmap = 'RdBu_r'
        else:
            cmap = 'viridis'

    # Create heatmap
    mesh = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cmap,
                        vmin=vmin, vmax=vmax, **kwargs)

    # Add colorbar
    if colorbar:
        label = get_label_from_filename(filename)
        cbar_label = f"{comp_label}({label})" if comp_label else label
        fig.colorbar(mesh, ax=ax, label=cbar_label)

    # Set labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel('$k_x$')

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel('$k_y$')

    if title is not None:
        ax.set_title(title)
    else:
        label = get_label_from_filename(filename)
        ax.set_title(f"{label}")

    ax.set_aspect('equal')

    return fig, ax, [mesh]
