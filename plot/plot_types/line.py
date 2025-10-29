"""
Line plot type for FFirefly plotting package.
"""

from typing import Optional, Tuple, Any
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ..loaders import FieldData, load, get_label_from_filename


def plot_line(
    filename: str,
    ax: Optional[Axes] = None,
    fig: Optional[Figure] = None,
    label: Optional[str] = None,
    plot_real: bool = True,
    plot_imag: bool = True,
    plot_magnitude: bool = False,
    plot_phase: bool = False,
    use_w_points: bool = True,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    **kwargs
) -> Tuple[Figure, Axes, list]:
    """Create a line plot from field data.

    Args:
        filename: Path to data file
        ax: Existing axes to plot on (optional)
        fig: Existing figure (optional)
        label: Label for the plot
        plot_real: Plot real part (for complex fields)
        plot_imag: Plot imaginary part (for complex fields)
        plot_magnitude: Plot magnitude (for complex fields)
        plot_phase: Plot phase (for complex fields)
        use_w_points: Use frequency points if available
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        **kwargs: Additional arguments passed to ax.plot()

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
        # Text file: assume first column is x, rest are y values
        data = field_data.data
        x = data[:, 0]

        for col_idx in range(1, data.shape[1]):
            y = data[:, col_idx]
            col_label = f"{label} col {col_idx}" if data.shape[1] > 2 else label
            line, = ax.plot(x, y, label=col_label, **kwargs)
            artists.append(line)

    elif field_data.is_matrix:
        raise NotImplementedError("Matrix field line plots not yet supported")

    else:
        # Scalar field
        if use_w_points and field_data.has_frequency:
            # Use frequency points
            x = field_data.w_points
            y = field_data(x)

            if field_data.is_complex:
                if plot_real:
                    line, = ax.plot(x, y.real, label=f"Re({label})", **kwargs)
                    artists.append(line)
                if plot_imag:
                    line, = ax.plot(x, y.imag, label=f"Im({label})", **kwargs)
                    artists.append(line)
                if plot_magnitude:
                    line, = ax.plot(x, np.abs(y), label=f"|{label}|", **kwargs)
                    artists.append(line)
                if plot_phase:
                    line, = ax.plot(x, np.angle(y), label=f"arg({label})", **kwargs)
                    artists.append(line)
            else:
                line, = ax.plot(x, y, label=label, **kwargs)
                artists.append(line)

        else:
            # 0D or 1D field without explicit frequency axis
            if field_data.dimension == 0:
                # 0D field with frequency
                if field_data.has_frequency:
                    x = field_data.w_points
                    y = field_data(x)

                    if field_data.is_complex:
                        if plot_real:
                            line, = ax.plot(x, y.real, label=f"Re({label})", **kwargs)
                            artists.append(line)
                        if plot_imag:
                            line, = ax.plot(x, y.imag, label=f"Im({label})", **kwargs)
                            artists.append(line)
                    else:
                        line, = ax.plot(x, y, label=label, **kwargs)
                        artists.append(line)
                else:
                    raise ValueError("Cannot plot 0D field without frequency points")

            elif field_data.dimension == 1:
                # 1D field: plot along k-direction
                nx = field_data.mesh[0] if field_data.mesh else 100
                x = np.linspace(0, 1, nx)

                if len(field_data.domain) > 0:
                    # Convert to real space if domain available
                    domain_vec = field_data.domain[0] if len(field_data.domain.shape) > 1 else field_data.domain
                    x_real = x * np.linalg.norm(domain_vec)
                else:
                    x_real = x

                # Evaluate field
                points = [[xi] for xi in x]
                y = field_data(points)

                if field_data.is_complex:
                    if plot_real:
                        line, = ax.plot(x_real, y.real, label=f"Re({label})", **kwargs)
                        artists.append(line)
                    if plot_imag:
                        line, = ax.plot(x_real, y.imag, label=f"Im({label})", **kwargs)
                        artists.append(line)
                else:
                    line, = ax.plot(x_real, y, label=label, **kwargs)
                    artists.append(line)
            else:
                raise ValueError(f"Cannot create line plot from {field_data.dimension}D field")

    # Set labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    elif use_w_points and field_data.has_frequency:
        ax.set_xlabel(r'$\omega$')

    if ylabel is not None:
        ax.set_ylabel(ylabel)

    if title is not None:
        ax.set_title(title)

    # Add legend if we have multiple artists
    if len(artists) > 1 or (len(artists) == 1 and artists[0].get_label()):
        ax.legend()

    # Add grid
    ax.grid(True, alpha=0.2)

    return fig, ax, artists
