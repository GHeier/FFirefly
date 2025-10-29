"""
Frequency plot type for FFirefly plotting package.
Specialized for Green's functions, self-energies, and other frequency-dependent quantities.
"""

from typing import Optional, Tuple, List
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ..loaders import load, get_label_from_filename


def plot_frequency(
    filename: str,
    ax: Optional[Axes] = None,
    fig: Optional[Figure] = None,
    label: Optional[str] = None,
    plot_real: bool = True,
    plot_imag: bool = True,
    plot_magnitude: bool = False,
    plot_phase: bool = False,
    positive_only: bool = True,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    marker: str = 'o',
    markersize: float = 3,
    **kwargs
) -> Tuple[Figure, Axes, list]:
    """Create frequency-dependent plots (Green's function, self-energy, etc.).

    Args:
        filename: Path to data file
        ax: Existing axes to plot on (optional)
        fig: Existing figure (optional)
        label: Label for the plot
        plot_real: Plot real part
        plot_imag: Plot imaginary part
        plot_magnitude: Plot magnitude
        plot_phase: Plot phase
        positive_only: Only plot positive frequencies
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        marker: Marker style
        markersize: Marker size
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

    if not field_data.has_frequency:
        raise ValueError(f"Field does not have frequency data")

    if label is None:
        label = get_label_from_filename(filename)

    # Get frequency points
    w = field_data.w_points
    if positive_only:
        w = w[w > 0]

    # Evaluate field
    if field_data.dimension == 0:
        # 0D field: frequency only
        y = field_data(w)
    else:
        # Evaluate at k=0 by default
        k_zero = [0.0] * field_data.dimension
        y = np.array([field_data(k_zero, wi) for wi in w])

    artists = []

    # Plot components
    if field_data.is_complex:
        if plot_real and np.max(np.abs(y.real)) > 1e-8:
            line, = ax.plot(w, y.real, marker=marker, markersize=markersize,
                          label=f'Re({label})', **kwargs)
            artists.append(line)

        if plot_imag and np.max(np.abs(y.imag)) > 1e-8:
            line, = ax.plot(w, y.imag, marker=marker, markersize=markersize,
                          label=f'Im({label})', **kwargs)
            artists.append(line)

        if plot_magnitude:
            line, = ax.plot(w, np.abs(y), marker=marker, markersize=markersize,
                          label=f'|{label}|', **kwargs)
            artists.append(line)

        if plot_phase:
            line, = ax.plot(w, np.angle(y), marker=marker, markersize=markersize,
                          label=f'arg({label})', **kwargs)
            artists.append(line)
    else:
        line, = ax.plot(w, y, marker=marker, markersize=markersize,
                       label=label, **kwargs)
        artists.append(line)

    # Set labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    else:
        ax.set_xlabel(r'$\omega$')

    if ylabel is not None:
        ax.set_ylabel(ylabel)

    if title is not None:
        ax.set_title(title)

    # Add legend
    if len(artists) > 0:
        ax.legend()

    # Add grid
    ax.grid(True, alpha=0.2)

    return fig, ax, artists


def plot_frequency_multi(
    filenames: List[str],
    ax: Optional[Axes] = None,
    fig: Optional[Figure] = None,
    labels: Optional[List[str]] = None,
    plot_real: bool = True,
    plot_imag: bool = True,
    positive_only: bool = True,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    **kwargs
) -> Tuple[Figure, Axes, list]:
    """Plot multiple frequency-dependent datasets on same axes.

    Args:
        filenames: List of file paths
        ax: Existing axes to plot on (optional)
        fig: Existing figure (optional)
        labels: Labels for each dataset
        plot_real: Plot real parts
        plot_imag: Plot imaginary parts
        positive_only: Only plot positive frequencies
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

    if labels is None:
        labels = [get_label_from_filename(f) for f in filenames]

    artists = []

    for filename, label in zip(filenames, labels):
        _, _, file_artists = plot_frequency(
            filename, ax=ax, fig=fig, label=label,
            plot_real=plot_real, plot_imag=plot_imag,
            positive_only=positive_only,
            xlabel=xlabel, ylabel=ylabel,
            **kwargs
        )
        artists.extend(file_artists)

    if title is not None:
        ax.set_title(title)

    return fig, ax, artists
