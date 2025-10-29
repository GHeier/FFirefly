"""
Path plot type for band structure visualization in FFirefly plotting package.
"""

from typing import Optional, Tuple, List
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from ..loaders import load, get_label_from_filename


def BZ_point_to_q(letter: str) -> List[float]:
    """Convert Brillouin zone letter to fractional coordinates.

    Args:
        letter: BZ point letter (g=Gamma, x, m, r)

    Returns:
        Fractional coordinates [qx, qy, qz]
    """
    letter = letter.lower()
    if letter == 'g':
        return [0, 0, 0]
    elif letter == 'x':
        return [1, 0, 0]
    elif letter == 'm':
        return [1, 1, 0]
    elif letter == 'r':
        return [1, 1, 1]
    else:
        raise ValueError(f"BZ point letter '{letter}' not recognized")


def plot_section(
    field_data,
    qi: np.ndarray,
    qf: np.ndarray,
    section: int,
    ax: Axes,
    color: str,
    label: Optional[str] = None,
    N: int = 100
) -> float:
    """Plot one section of a k-path.

    Args:
        field_data: Loaded field data
        qi: Initial k-point (fractional coordinates)
        qf: Final k-point (fractional coordinates)
        section: Section number (for x-axis positioning)
        ax: Axes to plot on
        color: Line color
        label: Plot label
        N: Number of points along path

    Returns:
        Minimum value along this section
    """
    # Default Brillouin zone (cubic)
    BZ = np.array([[2 * np.pi, 0, 0],
                   [0, 2 * np.pi, 0],
                   [0, 0, 2 * np.pi]])

    # Generate path
    t = np.linspace(0, 1, N)
    q = qi[None, :] + t[:, None] * (qf - qi)[None, :]  # shape (N, 3)
    q_cart = (BZ @ q.T).T

    # Evaluate field
    y = field_data(q_cart.tolist())

    # Take real part if complex
    if np.iscomplexobj(y):
        y = np.real(y)

    # X-axis for this section
    x = np.linspace(section - 1, section, N)

    ax.plot(x, y, color=color, label=label)

    return np.min(y)


def plot_path(
    filename: str,
    ax: Optional[Axes] = None,
    fig: Optional[Figure] = None,
    path: str = "gxmg",
    label: Optional[str] = None,
    color: Optional[str] = None,
    hline: bool = False,
    title: Optional[str] = None,
    ylabel: Optional[str] = None,
    **kwargs
) -> Tuple[Figure, Axes, list]:
    """Plot band structure along a k-path.

    Args:
        filename: Path to data file
        ax: Existing axes to plot on (optional)
        fig: Existing figure (optional)
        path: K-path specification (e.g., "gxmg" for Γ→X→M→Γ)
        label: Label for the plot
        color: Line color
        hline: Add horizontal line at y=0
        title: Plot title
        ylabel: Y-axis label
        **kwargs: Additional plotting arguments

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

    if color is None:
        # Get next color from cycler
        color = ax._get_lines.get_next_color()

    artists = []

    # Plot each section of the path
    path_points = list(path.lower())
    qi = np.array(BZ_point_to_q(path_points[0])) / 2.0

    for i, letter in enumerate(path_points[1:], start=1):
        qf = np.array(BZ_point_to_q(letter)) / 2.0
        section_label = label if i == 1 else None
        plot_section(field_data, qi, qf, i, ax, color, label=section_label)
        qi = qf

    # Add horizontal line at zero if requested
    if hline:
        line = ax.axhline(y=0, color='gray', linestyle='--', linewidth=1)
        artists.append(line)

    # Set x-axis to show high-symmetry points
    ax.set_xlim(0, len(path) - 1)
    ax.set_xticks(range(len(path)))

    # Format labels (Γ for G)
    labels = []
    for letter in path:
        if letter.lower() == 'g':
            labels.append('Γ')
        else:
            labels.append(letter.upper())

    ax.set_xticklabels(labels)

    # Set labels
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel('Energy')

    if title is not None:
        ax.set_title(title)

    # Add legend
    if label:
        ax.legend()

    # Add grid
    ax.grid(True, alpha=0.2, axis='y')

    return fig, ax, artists


def plot_path_multi(
    filenames: List[str],
    ax: Optional[Axes] = None,
    fig: Optional[Figure] = None,
    path: str = "gxmg",
    labels: Optional[List[str]] = None,
    colors: Optional[List[str]] = None,
    hline: bool = False,
    title: Optional[str] = None,
    ylabel: Optional[str] = None,
    **kwargs
) -> Tuple[Figure, Axes, list]:
    """Plot multiple band structures on same k-path.

    Args:
        filenames: List of file paths
        ax: Existing axes to plot on (optional)
        fig: Existing figure (optional)
        path: K-path specification
        labels: Labels for each dataset
        colors: Colors for each dataset
        hline: Add horizontal line at y=0
        title: Plot title
        ylabel: Y-axis label
        **kwargs: Additional plotting arguments

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

    for i, (filename, label) in enumerate(zip(filenames, labels)):
        color = colors[i] if colors else None
        _, _, file_artists = plot_path(
            filename, ax=ax, fig=fig, path=path, label=label,
            color=color, hline=(hline and i == 0),
            ylabel=ylabel, **kwargs
        )
        artists.extend(file_artists)

    if title is not None:
        ax.set_title(title)

    return fig, ax, artists
