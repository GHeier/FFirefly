"""
Plotter class for building complex multi-panel plots.
"""

from typing import Optional, Tuple, Any, Dict, List
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from .themes import get_theme, set_theme
from .plot_types import line, scatter, histogram, bar, heatmap, path, frequency


class Plotter:
    """
    A class for building complex multi-panel plots.

    Example:
        >>> p = Plotter(theme='dark')
        >>> p.add_subplot(2, 2, 1, plot_type='line', file='data1.h5')
        >>> p.add_subplot(2, 2, 2, plot_type='heatmap', file='data2.h5')
        >>> fig, axes, artists = p.build()
        >>> plt.show()
    """

    def __init__(
        self,
        figsize: Optional[Tuple[float, float]] = None,
        theme: str = 'default',
        dpi: Optional[int] = None
    ):
        """Initialize Plotter.

        Args:
            figsize: Figure size (width, height) in inches
            theme: Theme name to use
            dpi: Dots per inch for figure
        """
        self.theme_name = theme
        self.theme = get_theme(theme)

        # Apply theme
        self.theme.apply()

        # Override figsize and dpi if provided
        if figsize is not None:
            self.figsize = figsize
        else:
            self.figsize = self.theme.figsize

        if dpi is not None:
            self.dpi = dpi
        else:
            self.dpi = self.theme.dpi

        self.subplots = []
        self.fig: Optional[Figure] = None
        self.axes: List[Axes] = []
        self.artists: List[Any] = []

    def add_subplot(
        self,
        nrows: int,
        ncols: int,
        index: int,
        plot_type: str,
        file: str,
        **kwargs
    ) -> 'Plotter':
        """Add a subplot to the figure.

        Args:
            nrows: Number of rows in subplot grid
            ncols: Number of columns in subplot grid
            index: Index of this subplot (1-based)
            plot_type: Type of plot ('line', 'scatter', 'heatmap', etc.)
            file: Path to data file
            **kwargs: Additional arguments for the plot type

        Returns:
            self (for method chaining)
        """
        self.subplots.append({
            'nrows': nrows,
            'ncols': ncols,
            'index': index,
            'plot_type': plot_type,
            'file': file,
            'kwargs': kwargs
        })
        return self

    def build(self) -> Tuple[Figure, List[Axes], List[Any]]:
        """Build the figure with all subplots.

        Returns:
            (fig, axes, artists) tuple
        """
        if not self.subplots:
            raise ValueError("No subplots added. Use add_subplot() first.")

        # Determine grid size from subplots
        max_rows = max(sp['nrows'] for sp in self.subplots)
        max_cols = max(sp['ncols'] for sp in self.subplots)

        # Create figure
        self.fig = plt.figure(figsize=self.figsize, dpi=self.dpi)

        # Create each subplot
        for subplot_spec in self.subplots:
            ax = self.fig.add_subplot(
                subplot_spec['nrows'],
                subplot_spec['ncols'],
                subplot_spec['index']
            )

            # Create the plot
            _, _, subplot_artists = self._create_plot(
                ax,
                subplot_spec['plot_type'],
                subplot_spec['file'],
                **subplot_spec['kwargs']
            )

            self.axes.append(ax)
            self.artists.extend(subplot_artists)

        # Adjust layout
        self.fig.tight_layout()

        return self.fig, self.axes, self.artists

    def _create_plot(
        self,
        ax: Axes,
        plot_type: str,
        file: str,
        **kwargs
    ) -> Tuple[Figure, Axes, list]:
        """Create a single plot on given axes.

        Args:
            ax: Axes to plot on
            plot_type: Type of plot
            file: Path to data file
            **kwargs: Plot-specific arguments

        Returns:
            (fig, ax, artists) tuple
        """
        plot_type = plot_type.lower()

        if plot_type == 'line':
            return line.plot_line(file, ax=ax, **kwargs)
        elif plot_type == 'scatter':
            return scatter.plot_scatter(file, ax=ax, **kwargs)
        elif plot_type == 'histogram' or plot_type == 'hist':
            return histogram.plot_histogram(file, ax=ax, **kwargs)
        elif plot_type == 'bar':
            return bar.plot_bar(file, ax=ax, **kwargs)
        elif plot_type == 'heatmap' or plot_type == 'colormap':
            return heatmap.plot_heatmap(file, ax=ax, **kwargs)
        elif plot_type == 'path' or plot_type == 'band':
            return path.plot_path(file, ax=ax, **kwargs)
        elif plot_type == 'frequency' or plot_type == 'freq':
            return frequency.plot_frequency(file, ax=ax, **kwargs)
        else:
            raise ValueError(f"Unknown plot type: {plot_type}")

    def show(self):
        """Build and show the plot."""
        if self.fig is None:
            self.build()
        plt.show()

    def save(self, filename: str, dpi: Optional[int] = None, **kwargs):
        """Build and save the plot to file.

        Args:
            filename: Output filename
            dpi: DPI for saved figure (default: 300)
            **kwargs: Additional arguments passed to fig.savefig()
        """
        if self.fig is None:
            self.build()

        if dpi is None:
            dpi = 300

        self.fig.savefig(filename, dpi=dpi, bbox_inches='tight', **kwargs)
        print(f"Saved to {filename}")
