import mimetypes
import firefly as fly
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.ticker import MaxNLocator

def load_theme():
    plt.rcParams["axes.prop_cycle"] = cycler(color=[
        "blue", # Blue
        "red", # Red
        "#8b00e6", # Purple
        "green", # Green
        "#ff7700", # Orange
        "#00c7a9", # Teal (adds cool contrast)
        "#7f7f7f", # Gray (neutral)
        "#f80af1", # Pink
        "#865522", # Brown
        "#C9C22A", # Yellow (bright mid tone)
        "black"
    ])
    plt.rcParams["lines.linewidth"] = 1.0

    plt.rcParams['lines.markersize'] = 5.0
    plt.rcParams["scatter.marker"] = '.'

    plt.rcParams["axes.grid"] = True
    plt.rcParams["grid.color"] = "0.85"
    plt.rcParams["grid.linewidth"] = 0.8

    plt.rcParams.update({
        "axes.titlesize": 18,
        "axes.labelsize": 16,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
    })

    _old_axes_init = plt.Axes.__init__

    def _axes_init_with_sparse_grid(self, *args, **kwargs):
        _old_axes_init(self, *args, **kwargs)
        self.xaxis.set_major_locator(MaxNLocator(nbins=5))
        self.yaxis.set_major_locator(MaxNLocator(nbins=5))
        self.minorticks_off()

    plt.Axes.__init__ = _axes_init_with_sparse_grid


def get_file_type(file_path):
    mime_type, encoding = mimetypes.guess_type(file_path)
    return encoding


def main(files, plot_type, ax=None):
    if get_file_type == "h5":
        data_field = fly.Field()
    pass

