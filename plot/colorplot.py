import numpy as np
import matplotlib.pyplot as plt
import firefly as fly
import firefly.config as cfg

def plot_colorgrid(files, **kwargs):
    field = fly.Field_C(files[0])
    return colorplot_field(field, **kwargs)


def colorplot_field(field, xlim=(-np.pi, np.pi), ylim=(-np.pi, np.pi), resolution=100,
                    cmap='viridis', title="f(x, y)"):
    x = np.linspace(*xlim, resolution)
    y = np.linspace(*ylim, resolution)
    X, Y = np.meshgrid(x, y)

    Z = np.empty((resolution, resolution))
    for i in range(resolution):
        for j in range(resolution):
            Z[i, j] = field([X[i, j], Y[i, j]]).real

    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cmap)
    fig.colorbar(c, ax=ax, label="Re(f(x, y))")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    fig.tight_layout()

    return fig, ax
