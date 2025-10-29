import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from typing import Optional
import matplotlib.tri as tri
from scipy.spatial import ConvexHull
from matplotlib.collections import LineCollection
import firefly as fly
import firefly.config as cfg

def plot_colorgrid(files, ax: Optional[Axes] = None, **kwargs) -> Axes:
    field = fly.Field_C(files[0])
    return colorplot_field(field, ax=ax, **kwargs)

def plot_colorline(files, ax: Optional[Axes] = None, **kwargs) -> Axes:
    return colorplot_surface(files[0], ax=ax, **kwargs)

def colorplot_field(field, ax: Optional[Axes] = None, xlim=(-np.pi, np.pi), ylim=(-np.pi, np.pi),
                    resolution=500, cmap='bwr', title="f(x, y)") -> Axes:
    x = np.linspace(*xlim, resolution)
    y = np.linspace(*ylim, resolution)
    X, Y = np.meshgrid(x, y)

    maxZ = -1000
    minZ = 1000
    Z = np.empty((resolution, resolution))
    for i in range(resolution):
        for j in range(resolution):
            Z[i, j] = field([X[i, j], Y[i, j]]).real
            maxZ = max(maxZ, Z[i, j])
            minZ = min(minZ, Z[i, j])

    if (maxZ * minZ > 0):
        cmap = "viridis"

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    c = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cmap)
    fig.colorbar(c, ax=ax, label="Re(f(x, y))")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)

    if ax is None or ax.get_figure().get_axes() == [ax]:
        fig.tight_layout()

    return ax


def colorplot_surface(file, ax: Optional[Axes] = None, cmap='bwr') -> Axes:
    df = pd.read_csv(file, sep=None, engine='python')
    x = df.iloc[:, 0].values
    y = df.iloc[:, 1].values
    f = df.iloc[:, 2].values

    points = np.column_stack((x, y))
    hull = ConvexHull(points)
    ordered_idx = hull.vertices
    ordered_points = points[ordered_idx]
    ordered_f_vals = f[ordered_idx]

    ordered_points = np.vstack([ordered_points, ordered_points[0]])
    ordered_f_vals = np.append(ordered_f_vals, ordered_f_vals[0])

    segments = np.array([
        [ordered_points[i], ordered_points[i + 1]]
        for i in range(len(ordered_points) - 1)
    ])
    segment_vals = 0.5 * (ordered_f_vals[:-1] + ordered_f_vals[1:])

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.get_figure()

    lc = LineCollection(segments, cmap=cmap, array=segment_vals, linewidths=5)
    ax.add_collection(lc)
    fig.colorbar(lc, ax=ax)

    # ✅ Fix: Set axis limits to match the data
    ax.set_xlim(ordered_points[:, 0].min() * 1.03, ordered_points[:, 0].max() * 1.03)
    ax.set_ylim(ordered_points[:, 1].min() * 1.03, ordered_points[:, 1].max() * 1.03)

    ax.set_aspect('equal')

    if ax is None or ax.get_figure().get_axes() == [ax]:
        fig.tight_layout()

    return ax
