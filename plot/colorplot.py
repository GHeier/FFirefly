import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.spatial import ConvexHull
from matplotlib.collections import LineCollection
import firefly as fly
import firefly.config as cfg

def plot_colorgrid(files, **kwargs):
    field = fly.Field_C(files[0])
    return colorplot_field(field, **kwargs)

def plot_colorline(files, **kwargs):
    return colorplot_surface(files[0], **kwargs)

def colorplot_field(field, xlim=(-np.pi, np.pi), ylim=(-np.pi, np.pi), resolution=100,
                    cmap='bwr', title="f(x, y)"):
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
    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.pcolormesh(X, Y, Z, shading='auto', cmap=cmap)
    fig.colorbar(c, ax=ax, label="Re(f(x, y))")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    fig.tight_layout()

    return fig, ax


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.spatial import ConvexHull

def colorplot_surface(file, cmap='bwr'):
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

    fig, ax = plt.subplots(figsize=(8, 6))
    lc = LineCollection(segments, cmap=cmap, array=segment_vals, linewidths=5)
    ax.add_collection(lc)
    fig.colorbar(lc, ax=ax, label='f(x, y)')

    # ✅ Fix: Set axis limits to match the data
    ax.set_xlim(ordered_points[:, 0].min(), ordered_points[:, 0].max())
    ax.set_ylim(ordered_points[:, 1].min(), ordered_points[:, 1].max())

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('2D Surface Colored by f(x, y)')
    ax.set_aspect('equal')
    fig.tight_layout()

    return fig, ax
