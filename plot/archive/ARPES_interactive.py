import firefly as fly
import numpy as np
from typing import List
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from collections import defaultdict
from skimage import measure

BZ =  [[1.63882615, 0.        , 0.        ],
 [0.         ,1.61542976 ,0.        ],
 [0.         ,0.         ,0.53855993]]
starting_band = 6
nbnd = 8
mu = 9.6851

SLICE_RESOLUTION = 40
# === CONFIGURATION ===
nz = 4
ntheta = 8
nphi = 8

def generate_plane_mesh(height: float, theta: float, phi: float, BZ: list, resolution: int):
    xp = np.linspace(-BZ[0][0]/2, BZ[0][0]/2, resolution)
    yp = np.linspace(-BZ[1][1]/2, BZ[1][1]/2, resolution)
    X, Y = np.meshgrid(xp, yp, indexing='xy')
    Z = np.full_like(X, height)
    pts = np.stack([X.flatten(), Y.flatten(), Z.flatten()], axis=1)  # shape (N², 3)
    R = rotation_matrix_from_theta_phi(theta, phi)  # aligns n → z-axis
    pts_rotated = pts @ R.T
    # Reshape back to meshgrid arrays
    X_rot = pts_rotated[:, 0].reshape((resolution, resolution))
    Y_rot = pts_rotated[:, 1].reshape((resolution, resolution))
    Z_rot = pts_rotated[:, 2].reshape((resolution, resolution))

    return X_rot, Y_rot, Z_rot
    #normal = np.array([
    #    np.sin(theta) * np.cos(phi),
    #    np.sin(theta) * np.sin(phi),
    #    np.cos(theta)
    #])
    #center = height * normal
    ##lin = np.linspace(-size / 2, size / 2, resolution)
    #xp = np.linspace(-BZ[0][0]/2, BZ[0][0]/2, resolution)
    #yp = np.linspace(-BZ[1][1]/2, BZ[1][1]/2, resolution)
    #xx, yy = np.meshgrid(xp, yp)
    #if np.allclose(normal, [0, 0, 1]):
    #    u = np.array([1, 0, 0])
    #else:
    #    u = np.cross(normal, [0, 0, 1])
    #    u /= np.linalg.norm(u)
    #v = np.cross(normal, u)
    #X = center[0] + xx * u[0] + yy * v[0]
    #Y = center[1] + xx * u[1] + yy * v[1]
    #Z = center[2] + xx * u[2] + yy * v[2]
    #return X, Y, Z

def rotation_matrix_from_theta_phi(theta, phi):
    # Original normal from spherical angles
    n = np.array([
        np.sin(phi) * np.cos(theta),
        np.sin(phi) * np.sin(theta),
        np.cos(phi)
    ])
    n = n / np.linalg.norm(n)

    # Target direction: z-axis
    k = np.array([0, 0, 1])

    # Cross product and angle between
    v = np.cross(n, k)
    s = np.linalg.norm(v)
    c = np.dot(n, k)

    if np.isclose(s, 0):  # already aligned
        return np.eye(3)

    vx = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

    R = np.eye(3) + vx + vx @ vx * ((1 - c) / s**2)
    return R

def generate_flat_plane(height, size, resolution):
    lin = np.linspace(-size / 2, size / 2, resolution)
    X, Y = np.meshgrid(lin, lin)
    Z = np.full_like(X, fill_value=height)
    return X, Y, Z

def fold(x, y, z, BZ):
    BZ = np.array(BZ)
    k_cart = [x, y, z]
    BZ_inv = np.linalg.inv(BZ.T)  # convert Cartesian to fractional
    k_frac = np.dot(BZ_inv, k_cart)

    # Fold to range [-0.5, 0.5)
    k_frac_folded = k_frac - np.round(k_frac)

    # Convert back to Cartesian
    k_folded = np.dot(BZ.T, k_frac_folded)

    return k_folded[0], k_folded[1], k_folded[2]

def generate_slice_surface(z, theta, phi):
    X, Y, Z = generate_plane_mesh(z, theta, phi, BZ, SLICE_RESOLUTION)
    bands = fly.Field_R("bands.dat")

    zero_points = []

    for n in range(starting_band, nbnd+1):  # Loop over bands
        values = np.zeros_like(X)

        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                x, y, z_folded = fold(X[i, j], Y[i, j], Z[i, j], BZ)
                values[i, j] = bands(n, [x, y, z_folded])

        contours = measure.find_contours(values, level=mu)

        for contour in contours:
            for i_float, j_float in contour:
                i_int, j_int = int(i_float), int(j_float)
                if 0 <= i_int < X.shape[0] and 0 <= j_int < X.shape[1]:
                    zero_points.append([
                        X[i_int, j_int],
                        Y[i_int, j_int],
                        Z[i_int, j_int]
                    ])

    return np.array(zero_points)

def generate_surfaces(nz, ntheta, nphi):
    points, thetas, phis, centers = [], [], [], []

    for i in range(nz):
        for j in range(ntheta):
            for k in range(nphi):
                z = get_z(i, nz)
                theta = get_theta(j, ntheta)
                phi = get_phi(k, nphi)
                name = f"theta_{theta}_phi_{phi}_z_{z}"
                zero_points = generate_slice_surface(z, theta, phi)
                thetas.append(theta)
                phis.append(phi)
                centers.append(z)
                points.append(zero_points)

    return points, thetas, phis, centers

def apply_rotation_to_mesh(X, Y, Z, R):
    # Flatten and stack into (N^2, 3) array of points
    pts = np.stack([X.flatten(), Y.flatten(), Z.flatten()], axis=1)  # shape (N², 3)

    # Apply rotation
    rotated = pts @ R.T  # shape (N², 3)

    # Reshape back to meshgrid shape
    X_rot = rotated[:, 0].reshape(X.shape)
    Y_rot = rotated[:, 1].reshape(Y.shape)
    Z_rot = rotated[:, 2].reshape(Z.shape)

    return X_rot, Y_rot, Z_rot

# --- Parametric mapping ---
def get_z(i, n): return round(0 + BZ[2][2] / 2 * i / n, 2)
def get_theta(i, n): return round(0 + 0.5 * np.pi * i / n, 2)
def get_phi(i, n): return round(0 + 0.5 * np.pi * i / n, 2)

# --- Build animated plane surface ---
def generate_planes(thetas, phis, centers):
    Xs, Ys, Zs, names = [], [], [], []
    Xf, Yf, Zf = [], [], []
    z_vals, frame_groups = set(), defaultdict(list)

    for theta, phi, z in zip(thetas, phis, centers):
        name = f"theta_{theta}_phi_{phi}_z_{z}"
        X_flat, Y_flat, Z_flat = generate_flat_plane(z, 2*np.pi, SLICE_RESOLUTION)
        R = rotation_matrix_from_theta_phi(theta, phi)
        X, Y, Z = apply_rotation_to_mesh(X_flat, Y_flat, Z_flat, R)
        Xf.append(X_flat)
        Yf.append(Y_flat)
        Zf.append(Z_flat)
        Xs.append(X)
        Ys.append(Y)
        Zs.append(Z)
        names.append(name)
        z_vals.add(z)
        frame_groups[z].append(name)

    return Xs, Ys, Zs, names, sorted(z_vals), frame_groups, Xf, Yf, Zf

# --- Generate Plotly surface figure + frames ---
def plotly_animated_surface(Xs, Ys, Zs, names):
    frames = [
        go.Frame(
            name=names[i],
            data=[go.Surface(
                z=Zs[i], x=Xs[i], y=Ys[i],
                colorscale=[[0, 'gray'], [1, 'gray']],
                showscale=False  # ✅ Fix here
            )]
        )
        for i in range(len(Zs))
    ]
    fig = go.Figure(data=[frames[0].data[0]], frames=frames)
    return fig


def plotly_animated_cut(points_2d: List[np.ndarray], names: List[str]) -> go.Figure:
    """
    Creates an animated 2D point plot where each frame shows a slice of zero-crossing points.
    Adds an invisible anchor trace to fix axis limits across all frames.
    """
    # Create frames
    frames = []
    for i, pts in enumerate(points_2d):
        if pts.size == 0:
            trace = go.Scatter(x=[], y=[], mode="markers",
                               marker=dict(size=12, color='#e5690e', symbol='circle'),
                               name="cut")
        else:
            trace = go.Scatter(x=pts[:, 0], y=pts[:, 1], mode="markers",
                               marker=dict(size=12, color='#e5690e', symbol='circle'),
                               name="cut")
        
        # Invisible anchor to fix axes
        anchor = go.Scatter(
            x=[-1, 1], y=[-1, 1],
            mode="markers",
            marker=dict(opacity=0),
            showlegend=False
        )

        frames.append(go.Frame(name=names[i], data=[trace, anchor]))

    # Initial frame
    if points_2d[0].size == 0:
        initial_trace = go.Scatter(x=[], y=[], mode="markers",
                                   marker=dict(size=12, color='#e5690e', symbol='circle'),
                                   name="cut")
    else:
        initial_trace = go.Scatter(x=points_2d[0][:, 0], y=points_2d[0][:, 1], mode="markers",
                                   marker=dict(size=12, color='#e5690e', symbol='circle'),
                                   name="cut")

    # Initial invisible anchor
    initial_anchor = go.Scatter(
        x=[-1, 1], y=[-1, 1],
        mode="markers",
        marker=dict(opacity=0),
        showlegend=False
    )

    # Construct figure
    fig = go.Figure(data=[initial_trace, initial_anchor], frames=frames)

    return fig

# --- Group static + animated traces per frame ---
def get_frames(static_traces, plane_frames, slice_frames):
    frames = []
    for i in range(len(plane_frames)):
        plane = plane_frames[i].data[0]
        slice_ = slice_frames[i].data[0]
        name = plane_frames[i].name  # assuming slice_frames[i].name == plane_frames[i].name
        frames.append(go.Frame(
            name=name,
            data=static_traces + [plane, slice_]
        ))
    return frames

def make_flat(points: list[np.ndarray], thetas: list[float], phis: list[float], centers: list[np.ndarray]) -> list[np.ndarray]:
    flat = []

    for p, theta, phi, center in zip(points, thetas, phis, centers):
        # Rebuild rotation matrix from theta, phi
        R = rotation_matrix_from_theta_phi(theta, phi)  # aligns n → z-axis

        # Translate points to rotate around origin
        p_centered = p - center

        # Apply inverse rotation (rotate plane → flat)
        p_flat = p_centered @ R.T  # R.T = R⁻¹

        # Take x, y only
        flat.append(p_flat[:, :2])

    return flat


# === Generate surface frames ===
points, thetas, phis, centers = generate_surfaces(nz, ntheta, nphi)
Xs, Ys, Zs, frame_names, z_vals, frames_by_z, Xf, Yf, Zf = generate_planes(thetas, phis, centers)
flat_points = make_flat(points, thetas, phis, centers)
#slice_fig = plotly_animated_pointcloud(flat_points, frame_names)
slice_fig = plotly_animated_cut(flat_points, frame_names)
plane_fig = plotly_animated_surface(Xs, Ys, Zs, frame_names)
plane_frames = plane_fig.frames
slice_frames = slice_fig.frames

bands = fly.Field_R("bands.dat")

# === Generate isosurface data ===
pi = np.pi
x, y, z = np.mgrid[-BZ[0][0]/2:BZ[0][0]:20j, -BZ[1][1]/2:BZ[1][1]:20j, -BZ[2][2]/2:BZ[2][2]:20j] 

isosurface_traces = []
for n in range(starting_band, nbnd+1):
    f = np.zeros_like(x)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            for k in range(x.shape[2]):
                x_v, y_v, z_v = fold(x[i, j, k], y[i, j, k], z[i, j, k], BZ)
                f[i, j, k] = bands(n, [x_v, y_v, z_v])

    trace = go.Isosurface(
        x=x.flatten(), y=y.flatten(), z=z.flatten(), value=f.flatten(),
        isomin=mu, isomax=mu,
        surface_count=1,
        showscale=False,
        colorscale=[[0, '#e5800e'], [1, '#e5800e']],
        opacity=0.8,
        caps=dict(x_show=False, y_show=False, z_show=False)
    )
    isosurface_traces.append(trace)

# === Create figure ===
fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scene'}, {'type': 'xy'}]],
    column_widths=[0.6, 0.4]
)
for trace in isosurface_traces:
    fig.add_trace(trace, row=1, col=1)

#fig.add_trace(go.Scatter(
#    x=[-1, 1], y=[-1, 1],
#    mode='markers',
#    marker=dict(opacity=0),
#    showlegend=False
#), row=1, col=2)
fig.add_trace(plane_frames[0].data[0], row=1, col=1)
fig.add_trace(slice_frames[0].data[0], row=1, col=2)
fig.frames = get_frames(isosurface_traces, plane_frames, slice_frames)

# === Build z-value dropdown to trigger animation ===
dropdown_buttons = [
    {
        "label": f"z = {z}",
        "method": "animate",
        "args": [frames_by_z[z], {
            "frame": {"duration": 0, "redraw": True},
            "mode": "immediate",
            "transition": {"duration": 0}
        }]
    }
    for z in z_vals
]

# === Layout with clean dropdown and controls ===
fig.update_layout(
    height=800,
    width=1500,
    scene=dict(
        xaxis=dict(range=[-BZ[0][0]/2, BZ[0][0]/2], autorange=False),
        yaxis=dict(range=[-BZ[1][1]/2, BZ[1][1]/2], autorange=False),
        zaxis=dict(range=[-BZ[2][2]/2, BZ[2][2]/2], autorange=False),
        aspectmode='manual',
        aspectratio=dict(x=1, y=1, z=1)
    ),
    xaxis2=dict(
        range=[-1.0,1.0],
        autorange=False,
        fixedrange=True,
        scaleanchor='y2',      # 🔒 lock aspect ratio to y2
        scaleratio=1           # 1:1 aspect
    ),
    yaxis2=dict(
        range=[-1.0,1.0],
        autorange=False,
        fixedrange=True,
        scaleanchor='y2',      # 🔒 lock aspect ratio to y2
        scaleratio=1           # 1:1 aspect
    ),
updatemenus=[
    {
        "type": "buttons",
        "buttons": [
            {
                "label": "Play",
                "method": "animate",
                "args": [None, {
                    "frame": {"duration": 800, "redraw": True},  # ⬅️ smoother pacing
                    "fromcurrent": True,
                    "transition": {"duration": 500},
                    "mode": "afterall"  # ✅ wait for durations
                }]
            },
            {
                "label": "Pause",
                "method": "animate",
                "args": [[None], {
                    "mode": "immediate",
                    "frame": {"duration": 0},
                    "transition": {"duration": 0}
                }]
            }
        ],
        "pad": {"r": 10, "t": 10}
    },
    {
        "buttons": [
            {
                "args": [[name], {
                    "frame": {"duration": 800, "redraw": True},  # ⬅️ dropdown transitions
                    "mode": "afterall",
                    "transition": {"duration": 0}
                }],
                "label": name,
                "method": "animate"
            } for name in frame_names
        ],
        "direction": "down",
        "showactive": True,
        "x": 0.1,
        "xanchor": "left",
        "y": 1.15,
        "yanchor": "top"
    }
]
)

# === Export ===
fig.write_html("interactive.html", include_plotlyjs='cdn', full_html=True, auto_play=False)

#h, p, t = 0, 1, 0               # Center, phi, theta
#X, Y, Z = generate_plane_mesh(h, p, t, BZ, 3)
#points = generate_slice_surface(h, p, t)
#print("X, Y, Z: ", X, Y, Z)
#print("points: ", points)
#R = rotation_matrix_from_theta_phi(t, p)  # aligns n → z-axis
#print("R: ", R)
#pts = np.stack([X.flatten(), Y.flatten(), Z.flatten()], axis=1)  # shape (N², 3)
#print("pts: ", pts)
#point_centered = pts - h
#point_flat = point_centered @ R  # R.T = R⁻¹
#print("point_flat", point_flat)
