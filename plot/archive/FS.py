import numpy as np
import pyvista as pv
import firefly

# Color palette
colors = ['purple', 'blue', 'orange', 'red', 'green']
hex_colors = ['#b106fc', '#067dfc', '#fc8806', '#fc2006', '#19fc06']

# Tight-binding parameters
t = 1.0
a = 1.0
k_max = np.pi  # Brillouin zone extent
n = 100  # Grid resolution

# List of Fermi levels to plot
fermi_levels = [-4.0, -2.0, -1.0, 0.0, 1.0]  # Adjust or expand this list as needed

# Generate k-space grid
kx = np.linspace(-k_max, k_max, n)
ky = np.linspace(-k_max, k_max, n)
kz = np.linspace(-k_max, k_max, n)
kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing='ij')

# 3D tight-binding dispersion relation
E = -2 * t * (np.cos(kx3d * a) + np.cos(ky3d * a) + np.cos(kz3d * a))

# Create PyVista grid using ImageData (vtkImageData)
grid = pv.ImageData()
grid.dimensions = E.shape
grid.origin = (kx[0], ky[0], kz[0])
grid.spacing = (
    (kx[-1] - kx[0]) / (n - 1),
    (ky[-1] - ky[0]) / (n - 1),
    (kz[-1] - kz[0]) / (n - 1)
)
grid.point_data["energy"] = E.flatten(order="F")

# Start PyVista plotter
plotter = pv.Plotter()

# Plot a contour surface for each Fermi level
for idx, E_F in enumerate(fermi_levels):
    contours = grid.contour([E_F], scalars="energy")
    color = hex_colors[idx % len(hex_colors)]  # Wrap color index if too many levels
    plotter.add_mesh(contours, color=color, opacity=0.65, label=f"E_F = {E_F:.2f}")

# Add legend and axes
plotter.add_axes()
plotter.add_legend()
plotter.show()
