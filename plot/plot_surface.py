import ffireefly
import sys
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np

surf = ffirefly.Surface(sys.argv[1])

# Example: Load vertices and faces
vertices = surf.vertices
faces = surf.faces

# Build the 3D polygon collection
mesh = Poly3DCollection(vertices[faces], alpha=0.7)
mesh.set_facecolor('lightblue')
mesh.set_edgecolor('k')

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.add_collection3d(mesh)

# Auto scale to the mesh size
ax.auto_scale_xyz(vertices[:, 0], vertices[:, 1], vertices[:, 2])

plt.show()
