import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N = 400
x = np.linspace(-2, 2, 400)
y = np.linspace(-2, 2, 400)
X, Y = np.meshgrid(x, y)

Sigma = ffirefly.self_energy(sys.argv[1])
BZ = Sigma.cmf.domain
mesh_r = np.array(N, N, N)
mesh_i = np.array(N, N, N)

for i in range(N):
    for j in range(N):
        for k in range(N):
            q = BZ * [i/N - 0.5, j/N - 0.5, k/N - 0.5]
            mesh_r[i][j][k] = real(Sigma(q))
            mesh_i[i][j][k] = imag(Sigma(q))

# Imaginary part for color strength (closer to 0 = brighter)
imag_magnitude = np.abs(mesh_i)

# Normalize to [0, 1] and invert (low imaginary = high intensity)
color_strength = 1 - (imag_magnitude / np.max(imag_magnitude))
colors = plt.cm.viridis(color_strength)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(
    X, Y, mesh_r,
    facecolors=colors,
    rstride=1,
    cstride=1,
    linewidth=0,
    antialiased=False,
    shade=False
)

# Hide the default color bar (or show your own with custom logic)
plt.show()

