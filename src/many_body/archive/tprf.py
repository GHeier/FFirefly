import matplotlib.pyplot as plt
import numpy as np
from triqs_tprf.tight_binding import TBLattice

t = 1.0
H = TBLattice(
    units = [(1, 0, 0), (0, 1, 0)],
    hopping = {
        # nearest neighbour hopping -t
        ( 0,+1): -t * np.eye(2),
        ( 0,-1): -t * np.eye(2),
        (+1, 0): -t * np.eye(2),
        (-1, 0): -t * np.eye(2),
        },
    orbital_positions = [(0,0,0)]*2,
    orbital_names = ['up', 'do'],
    )
kmesh = H.get_kmesh(n_k=(32, 32, 1))

e_k = H.fourier(kmesh)

G = [0.0, 0.0, 0.0]
X = [0.5, 0.0, 0.0]
M = [0.5, 0.5, 0.0]

paths = [(G, X), (X, M), (M, G)]

from triqs_tprf.lattice_utils import k_space_path
k_vecs, k_plot, k_ticks = k_space_path(paths, bz=H.bz)

e_k_interp = np.vectorize(lambda k : e_k(k)[0, 0].real, signature='(n)->()')

plt.plot(k_plot, e_k_interp(k_vecs))
plt.xticks(k_ticks, [r'$\Gamma$',r'$X$',r'$M$',r'$\Gamma$'])
plt.ylabel(r'$\epsilon(\mathbf{k})$')
plt.grid(True)
e_k = H.fourier(H.get_kmesh(n_k=(32, 32, 1)))
k = np.linspace(-0.5, 0.5, num=200) * 2. * np.pi
Kx, Ky = np.meshgrid(k, k)

e_k_interp = np.vectorize(lambda kx, ky : e_k([kx, ky, 0])[0, 0].real)
e_k_interp = e_k_interp(Kx, Ky)

plt.figure()
extent = (k.min(), k.max(), k.min(), k.max())
plt.imshow(e_k_interp, cmap=plt.get_cmap('RdBu'),
           extent=extent, origin='lower')
plt.colorbar()
plt.contour(Kx, Ky, e_k_interp, levels=[0])
plt.xlabel(r'$k_x$'); plt.ylabel(r'$k_y$');
