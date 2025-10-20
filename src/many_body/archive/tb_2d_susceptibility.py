from triqs.gf import *
import numpy as np
from math import pi
import matplotlib.pyplot as plt
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.tight_binding import TBLattice

print("1")
# Parameters
t = 1.0  # hopping parameter
mu = 0.0  # chemical potential
beta = 4.0
w_max = 10.0
eps = 1e-10

print("2")
H_r = TBLattice(
    units=[
        (1,0,0), # basis vector in the x-direction
        (0,1,0), # basis vector in the y-direction
    ],
    hoppings={
        (+1,0) : [[-t]], # hopping in the +x direction
        (-1,0) : [[-t]], # hopping in the -x direction
        (0,+1) : [[-t]], # hopping in the +y direction
        (0,-1) : [[-t]], # hopping in the -y direction
    })
print("3")
# Lattice parameters
Nk = 60
Nx, Ny = 60, 60  # Number of k-points in each direction
Lx, Ly = Nx, Ny  # Real space lattice size (same as k-grid for clean FFTs)

def epsilon_k(kx, ky):
    """2D tight binding dispersion on square lattice: ε(k) = -2t(cos(kx) + cos(ky))"""
    return -2.0 * t * (np.cos(kx) + np.cos(ky))

# Create k-space mesh
kx_vals = np.linspace(-pi, pi, Nx, endpoint=False)
ky_vals = np.linspace(-pi, pi, Ny, endpoint=False)
kx_grid, ky_grid = np.meshgrid(kx_vals, ky_vals, indexing='ij')

# Create DLR frequency mesh

print("4")
# Griffin code
dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
print("4.1")
kmesh = H_r.get_kmesh(n_k=Nk)
print("4.2")
e_k = H_r.fourier(kmesh)
print("4.3")
wkmesh = MeshProduct(dlr_iw_mesh, kmesh)
print("5")
g0_wk = lattice_dyson_g0_wk(mu=0., e_k=e_k, mesh=dlr_iw_mesh)
print("6")
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
print("7")
chi0_wk = 2 * imtime_bubble_chi0_wk(g0_wk, nw=100) # Factor of 2 for spin
print(chi0_wk)
start = int(len(chi0_wk.data) / 2)
print(len(chi0_wk.data))
print(len(chi0_wk.data[0]))
print(len(chi0_wk.data[0][0]))
print(chi0_wk.data.shape)
chidata = np.reshape(chi0_wk.data[start], (Nk, Nk))
print("8")
k = np.linspace(0, 2*np.pi, num=100, endpoint=True)
kx, ky = np.meshgrid(k, k)

#chi_interp = np.vectorize(lambda kx, ky: chi0_wk(0, (kx, ky, 0)).real)

plt.pcolormesh(chidata.T.real)
#plt.pcolormesh(kx, ky, chi_interp(kx, ky), rasterized=True)

plt.title('Static susceptibility $\chi_0(\mathbf{q}, \omega=0)$')
ticks, labels = [0, np.pi, 2*np.pi], [r"0",r"$\pi$",r"$2\pi$"]
plt.xticks(ticks, labels); plt.yticks(ticks, labels);
plt.xlabel(r'$q_x$'); plt.ylabel(r'$q_y$')
plt.colorbar();
plt.show()
