from triqs.gf import *
import numpy as np
from math import pi
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.tight_binding import TBLattice
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

# Parameters
t = 1.0  # hopping parameter
mu = 0.0  # chemical potential
beta = 4.0
w_max = 10.0
eps = 1e-10

# Create 3D tight binding Hamiltonian
H_r = TBLattice(
    units=[
        (1, 0, 0),  # basis vector in the x-direction
        (0, 1, 0),  # basis vector in the y-direction
        (0, 0, 1),  # basis vector in the z-direction
    ],
    hoppings={
        (+1, 0, 0): [[-t]],  # hopping in the +x direction
        (-1, 0, 0): [[-t]],  # hopping in the -x direction
        (0, +1, 0): [[-t]],  # hopping in the +y direction
        (0, -1, 0): [[-t]],  # hopping in the -y direction
        (0, 0, +1): [[-t]],  # hopping in the +z direction
        (0, 0, -1): [[-t]],  # hopping in the -z direction
    })

# Lattice parameters
Nk = 32  # Number of k-points in each direction

print("Creating k-mesh...")
kmesh = H_r.get_kmesh(n_k=Nk)

print("Computing dispersion...")
e_k = H_r.fourier(kmesh)

print("Creating Matsubara frequency mesh...")
n_iw = 100
iw_mesh = MeshImFreq(beta=beta, S='Fermion', n_iw=n_iw)

print("Computing Green's function G(iw,k)...")
g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=iw_mesh)

print("Computing susceptibility chi0(iw,k) using convolution theorem...")
# Computing bare susceptibility (without explicit spin factor)
chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=100)

print(f"chi0_wk shape: {chi0_wk.data.shape}")

# Get maximum value of chi0
max_chi0_abs = np.max(np.abs(chi0_wk.data))
max_chi0_real = np.max(chi0_wk.data.real)

# Also check the static (w=0) susceptibility
start = int(len(chi0_wk.data) / 2)
max_chi0_static_abs = np.max(np.abs(chi0_wk.data[start]))
max_chi0_static_real = np.max(chi0_wk.data[start].real)

# Find momentum point of maximum
max_idx = np.unravel_index(np.argmax(np.abs(chi0_wk.data[start])), chi0_wk.data[start].shape)

# Check specific momentum points
# Get chi0 at q = (pi, pi, pi) - the zone boundary
# In the k-mesh with Nk points, pi corresponds to index Nk/2
q_idx = (Nk//2, Nk//2, Nk//2)
# Map 3D index to flat index for the k-mesh
flat_idx = q_idx[0] * Nk * Nk + q_idx[1] * Nk + q_idx[2]
chi0_at_pi = np.abs(chi0_wk.data[start, flat_idx, 0, 0, 0, 0])

# Print with different potential normalizations
print(f"\nResults for 3D tight binding at beta={beta}, t={t}, mu={mu}:")
print(f"="*60)
print(f"Max |chi0(iw,k)|: {max_chi0_abs:.6f}")
print(f"Max Re[chi0(iw,k)]: {max_chi0_real:.6f}")
print(f"\nStatic (w=0) susceptibility:")
print(f"Max |chi0(w=0,k)|: {max_chi0_static_abs:.6f}")
print(f"Max Re[chi0(w=0,k)]: {max_chi0_static_real:.6f}")
print(f"|chi0(w=0, q=(π,π,π))|: {chi0_at_pi:.6f}")
print(f"\nWith spin factor of 2:")
print(f"Max |chi0(w=0,k)| * 2: {max_chi0_static_abs * 2:.6f}")
print(f"="*60)
