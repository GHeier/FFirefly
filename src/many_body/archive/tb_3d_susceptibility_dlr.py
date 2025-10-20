from triqs.gf import *
import numpy as np
from math import pi
from triqs_tprf.lattice import *
from triqs_tprf.tight_binding import TBLattice

# Parameters
t = 1.0  # hopping parameter
mu = 0.0  # chemical potential
beta = 4.0  # Inverse temperature
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

print("Creating imaginary frequency mesh...")
# NOTE: DLR mesh causes segfault in this TRIQS version, using Matsubara instead
n_iw = 100
iw_mesh = MeshImFreq(beta=beta, S='Fermion', n_iw=n_iw)
dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
dlr_iv_mesh = MeshDLRImFreq(beta=beta, statistic='Boson', w_max=w_max, eps=eps)

print("Computing Green's function G(iw,k)...")
g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=dlr_iw_mesh)
print(f"G(iw,k) shape: {g0_wk.data.shape}")

print("Step 1: Fourier transform G(iw,k) -> G(iw,r)...")
g0_wr = fourier_wk_to_wr(g0_wk)
print(f"G(iw,r) shape: {g0_wr.data.shape}")

print("Step 2: Fourier transform G(iw,r) -> G(tau,r)...")
g0_tr = fourier_wr_to_tr(g0_wr)
print(g0_tr)
print(f"G(tau,r) shape: {g0_tr.data.shape}")

print("Step 3: Computing chi0(tau,r) using convolution theorem (bubble diagram)...")
chi0_tr = chi0_tr_from_grt_PH(g0_tr)
print(chi0_tr)
print(f"chi0(tau,r) shape: {chi0_tr.data.shape}")

print("Step 4: Fourier transform chi0(tau,r) -> chi0(iw,r)...")
chi0_wr = fourier_tr_to_wr(chi0_tr)
print(chi0_wr)
#chi0_wr = chi_wr_from_chi_tr(chi0_tr, nw=n_iw)
print(f"chi0(iw,r) shape: {chi0_wr.data.shape}")

print("Step 5: Fourier transform chi0(iw,r) -> chi0(iw,k)...")
chi0_wk = fourier_wr_to_wk(chi0_wr)
#chi0_wk = chi_wk_from_chi_wr(chi0_wr)
print(f"chi0(iw,k) shape: {chi0_wk.data.shape}")

# Get maximum value of chi0
max_chi0_abs = np.max(np.abs(chi0_wk.data))
max_chi0_real = np.max(chi0_wk.data.real)

# Also check the static (w=0) susceptibility
start = int(len(chi0_wk.data) / 2)
max_chi0_static_abs = np.max(np.abs(chi0_wk.data[start]))
max_chi0_static_real = np.max(chi0_wk.data[start].real)

# Check specific momentum points
# Get chi0 at q = (pi, pi, pi) - the zone boundary
q_idx = (Nk//2, Nk//2, Nk//2)
flat_idx = q_idx[0] * Nk * Nk + q_idx[1] * Nk + q_idx[2]
chi0_at_pi = np.abs(chi0_wk.data[start, flat_idx, 0, 0, 0, 0])

# Print results
print(f"\n" + "="*60)
print(f"Results for 3D tight binding at beta={beta}, t={t}, mu={mu}")
print(f"Using explicit Fourier transforms (fourier_wk_to_wr, fourier_wr_to_tr)")
print(f"="*60)
print(f"Max |chi0(iw,k)|: {max_chi0_abs:.6f}")
print(f"Max Re[chi0(iw,k)]: {max_chi0_real:.6f}")
print(f"\nStatic (w=0) susceptibility:")
print(f"Max |chi0(w=0,k)|: {max_chi0_static_abs:.6f}")
print(f"Max Re[chi0(w=0,k)]: {max_chi0_static_real:.6f}")
print(f"|chi0(w=0, q=(π,π,π))|: {chi0_at_pi:.6f}")
print(f"="*60)
