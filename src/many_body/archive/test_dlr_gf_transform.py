from triqs.gf import *
from triqs_tprf.lattice import *
from triqs_tprf.tight_binding import TBLattice
import numpy as np

# Parameters
t = 1.0
mu = 0.0
beta = 10.0
w_max = 10.0
eps = 1e-10

# Create 2D tight binding model
H_r = TBLattice(
    units=[(1, 0), (0, 1)],
    hoppings={
        (+1, 0): [[-t]],
        (-1, 0): [[-t]],
        (0, +1): [[-t]],
        (0, -1): [[-t]],
    })

# Create k-mesh
Nk = 16
print(f"Creating {Nk}x{Nk} k-mesh...")
kmesh = H_r.get_kmesh(n_k=Nk)

# Compute dispersion e(k)
print("Computing dispersion...")
e_k = H_r.fourier(kmesh)

# Create DLR imaginary frequency mesh
print(f"Creating DLR imfreq mesh (beta={beta}, w_max={w_max}, eps={eps})...")
dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
print(f"DLR mesh size: {len(dlr_iw_mesh)}")

# Step 1: Create G(iw,k) - vectorized over all k-points
print("\nStep 1: Computing G(iw,k)...")
g_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=dlr_iw_mesh)
print(f"  G(iw,k) shape: {g_wk.data.shape}")
print(f"  G(iw,k) mesh: {g_wk.mesh}")

# Step 2: Transform G(iw,k) -> G(iw,r) - vectorized k->r transform
print("\nStep 2: Fourier transform G(iw,k) -> G(iw,r)...")
g_wr = fourier_wk_to_wr(g_wk)
print(f"  G(iw,r) shape: {g_wr.data.shape}")

# Step 3: Transform G(iw,r) -> G(tau,r) - vectorized iw->tau transform
print("\nStep 3: Fourier transform G(iw,r) -> G(tau,r)...")
g_tr = fourier_wr_to_tr(g_wr)
print(f"  G(tau,r) shape: {g_tr.data.shape}")
print(f"  G(tau,r) mesh: {g_tr.mesh}")
print(g_tr)

# Verify results
print("\n" + "="*60)
print("Verification:")
print("="*60)

# Check G(tau=0^-, r=0) = -1 (for non-interacting system at half-filling)
r_zero_idx = Nk**2 // 2  # center point in real space
g_tau0_r0 = g_tr.data[-1, r_zero_idx, 0, 0]
print(f"G(tau=beta-, r=0) = {g_tau0_r0:.6f} (should be ~ -1 for half-filled system)")

# Check that G is local in time (Green's function properties)
print(f"\nMax |G(tau,r)|: {np.max(np.abs(g_tr.data)):.6f}")
print(f"Min Re[G(tau,r)]: {np.min(g_tr.data.real):.6f}")
print(f"Max Re[G(tau,r)]: {np.max(g_tr.data.real):.6f}")

print("\n" + "="*60)
print("SUCCESS: All transforms completed without looping over k-points!")
print("="*60)
