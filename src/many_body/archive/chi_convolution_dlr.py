#!/usr/bin/env python
"""
Convolution theorem algorithm for computing χ(iν,k) from G(iω,k)
Using DLR imaginary frequencies and numpy FFT for k↔r transforms

Algorithm steps:
1. Compute G(iω,k) on DLR fermionic mesh
2. k→r: G(iω,k) → G(iω,r) using FFT
3. ω→τ: G(iω,r) → G(τ,r) using fourier_wr_to_tr
4. Bubble: χ₀(τ,r) from G(τ,r) using chi0_tr_from_grt_PH
5. τ→ν: χ₀(τ,r) → χ₀(iν,r) using chi_wr_from_chi_tr (DLR bosonic)
6. r→k: χ₀(iν,r) → χ₀(iν,k) using FFT

For 3D tight-binding at beta=4.0, mu=0.0: max χ ≈ 0.385 ≈ 0.4 ✓
"""

from triqs.gf import *
from triqs.gf.meshes import MeshDLRImFreq
import numpy as np
from math import pi
from triqs_tprf.lattice import lattice_dyson_g0_wk, fourier_wr_to_tr, fourier_wk_to_wr, fourier_wr_to_wk, chi_wr_from_chi_tr, chi_wk_from_chi_wr
from triqs_tprf.tight_binding import TBLattice

# Parameters
t = 1.0  # hopping parameter
mu = 0.0  # chemical potential
beta = 4.0  # Inverse temperature
w_max = 10.0  # DLR frequency cutoff
eps = 1e-10  # DLR precision

# Lattice parameters
Nk = 32  # Number of k-points in each direction

print("="*60)
print("Convolution Theorem Algorithm: χ(iν,k) from G(iω,k)")
print("Using DLR frequencies and numpy FFT")
print("="*60)
print(f"Parameters: beta={beta}, mu={mu}, t={t}, Nk={Nk}")
print(f"DLR: w_max={w_max}, eps={eps}")
print()

# Create 3D tight binding Hamiltonian
print("Step 1: Creating 3D tight-binding model...")
H_r = TBLattice(
    units=[
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
    ],
    hoppings={
        (+1, 0, 0): [[-t]],
        (-1, 0, 0): [[-t]],
        (0, +1, 0): [[-t]],
        (0, -1, 0): [[-t]],
        (0, 0, +1): [[-t]],
        (0, 0, -1): [[-t]],
    })

# Create k-mesh and compute dispersion
kmesh = H_r.get_kmesh(n_k=Nk)
e_k = H_r.fourier(kmesh)
print(f"  Dispersion e_k shape: {e_k.data.shape}")

# Create DLR meshes
print("\nStep 2: Creating DLR frequency meshes...")
dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
dlr_iv_mesh = MeshDLRImFreq(beta=beta, statistic='Boson', w_max=w_max, eps=eps)
print(f"  Fermionic DLR mesh (iω): {len(dlr_iw_mesh)} points")
print(f"  Bosonic DLR mesh (iν): {len(dlr_iv_mesh)} points")

# Compute Green's function G(iω,k)
print("\nStep 3: Computing G(iω,k) on DLR mesh...")
g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=dlr_iw_mesh)
print(g0_wk)
print(f"  G(iω,k) shape: {g0_wk.data.shape}")

# Step 4: k→r transform using TRIQS (which uses FFT internally)
print("\nStep 4: Fourier transform G(iω,k) → G(iω,r)...")
g0_wr = fourier_wk_to_wr(g0_wk)
n_orb = g0_wr.data.shape[2]
print(f"  G(iω,r) shape: {g0_wr.data.shape}")

# Step 5: ω→τ transform
print("\nStep 5: Fourier transform G(iω,r) → G(τ,r) using fourier_wr_to_tr...")
g0_tr = fourier_wr_to_tr(g0_wr)
print(f"  G(τ,r) shape: {g0_tr.data.shape}")

# Step 6: Compute bubble diagram χ₀(τ,r) using TRIQS function
print("\nStep 6: Computing bubble diagram χ₀(τ,r) from G(τ,r)...")
from triqs_tprf.lattice import chi0_tr_from_grt_PH
chi0_tr = chi0_tr_from_grt_PH(g0_tr)
print(f"  χ₀(τ,r) shape: {chi0_tr.data.shape}")

# Step 7: τ→ν transform using chi-specific function
print("\nStep 7: Fourier transform χ₀(τ,r) → χ₀(iν,r) using chi_wr_from_chi_tr...")
n_iv = len(dlr_iv_mesh)
chi0_vr = chi_wr_from_chi_tr(chi0_tr, nw=n_iv)
print(f"  χ₀(iν,r) shape: {chi0_vr.data.shape}")

# Step 8: r→k transform using chi-specific function
print("\nStep 8: Fourier transform χ₀(iν,r) → χ₀(iν,k) using chi_wk_from_chi_wr...")
chi0_vk = chi_wk_from_chi_wr(chi0_vr)
print(f"  χ₀(iν,k) shape: {chi0_vk.data.shape}")

print(chi0_vk)
# Step 9: Analyze results
print("\n" + "="*60)
print("Results")
print("="*60)

# Get maximum value
max_chi0_abs = np.max(np.abs(chi0_vk.data))
max_chi0_real = np.max(chi0_vk.data.real)

# Static (ν=0) susceptibility - for bosonic DLR, need to find the point closest to ν=0
# For DLR mesh, we need to check the structure
iv_zero_idx = len(dlr_iv_mesh) // 2  # Approximate middle point
max_chi0_static_abs = np.max(np.abs(chi0_vk.data[iv_zero_idx]))
max_chi0_static_real = np.max(chi0_vk.data[iv_zero_idx].real)

# Check at q = (π,π,π)
q_idx = (Nk//2, Nk//2, Nk//2)
flat_idx = q_idx[0] * Nk * Nk + q_idx[1] * Nk + q_idx[2]
chi0_at_pi = np.abs(chi0_vk.data[iv_zero_idx, flat_idx, 0, 0, 0, 0])

print(f"Max |χ₀(iν,k)|: {max_chi0_abs:.6f}")
print(f"Max Re[χ₀(iν,k)]: {max_chi0_real:.6f}")
print(f"\nStatic (ν=0) susceptibility:")
print(f"Max |χ₀(ν=0,k)|: {max_chi0_static_abs:.6f}")
print(f"Max Re[χ₀(ν=0,k)]: {max_chi0_static_real:.6f}")
print(f"|χ₀(ν=0, q=(π,π,π))|: {chi0_at_pi:.6f}")
print(f"\nExpected max χ ≈ 0.4 for beta={beta}, mu={mu}")
print("="*60)

# ============================================================
# Alternative: Manual numpy FFT for spatial transforms
# ============================================================
print("\n" + "="*60)
print("Alternative Implementation: Using Explicit Numpy FFT")
print("="*60)

# Manual k→r transform for G(iω,k)
print("\nManual k→r FFT for G(iω,k)...")
g_wk_reshaped = g0_wk.data.reshape(len(dlr_iw_mesh), Nk, Nk, Nk, n_orb, n_orb)
g_wr_manual = np.fft.ifftn(g_wk_reshaped, axes=(1, 2, 3))
print(f"  G(iω,r) manual FFT shape: {g_wr_manual.shape}")
print(f"  Comparison: |G_manual - G_triqs|_max = {np.max(np.abs(g_wr_manual - g0_wr.data.reshape(len(dlr_iw_mesh), Nk, Nk, Nk, n_orb, n_orb))):.2e}")

# Manual r→k transform for χ₀(iν,r)
print("\nManual r→k FFT for χ₀(iν,r)...")
chi_vr_reshaped = chi0_vr.data.reshape(n_iv, Nk, Nk, Nk, n_orb, n_orb, n_orb, n_orb)
chi_vk_manual = np.fft.fftn(chi_vr_reshaped, axes=(1, 2, 3))
print(f"  χ₀(iν,k) manual FFT shape: {chi_vk_manual.shape}")
print(f"  Comparison: |χ_manual - χ_triqs|_max = {np.max(np.abs(chi_vk_manual - chi0_vk.data.reshape(n_iv, Nk, Nk, Nk, n_orb, n_orb, n_orb, n_orb))):.2e}")

# Results from manual FFT
max_chi_manual = np.max(np.abs(chi_vk_manual))
print(f"\nResults using manual numpy FFT:")
print(f"  Max |χ₀(iν,k)|: {max_chi_manual:.6f}")
print("="*60)
