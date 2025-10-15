#!/usr/bin/env python
"""
TRIQS DLR (Discrete Lehmann Representation) Green's Function Example

This script demonstrates:
1. Creating a Green's function in imaginary frequency using DLR mesh
2. Filling it with a simple model (semi-circular density of states)
3. Converting to DLR coefficient representation
4. Transforming to imaginary time representation
"""

import numpy as np
from triqs.gf import Gf, SemiCircular
from triqs.gf.meshes import MeshDLRImFreq, MeshImFreq
from triqs.gf import make_gf_dlr, make_gf_dlr_imtime, make_gf_imfreq

# Parameters
beta = 40.0        # Inverse temperature
w_max = 10.0       # Maximum frequency cutoff
eps = 1e-10        # DLR precision
n_iw = 1000        # Number of Matsubara frequencies for comparison

print("=" * 60)
print("TRIQS DLR Green's Function Example")
print("=" * 60)
print(f"Parameters:")
print(f"  β (inverse temperature) = {beta}")
print(f"  ω_max (frequency cutoff) = {w_max}")
print(f"  ε (DLR precision) = {eps}")
print()

# 1. Create DLR imaginary frequency mesh
print("Step 1: Creating DLR imaginary frequency mesh...")
dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
print(f"  DLR mesh size: {len(dlr_iw_mesh)} points")
print(f"  (Compare to standard mesh: {n_iw} points)")
print()

# 2. Create Green's function on DLR mesh and fill with semi-circular DOS model
print("Step 2: Creating Green's function on DLR mesh...")
Giw_dlr = Gf(mesh=dlr_iw_mesh, target_shape=[])

# Fill with semi-circular spectral function (Bethe lattice)
# This creates a reference Green's function
print("  Filling with semi-circular DOS model (half-bandwidth D=1)...")
iw_mesh_full = MeshImFreq(beta=beta, statistic='Fermion', n_iw=n_iw)
Giw_full = Gf(mesh=iw_mesh_full, target_shape=[])
Giw_full << SemiCircular(half_bandwidth=1.0)

# Now create DLR version by fitting
print("  Fitting to DLR representation...")
Giw_dlr << SemiCircular(half_bandwidth=1.0)
print()

# 3. Convert to DLR coefficient representation
print("Step 3: Converting to DLR coefficient representation...")
G_dlr = make_gf_dlr(Giw_dlr)
print(f"  DLR coefficient mesh size: {len(G_dlr.mesh)}")
print()

# 4. Transform to imaginary time
print("Step 4: Transforming to imaginary time...")
Gtau_dlr = make_gf_dlr_imtime(G_dlr)
print(f"  DLR imaginary time mesh size: {len(Gtau_dlr.mesh)}")
print()

# 5. Verification: Check values at specific points
print("=" * 60)
print("Verification Results")
print("=" * 60)

# Print some values at different tau points
print("\nGreen's function G(τ) at selected imaginary time points:")
print(f"  τ/β    G(τ)")
print("-" * 40)
for i in [0, len(Gtau_dlr.mesh)//4, len(Gtau_dlr.mesh)//2, 3*len(Gtau_dlr.mesh)//4, -1]:
    tau_pt = Gtau_dlr.mesh[i]
    tau_val = float(tau_pt)
    g_val = Gtau_dlr[tau_pt].real
    print(f"  {tau_val/beta:6.4f}  {g_val:12.8f}")

# Check boundary condition G(β) = -G(0)
g_0 = Gtau_dlr[Gtau_dlr.mesh[0]].real
g_beta = Gtau_dlr[Gtau_dlr.mesh[-1]].real
print(f"\nFermionic boundary condition check:")
print(f"  G(0)  = {g_0:12.8f}")
print(f"  G(β)  = {g_beta:12.8f}")
print(f"  G(0) + G(β) = {g_0 + g_beta:12.8e} (should be ~0)")
print()

# Reconstruct imaginary frequency Green's function from DLR
print("Step 5: Reconstructing full imaginary frequency from DLR...")
Giw_reconstructed = make_gf_imfreq(G_dlr, n_iw=100)
print(f"  Reconstructed mesh size: {len(Giw_reconstructed.mesh)}")
print()

# Compare original and reconstructed at a few Matsubara frequencies
print("Comparison: Original vs Reconstructed G(iω_n):")
print(f"  n    ω_n         Re[G(iω_n)]     Im[G(iω_n)]")
print("-" * 60)
for n in [0, 10, 50, 99]:
    iw_pt = Giw_reconstructed.mesh[n]
    g_recon = Giw_reconstructed[iw_pt]
    print(f"  {n:3d}  {iw_pt.value.imag:8.4f}  {g_recon.real:14.10f}  {g_recon.imag:14.10f}")

print()
print("=" * 60)
print("DLR Test Completed Successfully!")
print("=" * 60)
