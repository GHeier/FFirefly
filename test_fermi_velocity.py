#!/usr/bin/env python3
"""Test Hamiltonian.get_fermi_velocity() in Python"""

import sys
sys.path.insert(0, '/home/g/Research/FFirefly/pypkg')

from firefly.src.module.imports.cpp_imports import Hamiltonian, load_config
import numpy as np

print("="*60)
print("Testing Hamiltonian.get_fermi_velocity in Python")
print("="*60)

# Load configuration
load_config("/home/g/Research/Materials/Tight_Binding/test/sample.cfg")

# Create Hamiltonian
print("\nInitializing Hamiltonian...")
H = Hamiltonian()
print(f"Hamiltonian loaded from file: {H.file_found}")

# Test 1: Single k-point
print("\n" + "="*60)
print("Test 1: Single k-point")
print("="*60)

k = [0.1, 0.2, 0.0]
print(f"k-point: {k}")

# Get bands
bands = H.get_bands(k)
print(f"Number of bands: {len(bands)}")
print(f"Bands: {bands}")

# Get Fermi velocity
vels = H.get_fermi_velocity(k)
print(f"\nFermi velocities (shape: {vels.shape}):")
for i in range(vels.shape[0]):
    v = vels[i, :]
    v_norm = np.linalg.norm(v)
    print(f"  Band {i+1}: v = [{v[0]:.6f}, {v[1]:.6f}, {v[2]:.6f}], |v| = {v_norm:.6f}")

# Test 2: Multiple k-points
print("\n" + "="*60)
print("Test 2: Multiple k-points")
print("="*60)

k_points = [
    [0.0, 0.0, 0.0],
    [0.5, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.0]
]
print(f"Number of k-points: {len(k_points)}")

# Get velocities for all k-points
vels_list = H.get_fermi_velocity(k_points)
print(f"Velocities array shape: {vels_list.shape}")

for i, k in enumerate(k_points):
    print(f"\nk[{i+1}] = {k}")
    for n in range(vels_list.shape[1]):
        v = vels_list[i, n, :]
        v_norm = np.linalg.norm(v)
        print(f"  Band {n+1}: v = [{v[0]:.4f}, {v[1]:.4f}, {v[2]:.4f}], |v| = {v_norm:.4f}")

# Test 3: Verify numerical derivative
print("\n" + "="*60)
print("Test 3: Verify numerical derivative")
print("="*60)

k0 = [0.3, 0.4, 0.0]
dk = 0.001

print(f"Testing at k = {k0}")
print(f"Using finite difference step dk = {dk}")

# Get Fermi velocity from our function
v_computed = H.get_fermi_velocity(k0)

# Compute numerical derivative manually
E0 = H.get_bands(k0)
nbands = len(E0)

kx_plus = [k0[0] + dk, k0[1], k0[2]]
ky_plus = [k0[0], k0[1] + dk, k0[2]]

Ex_plus = H.get_bands(kx_plus)
Ey_plus = H.get_bands(ky_plus)

v_manual = np.zeros((nbands, 3), dtype=np.float32)
for n in range(nbands):
    v_manual[n, 0] = (Ex_plus[n] - E0[n]) / dk
    v_manual[n, 1] = (Ey_plus[n] - E0[n]) / dk
    v_manual[n, 2] = 0.0

print("\nComparison:")
for n in range(nbands):
    v_comp = v_computed[n, :]
    v_man = v_manual[n, :]
    diff = np.linalg.norm(v_comp - v_man)
    print(f"  Band {n+1}:")
    print(f"    Computed: [{v_comp[0]:.6f}, {v_comp[1]:.6f}, {v_comp[2]:.6f}]")
    print(f"    Manual:   [{v_man[0]:.6f}, {v_man[1]:.6f}, {v_man[2]:.6f}]")
    print(f"    Difference: {diff:.8f}")

print("\n" + "="*60)
print("Python test completed successfully!")
print("="*60)
