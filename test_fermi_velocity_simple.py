#!/usr/bin/env python3
"""Simple test of Hamiltonian.get_fermi_velocity() in Python using ctypes directly"""

import ctypes
import numpy as np

# Load library
lib = ctypes.CDLL('/home/g/Research/FFirefly/build/lib/libfly.so')

# Load config
lib.load_config_export0.argtypes = [ctypes.c_char_p]
lib.load_config_export0.restype = None
lib.load_config_export0(b"/home/g/Research/Materials/Tight_Binding/test/sample.cfg")

print("="*60)
print("Testing Hamiltonian.get_fermi_velocity (Direct ctypes)")
print("="*60)

# Create Hamiltonian
lib.Hamiltonian_export0.argtypes = []
lib.Hamiltonian_export0.restype = ctypes.c_void_p
H_ptr = lib.Hamiltonian_export0()
print(f"\nHamiltonian created: {H_ptr != 0}")

# Test get_bands first
print("\n" + "="*60)
print("Test: get_bands")
print("="*60)

k = [0.1, 0.2, 0.0]
k_array = (ctypes.c_float * len(k))(*k)
k_len = ctypes.c_int(len(k))

max_bands = 100
eigenvalues_out = (ctypes.c_float * max_bands)()
num_bands = ctypes.c_int(0)

lib.Hamiltonian_get_bands_export0.argtypes = [
    ctypes.c_void_p,
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_int),
]
lib.Hamiltonian_get_bands_export0.restype = None

lib.Hamiltonian_get_bands_export0(
    H_ptr, k_array, k_len,
    eigenvalues_out, ctypes.byref(num_bands)
)

n = num_bands.value
print(f"Number of bands: {n}")
bands = [eigenvalues_out[i] for i in range(n)]
print(f"Bands at k={k}: {bands}")

# Test get_fermi_velocity
print("\n" + "="*60)
print("Test: get_fermi_velocity")
print("="*60)

velocities_out = (ctypes.c_float * (max_bands * 3))()
num_bands2 = ctypes.c_int(0)

lib.Hamiltonian_get_fermi_velocity_export0.argtypes = [
    ctypes.c_void_p,
    ctypes.POINTER(ctypes.c_float),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_float),
    ctypes.POINTER(ctypes.c_int),
]
lib.Hamiltonian_get_fermi_velocity_export0.restype = None

lib.Hamiltonian_get_fermi_velocity_export0(
    H_ptr, k_array, k_len,
    velocities_out, ctypes.byref(num_bands2)
)

n2 = num_bands2.value
print(f"Number of bands from velocity: {n2}")

vels = np.zeros((n2, 3), dtype=np.float32)
for i in range(n2):
    vels[i, 0] = velocities_out[i * 3 + 0]
    vels[i, 1] = velocities_out[i * 3 + 1]
    vels[i, 2] = velocities_out[i * 3 + 2]

print(f"\nFermi velocities:")
for i in range(n2):
    v_norm = np.linalg.norm(vels[i, :])
    print(f"  Band {i+1}: v = [{vels[i,0]:.6f}, {vels[i,1]:.6f}, {vels[i,2]:.6f}], |v| = {v_norm:.6f}")

print("\n" + "="*60)
print("Direct ctypes test completed successfully!")
print("="*60)
