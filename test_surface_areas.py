#!/usr/bin/env python3

import sys
sys.path.insert(0, '/home/g/Research/FFirefly/pypkg')

import firefly as ff
from firefly.src.module.imports.cpp_imports import epsilon
import numpy as np

print("="*60)
print("Testing Surface with get_faces_and_areas in Python")
print("="*60)

# Load configuration
from firefly.src.module.imports.cpp_imports import load_config
load_config("/home/g/Research/Materials/Tight_Binding/test/sample.cfg")

# Create epsilon function
def eps_func(k):
    """Epsilon function for tight binding model"""
    return epsilon(1, [k.x, k.y, k.z])

# Create Surface at Fermi energy
mu = -1.5
print(f"\nCreating surface at μ = {mu}")
from firefly.src.module.imports.cpp_imports import Surface
surf = Surface(eps_func, mu)

# Test old method (already loaded in __init__)
print("\nTesting faces attribute (old method):")
faces = surf.faces
print(f"  Number of faces: {len(faces)}")
print(f"  First face k-point: {faces[0]}")
print(f"  Dimension of first face: {len(faces[0])}")

# Test new method
print("\nTesting get_faces_and_areas (new method):")
kpoints, areas = surf.get_faces_and_areas()
print(f"  Number of k-points: {len(kpoints)}")
print(f"  Number of areas: {len(areas)}")
print(f"  First k-point: {kpoints[0]}")
print(f"  First area: {areas[0]}")

# Verify they match
print("\nVerification:")
kpoints_match = all(
    all(abs(k1 - k2) < 1e-6 for k1, k2 in zip(kp1, kp2))
    for kp1, kp2 in zip(kpoints, faces)
)
print(f"  K-points match: {kpoints_match}")
print(f"  Total area: {sum(areas)}")
print(f"  Min area: {min(areas)}")
print(f"  Max area: {max(areas)}")
print(f"  Mean area: {np.mean(areas)}")

print("\n" + "="*60)
print("Python test completed successfully!")
print("="*60)
