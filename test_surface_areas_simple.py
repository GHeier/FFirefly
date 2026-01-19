#!/usr/bin/env python3
"""
Simple test of Surface.get_faces_and_areas() in Python
Uses a pre-computed Surface to test the new method
"""

import sys
sys.path.insert(0, '/home/g/Research/FFirefly/pypkg')

from firefly.src.module.imports.cpp_imports import Surface, load_config, Vec
import numpy as np

print("="*60)
print("Testing Surface.get_faces_and_areas() in Python")
print("="*60)

# Load configuration
load_config("/home/g/Research/Materials/Tight_Binding/test/sample.cfg")

# Create a simple tight-binding epsilon function
def eps_func(k):
    """2D tight-binding model: epsilon = -2t(cos(kx) + cos(ky))"""
    t = 1.0
    return -2.0 * t * (np.cos(k.x) + np.cos(k.y))

# Create Surface at Fermi energy
mu = -1.5
print(f"\nCreating surface at μ = {mu}")
surf = Surface(eps_func, mu)

# Test old method (faces attribute loaded in __init__)
print("\nTesting faces attribute (old method):")
faces = surf.faces
n_faces = len(faces)
print(f"  Number of faces: {n_faces}")

if n_faces > 0:
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
    print(f"  Total area: {sum(areas):.6f}")
    print(f"  Min area: {min(areas):.6f}")
    print(f"  Max area: {max(areas):.6f}")
    print(f"  Mean area: {np.mean(areas):.6f}")

    # Check that areas are positive
    all_positive = all(a > 0 for a in areas)
    print(f"  All areas positive: {all_positive}")

    print("\n" + "="*60)
    print("Python test completed successfully!")
    print("="*60)
else:
    print("  ERROR: No faces found!")
