#!/usr/bin/env python3
# Test k-point ordering

import numpy as np

nk = 5

# Python approach
kx = np.linspace(-0.5, 0.5, nk)
ky = np.linspace(-0.5, 0.5, nk)
kgrid = np.meshgrid(kx, ky, indexing='ij')
kgrid_flat = np.stack(kgrid, axis=-1).reshape(-1, 2)

print("Python k-points (first 10):")
for i in range(10):
    print(f"  {i+1}: {kgrid_flat[i]}")
