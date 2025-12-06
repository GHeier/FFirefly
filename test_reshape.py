#!/usr/bin/env python3
import numpy as np

nk = 3

# Create flat array
vals = np.arange(9)
print("Flat values:", vals)

# Python reshape with C-order (row-major)
reshaped_c = vals.reshape(nk, nk, order='C')
print("\nPython reshape (C-order):")
print(reshaped_c)
print("Element [1,1]:", reshaped_c[1,1], "(should be 4)")
