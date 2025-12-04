#!/usr/bin/env python3
"""
Plot comparison of tetrahedron method vs numerical integration
for generalized susceptibility χ(q)
"""

import numpy as np
import matplotlib.pyplot as plt

# Read data
data = np.loadtxt('chi_comparison.dat')
qx = data[:, 0]  # qx/π
chi_tetra = data[:, 1]
chi_num = data[:, 2]
ratio = data[:, 3]

# Create figure with 3 subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10))

# Plot 1: Both methods
ax1.plot(qx, chi_tetra, 'o-', label='Tetrahedron Method', markersize=4)
ax1.plot(qx, chi_num, 's-', label='Numerical Integration', markersize=4)
ax1.set_xlabel('$q_x/\pi$')
ax1.set_ylabel('$\chi(q)$')
ax1.set_title('Generalized Susceptibility: 2D Tight-Binding Model')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Just numerical (easier to see structure)
ax2.plot(qx, chi_num, 's-', color='C1', markersize=4)
ax2.set_xlabel('$q_x/\pi$')
ax2.set_ylabel('$\chi(q)$ [Numerical]')
ax2.set_title('Numerical Integration Result')
ax2.grid(True, alpha=0.3)

# Plot 3: Ratio
ax3.plot(qx, ratio, 'd-', color='C2', markersize=4)
ax3.set_xlabel('$q_x/\pi$')
ax3.set_ylabel('Ratio (Tetrahedron/Numerical)')
ax3.set_title('Ratio of Methods')
ax3.grid(True, alpha=0.3)
ax3.axhline(y=1, color='k', linestyle='--', alpha=0.5, label='Perfect agreement')
ax3.legend()

plt.tight_layout()
plt.savefig('chi_comparison.png', dpi=150)
print("Plot saved to chi_comparison.png")

# Print some statistics
print(f"\nStatistics:")
print(f"  Mean ratio: {np.mean(ratio[ratio>0]):.2f}")
print(f"  Std ratio: {np.std(ratio[ratio>0]):.2f}")
print(f"  Min ratio: {np.min(ratio[ratio>0]):.2f}")
print(f"  Max ratio: {np.max(ratio[ratio>0]):.2f}")
