#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Read data
data = np.loadtxt('response_comparison.dat')
i = data[:, 0]
qx = data[:, 1]
ref = data[:, 2]
tet = data[:, 3]
ratio = data[:, 4]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Plot 1: Comparison
ax1.plot(qx, ref, 'o-', label='Reference (Nk=20000)', markersize=5)
ax1.plot(qx, tet, 's-', label='Tetrahedron (Nk=200)', markersize=4)
ax1.set_xlabel('$q_x$')
ax1.set_ylabel(r'Response $\chi(q)$')
ax1.set_title('Response Function Comparison')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_yscale('log')

# Plot 2: Ratio
ax2.plot(qx[1:], ratio[1:], 'd-', color='C2', markersize=4)
ax2.set_xlabel('$q_x$')
ax2.set_ylabel('Ratio (Tetrahedron/Reference)')
ax2.set_title('Method Ratio')
ax2.grid(True, alpha=0.3)
ax2.axhline(y=1, color='k', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('response_comparison.png', dpi=150)
print("Plot saved to response_comparison.png")

# Statistics
print(f"\nRatio statistics (excluding q=0):")
print(f"  Mean: {np.mean(ratio[1:]):.2f}")
print(f"  Std: {np.std(ratio[1:]):.2f}")
print(f"  Min: {np.min(ratio[1:]):.2f}")
print(f"  Max: {np.max(ratio[1:]):.2f}")
