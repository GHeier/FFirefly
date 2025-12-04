#!/usr/bin/env python3
"""
Compare Lindhard functions for different tight-binding parameters side-by-side.
Shows colorplots only (no high-symmetry path plots).
"""

import numpy as np
import matplotlib.pyplot as plt
from fermi_surface_integral import compute_lindhard_map


def plot_comparison():
    """
    Create side-by-side colorplots comparing two parameter sets.
    """
    print("\nGenerating side-by-side comparison plots...")

    # Common parameters
    t = 1.0
    mu = 0.0
    T = 0.1

    # Two different parameter sets
    params = [
        {'tp': 0.0, 'tpp': 0.0, 'label': "t'=0, t''=0"},
        {'tp': -0.2, 'tpp': 0.16, 'label': "t'=-0.2, t''=0.16"}
    ]

    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    for idx, param_set in enumerate(params):
        tp = param_set['tp']
        tpp = param_set['tpp']
        label = param_set['label']
        ax = axes[idx]

        print(f"\nComputing for {label}...")
        print(f"  Parameters: t={t}, t'={tp}, t''={tpp}, μ={mu}, T={T}")

        # Compute Lindhard function
        qx, qy, chi_q = compute_lindhard_map(
            q_max=np.pi, nq=50, t=t, tp=tp, mu=mu, T=T, nk=200
        )

        # Create meshgrid for plotting
        QX, QY = np.meshgrid(qx, qy, indexing='ij')

        # Create color plot
        im = ax.pcolormesh(QX, QY, chi_q, cmap='viridis', shading='auto')

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(r'$|\chi(q)|$', fontsize=12, rotation=270, labelpad=20)

        # Labels and title
        ax.set_xlabel(r'$q_x$', fontsize=14)
        ax.set_ylabel(r'$q_y$', fontsize=14)
        ax.set_title(f'{label}\nT = {T}', fontsize=14, fontweight='bold')
        ax.set_aspect('equal')

        # Add grid and reference lines
        ax.axhline(0, color='white', linewidth=1, alpha=0.5, linestyle='--')
        ax.axvline(0, color='white', linewidth=1, alpha=0.5, linestyle='--')

        # Mark high-symmetry points
        ax.plot(0, 0, 'r*', markersize=12, label=r'$\Gamma$ (0,0)')
        ax.plot(np.pi, 0, 'g*', markersize=12, label=r'X ($\pi$,0)')
        ax.plot(np.pi, np.pi, 'b*', markersize=12, label=r'M ($\pi$,$\pi$)')
        ax.legend(loc='upper right', fontsize=9)

        # Set ticks at multiples of π
        pi_ticks = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
        pi_labels = [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
        ax.set_xticks(pi_ticks)
        ax.set_xticklabels(pi_labels)
        ax.set_yticks(pi_ticks)
        ax.set_yticklabels(pi_labels)

        # Print statistics
        print(f"  Maximum value: {chi_q.max():.6f} at q = ({qx[np.unravel_index(chi_q.argmax(), chi_q.shape)[0]]:.3f}, {qy[np.unravel_index(chi_q.argmax(), chi_q.shape)[1]]:.3f})")
        print(f"  Minimum value: {chi_q.min():.6f}")

    plt.tight_layout()
    output_file = 'lindhard_comparison.png'
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    plt.close()


if __name__ == "__main__":
    # Set matplotlib backend
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend

    # Generate comparison plot
    try:
        plot_comparison()
    except Exception as e:
        print(f"Plotting failed: {e}")
        raise
