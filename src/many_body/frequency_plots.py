"""
Frequency plotting utilities for many-body calculations.
Handles plotting of Green's functions, self-energies, and other frequency-dependent quantities.
"""

import numpy as np
import matplotlib.pyplot as plt
from triqs.gf import Gf
from triqs.gf.meshes import MeshDLRImFreq


def plot_frequency_data(G_iw, title="", plot_real=True, plot_imag=True,
                       positive_only=True, indices=(0, 0), figsize=(10, 6),
                       save_path=None, show=True):
    """
    Plot frequency-dependent data (Green's functions, self-energy, etc.).

    Parameters:
    -----------
    G_iw : Gf object or Diagram object
        TRIQS Green's function or Diagram with frequency mesh
    title : str
        Plot title
    plot_real : bool
        Whether to plot real part
    plot_imag : bool
        Whether to plot imaginary part
    positive_only : bool
        If True, only plot positive frequencies
    indices : tuple
        Matrix indices to plot (for matrix-valued functions)
    figsize : tuple
        Figure size (width, height)
    save_path : str, optional
        If provided, save figure to this path
    show : bool
        Whether to display the plot
    """
    # Handle Diagram objects
    if hasattr(G_iw, 'obj_w'):
        gf_data = G_iw.obj_w
    else:
        gf_data = G_iw

    # Extract frequency points
    mesh = gf_data.mesh

    # Get Matsubara frequencies
    w_points = np.array([float(iw.imag) for iw in mesh])

    # Extract data
    if len(gf_data.data.shape) == 1:
        # Scalar Green's function
        data = gf_data.data
    else:
        # Matrix Green's function - extract specific indices
        data = gf_data.data[:, indices[0], indices[1]]

    # Filter to positive frequencies if requested
    if positive_only:
        mask = w_points >= 0
        w_points = w_points[mask]
        data = data[mask]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Plot real part
    if plot_real:
        ax.plot(w_points, data.real, 'o-', label='Re', markersize=3, linewidth=1.5)

    # Plot imaginary part
    if plot_imag:
        ax.plot(w_points, data.imag, 's-', label='Im', markersize=3, linewidth=1.5)

    # Formatting
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.8)
    ax.set_xlabel(r'$\omega_n$', fontsize=12)
    ax.set_ylabel('Value', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    if title:
        ax.set_title(title, fontsize=13)

    plt.tight_layout()

    # Save if requested
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    # Show if requested
    if show:
        plt.show()

    return fig, ax


def plot_frequency_comparison(G_list, labels, title="Comparison",
                              plot_real=True, plot_imag=True,
                              positive_only=True, indices=(0, 0),
                              figsize=(12, 5), save_path=None, show=True):
    """
    Plot multiple frequency-dependent quantities for comparison.

    Parameters:
    -----------
    G_list : list
        List of Gf objects or Diagram objects to compare
    labels : list of str
        Labels for each quantity
    title : str
        Plot title
    plot_real : bool
        Whether to plot real parts
    plot_imag : bool
        Whether to plot imaginary parts
    positive_only : bool
        If True, only plot positive frequencies
    indices : tuple
        Matrix indices to plot (for matrix-valued functions)
    figsize : tuple
        Figure size (width, height)
    save_path : str, optional
        If provided, save figure to this path
    show : bool
        Whether to display the plot
    """
    n_plots = int(plot_real) + int(plot_imag)
    if n_plots == 0:
        raise ValueError("Must plot at least one of real or imaginary part")

    fig, axes = plt.subplots(1, n_plots, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    plot_idx = 0

    # Plot real parts
    if plot_real:
        ax = axes[plot_idx]
        for G_iw, label in zip(G_list, labels):
            # Handle Diagram objects
            if hasattr(G_iw, 'obj_w'):
                gf_data = G_iw.obj_w
            else:
                gf_data = G_iw

            # Extract frequency points
            mesh = gf_data.mesh
            w_points = np.array([float(iw.imag) for iw in mesh])

            # Extract data
            if len(gf_data.data.shape) == 1:
                data = gf_data.data
            else:
                data = gf_data.data[:, indices[0], indices[1]]

            # Filter to positive frequencies
            if positive_only:
                mask = w_points >= 0
                w_points = w_points[mask]
                data = data[mask]

            ax.plot(w_points, data.real, 'o-', label=label, markersize=3, linewidth=1.5)

        ax.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.8)
        ax.set_xlabel(r'$\omega_n$', fontsize=12)
        ax.set_ylabel('Re', fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_title(f'{title} - Real Part', fontsize=11)
        plot_idx += 1

    # Plot imaginary parts
    if plot_imag:
        ax = axes[plot_idx]
        for G_iw, label in zip(G_list, labels):
            # Handle Diagram objects
            if hasattr(G_iw, 'obj_w'):
                gf_data = G_iw.obj_w
            else:
                gf_data = G_iw

            # Extract frequency points
            mesh = gf_data.mesh
            w_points = np.array([float(iw.imag) for iw in mesh])

            # Extract data
            if len(gf_data.data.shape) == 1:
                data = gf_data.data
            else:
                data = gf_data.data[:, indices[0], indices[1]]

            # Filter to positive frequencies
            if positive_only:
                mask = w_points >= 0
                w_points = w_points[mask]
                data = data[mask]

            ax.plot(w_points, data.imag, 's-', label=label, markersize=3, linewidth=1.5)

        ax.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.8)
        ax.set_xlabel(r'$\omega_n$', fontsize=12)
        ax.set_ylabel('Im', fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_title(f'{title} - Imaginary Part', fontsize=11)

    plt.tight_layout()

    # Save if requested
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    # Show if requested
    if show:
        plt.show()

    return fig, axes


def plot_dmft_results(S, save_dir=None):
    """
    Plot standard DMFT results: G_loc, Sigma_loc, G_Weiss.

    Parameters:
    -----------
    S : ManyBodySolver
        Solver object with DMFT results
    save_dir : str, optional
        Directory to save plots
    """
    import os

    # Plot G_loc
    title = "Local Green's Function $G_{loc}(i\\omega_n)$"
    save_path = os.path.join(save_dir, "G_loc.png") if save_dir else None
    plot_frequency_data(S.G_loc, title=title, positive_only=True,
                       save_path=save_path, show=False)

    # Plot Sigma_loc
    title = "Local Self-Energy $\\Sigma_{loc}(i\\omega_n)$"
    save_path = os.path.join(save_dir, "Sigma_loc.png") if save_dir else None
    plot_frequency_data(S.Sigma_loc, title=title, positive_only=True,
                       save_path=save_path, show=False)

    # Plot G_Weiss if it exists
    if hasattr(S, 'G_Weiss'):
        title = "Weiss Field $G_0(i\\omega_n)$"
        save_path = os.path.join(save_dir, "G_Weiss.png") if save_dir else None
        plot_frequency_data(S.G_Weiss, title=title, positive_only=True,
                           save_path=save_path, show=False)

    # Comparison plot
    G_list = [S.G_loc, S.Sigma_loc]
    labels = ['$G_{loc}$', '$\\Sigma_{loc}$']
    if hasattr(S, 'G_Weiss'):
        G_list.append(S.G_Weiss)
        labels.append('$G_0$')

    save_path = os.path.join(save_dir, "DMFT_comparison.png") if save_dir else None
    plot_frequency_comparison(G_list, labels, title="DMFT Results",
                             positive_only=True, save_path=save_path)
