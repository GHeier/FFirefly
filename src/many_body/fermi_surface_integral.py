#!/usr/bin/env python3
"""
Compute integral of f(ε(k+q)) - f(ε(k)) over k for 2D tight-binding model.

This is useful for calculating susceptibility and response functions in many-body theory.
The integral appears in Lindhard functions and related quantities.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# ==================== CONSTANTS ====================
T = 0.01      # Temperature in energy units
mu1 = -0.9      # Chemical potential (Fermi energy)
mu2 = -1.3      # Chemical potential (Fermi energy)
nk = 150      # Number of k-points per dimension
# ==================================================


def fermi_dirac(energy, mu=0.0, T=0.01):
    """
    Fermi-Dirac distribution function.

    Parameters:
    -----------
    energy : float or array
        Energy values
    mu : float
        Chemical potential (Fermi energy)
    T : float
        Temperature in energy units

    Returns:
    --------
    f : float or array
        Occupation number f(E) = 1/(1 + exp((E-mu)/T))
    """
    if T < 1e-10:
        # Zero temperature limit
        return np.where(energy < mu, 1.0, 0.0)

    # Prevent overflow in exponential
    x = (energy - mu) / T
    x = np.clip(x, -50, 50)
    return 1.0 / (1.0 + np.exp(x))


def tight_binding_2d(kx, ky, t=1.0, tp=0.0, tpp = 0.0, mu=0.0):
    """
    2D tight-binding dispersion relation.

    Parameters:
    -----------
    kx, ky : float or array
        Momentum components
    t : float
        Nearest-neighbor hopping
    tp : float
        Next-nearest-neighbor hopping
    tpp : float
        Third-nearest-neighbor hopping
    mu : float
        Chemical potential

    Returns:
    --------
    energy : float or array
        Band energy ε(k) = -2t(cos(kx) + cos(ky)) - 4tp*cos(kx)*cos(ky) - 2tpp*(cos(2kx) + cos(2ky)) - mu
    """
    energy = -2.0 * t * (np.cos(kx) + np.cos(ky))
    if tp != 0.0:
        energy -= 4.0 * tp * np.cos(kx) * np.cos(ky)
    if tpp != 0.0:
        energy -= 2.0 * tpp * (np.cos(2*kx) + np.cos(2*ky))
    return energy - mu


def integrand_2d(kx, ky, qx, qy, t=1.0, tp=0.0, tpp=0.0, mu=0.0, T=0.01):
    """
    Compute f(ε(k+q)) - f(ε(k)) for given k and q.

    Parameters:
    -----------
    kx, ky : float or array
        Momentum components
    qx, qy : float
        Momentum transfer
    t, tp, tpp : float
        Hopping parameters
    mu : float
        Chemical potential
    T : float
        Temperature

    Returns:
    --------
    difference : float or array
        f(ε(k+q)) - f(ε(k))
    """
    # Energy at k
    ek = tight_binding_2d(kx, ky, t, tp, tpp, mu)

    # Energy at k+q (with periodic boundary conditions)
    ekq = tight_binding_2d(kx + qx, ky + qy, t, tp, tpp, mu)

    # Fermi function difference
    fk = fermi_dirac(ek, mu=0.0, T=T)  # mu already subtracted in dispersion
    fkq = fermi_dirac(ekq, mu=0.0, T=T)

    return abs(fkq - fk)


def compute_integral_grid(qx, qy, t=1.0, tp=0.0, tpp=0.0, mu=0.0, T=0.01, nk=200):
    """
    Compute integral using uniform grid sampling.

    Integral = (1/(2π)²) ∫∫ dk_x dk_y [f(ε(k+q)) - f(ε(k))]

    Parameters:
    -----------
    qx, qy : float
        Momentum transfer
    t, tp, tpp : float
        Hopping parameters
    mu : float
        Chemical potential
    T : float
        Temperature
    nk : int
        Number of k-points per dimension

    Returns:
    --------
    integral : float
        Value of the integral
    """
    # Create k-space grid over first Brillouin zone [-π, π]
    kx = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    ky = np.linspace(-np.pi, np.pi, nk, endpoint=False)

    # 2D grid
    KX, KY = np.meshgrid(kx, ky, indexing='ij')

    # Compute integrand
    integrand_vals = integrand_2d(KX, KY, qx, qy, t, tp, tpp, mu, T)

    # Integrate using trapezoid rule
    # Volume element: dk_x dk_y / (2π)²
    dk = kx[1] - kx[0]
    integral = np.sum(integrand_vals) * dk**2 / (2 * np.pi)**2

    return integral


def compute_lindhard_map(q_max=np.pi, nq=50, t=1.0, tp=0.0, tpp=0.0, mu=0.0, T=0.01, nk=200):
    """
    Compute the integral as a function of momentum transfer q.
    This gives the Lindhard function (bare susceptibility).

    Parameters:
    -----------
    q_max : float
        Maximum q value
    nq : int
        Number of q points per dimension
    t, tp, tpp : float
        Hopping parameters
    mu : float
        Chemical potential
    T : float
        Temperature
    nk : int
        k-space grid size

    Returns:
    --------
    qx_array, qy_array : ndarray
        q-space grid
    chi_q : ndarray
        Lindhard function values
    """
    qx = np.linspace(-q_max, q_max, nq)
    qy = np.linspace(-q_max, q_max, nq)
    QX, QY = np.meshgrid(qx, qy, indexing='ij')

    chi_q = np.zeros_like(QX)

    print(f"Computing Lindhard function on {nq}×{nq} q-grid...")
    for i in range(nq):
        if i % 10 == 0:
            print(f"  Progress: {i}/{nq}")
        for j in range(nq):
            chi_q[i, j] = compute_integral_grid(QX[i, j], QY[i, j], t, tp, tpp, mu, T, nk)

    return qx, qy, chi_q


def plot_lindhard_function():
    """
    Visualize the Lindhard function in q-space.
    """
    print("\nGenerating Lindhard function plot...")

    # Parameters
    t = 1.0
    tp = 0.0
    tpp = 0.0

    # Compute on q-grid
    qx, qy, chi_q = compute_lindhard_map(
        q_max=np.pi, nq=30, t=t, tp=tp, tpp=tpp, mu=mu, T=T, nk=nk
    )

    # Create plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # 2D contour plot
    QX, QY = np.meshgrid(qx, qy, indexing='ij')
    levels = np.linspace(chi_q.min(), chi_q.max(), 20)
    cs = ax1.contourf(QX, QY, chi_q, levels=levels, cmap='RdBu_r')
    ax1.contour(QX, QY, chi_q, levels=levels, colors='k', linewidths=0.5, alpha=0.3)
    plt.colorbar(cs, ax=ax1, label=r'$\chi(q)$')
    ax1.set_xlabel(r'$q_x$')
    ax1.set_ylabel(r'$q_y$')
    ax1.set_title(f'Lindhard function (T={T}, μ={mu})')
    ax1.set_aspect('equal')
    ax1.axhline(0, color='k', linewidth=0.5, alpha=0.5)
    ax1.axvline(0, color='k', linewidth=0.5, alpha=0.5)

    # Cut along high-symmetry directions
    nq = len(qx)
    mid = nq // 2

    # (0,0) -> (π,0) -> (π,π) -> (0,0)
    q_path = []
    chi_path = []

    # Γ to X: (0,0) to (π,0)
    for i in range(mid, nq):
        q_path.append(qx[i] - qx[mid])
        chi_path.append(chi_q[i, mid])

    # X to M: (π,0) to (π,π)
    offset = qx[-1] - qx[mid]
    for j in range(mid, nq):
        q_path.append(offset + (qy[j] - qy[mid]))
        chi_path.append(chi_q[-1, j])

    # M to Γ: (π,π) to (0,0)
    offset = q_path[-1]
    for i in range(nq-1, mid-1, -1):
        dist = np.sqrt((qx[i] - qx[mid])**2 + (qy[i] - qy[mid])**2)
        q_path.append(offset + dist)
        chi_path.append(chi_q[i, i])

    ax2.plot(q_path, chi_path, 'b-', linewidth=2)
    ax2.set_xlabel('Path in q-space')
    ax2.set_ylabel(r'$\chi(q)$')
    ax2.set_title('High-symmetry path')
    ax2.grid(True, alpha=0.3)

    # Mark special points
    x_point = qx[-1] - qx[mid]
    m_point = offset + (qy[-1] - qy[mid])
    ax2.axvline(0, color='r', linestyle='--', alpha=0.5, label='Γ')
    ax2.axvline(x_point, color='g', linestyle='--', alpha=0.5, label='X')
    ax2.axvline(m_point, color='b', linestyle='--', alpha=0.5, label='M')
    ax2.legend()

    plt.tight_layout()
    plt.savefig('lindhard_function.png', dpi=150, bbox_inches='tight')
    print(f"Plot saved to: lindhard_function.png")
    plt.close()


def compute_filling(t=1.0, tp=0.0, tpp=0.0, mu=0.0, T=0.01, nk=200):
    """
    Compute the electron filling factor n = <n_k> via k-space integration.

    Parameters:
    -----------
    t, tp, tpp : float
        Hopping parameters
    mu : float
        Chemical potential
    T : float
        Temperature
    nk : int
        k-space grid size

    Returns:
    --------
    n : float
        Filling factor (average occupation per site)
    """
    # Create k-space grid
    kx = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    ky = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    KX, KY = np.meshgrid(kx, ky, indexing='ij')

    # Compute dispersion
    ek = tight_binding_2d(KX, KY, t, tp, tpp, mu)

    # Compute Fermi function
    fk = fermi_dirac(ek, mu=0.0, T=T)  # mu already subtracted in dispersion

    # Integrate over BZ
    dk = kx[1] - kx[0]
    n = np.sum(fk) * dk**2 / (2 * np.pi)**2

    return n


def plot_2d_colormap():
    """
    Create side-by-side comparison of Lindhard functions for two different band structures.
    """
    print("\nGenerating side-by-side 2D color plots...")

    # Band structure 1: nearest-neighbor only
    t1 = 1.0
    tp1 = 0.0
    tpp1 = 0.0

    # Band structure 2: with next-nearest-neighbor
    t2 = 1.0
    tp2 = -0.3
    tpp2 = 0.0

    # Compute filling factors
    print(f"\nBand structure 1: t={t1}, t'={tp1}, t''={tpp1}, μ={mu1}, T={T}, nk={nk}")
    n1 = compute_filling(t=t1, tp=tp1, tpp=tpp1, mu=mu1, T=T, nk=nk)
    print(f"  Filling n1 = {n1:.4f}")

    print(f"\nBand structure 2: t={t2}, t'={tp2}, t''={tpp2}, μ={mu2}, T={T}, nk={nk}")
    n2 = compute_filling(t=t2, tp=tp2, tpp=tpp2, mu=mu2, T=T, nk=nk)
    print(f"  Filling n2 = {n2:.4f}")

    # Compute Lindhard functions
    print(f"\nComputing Lindhard function for band structure 1...")
    qx1, qy1, chi_q1 = compute_lindhard_map(
        q_max=np.pi, nq=40, t=t1, tp=tp1, tpp=tpp1, mu=mu1, T=T, nk=nk
    )

    print(f"\nComputing Lindhard function for band structure 2...")
    qx2, qy2, chi_q2 = compute_lindhard_map(
        q_max=np.pi, nq=40, t=t2, tp=tp2, tpp=tpp2, mu=mu2, T=T, nk=nk
    )

    # Create meshgrids for plotting
    QX1, QY1 = np.meshgrid(qx1, qy1, indexing='ij')
    QX2, QY2 = np.meshgrid(qx2, qy2, indexing='ij')

    # Create figure with two panels
    fig = plt.figure(figsize=(18, 8))

    # Add title at the top with parameters
    fig.suptitle(f'T = {T},  μ1 = {mu1}, μ2 = {mu2}', fontsize=16, y=0.98, weight='bold')

    # Left panel: Band structure 1
    ax1 = plt.subplot(1, 2, 1)
    im1 = ax1.pcolormesh(QX1, QY1, chi_q1, cmap='viridis', shading='auto')
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_label(r'$|\chi(q)|$', fontsize=12, rotation=270, labelpad=20)

    ax1.set_xlabel(r'$q_x$', fontsize=14)
    ax1.set_ylabel(r'$q_y$', fontsize=14)
    ax1.set_title(f't = {t1},  t\' = {tp1},  t\'\' = {tpp1}\nn = {n1:.4f}',
                  fontsize=13, pad=10)
    ax1.set_aspect('equal')

    ax1.axhline(0, color='white', linewidth=1, alpha=0.5, linestyle='--')
    ax1.axvline(0, color='white', linewidth=1, alpha=0.5, linestyle='--')

    ax1.plot(0, 0, 'r*', markersize=12, label=r'$\Gamma$')
    ax1.plot(np.pi, 0, 'g*', markersize=12, label=r'X')
    ax1.plot(np.pi, np.pi, 'b*', markersize=12, label=r'M')
    ax1.legend(loc='upper right', fontsize=9)

    pi_ticks = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
    pi_labels = [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
    ax1.set_xticks(pi_ticks)
    ax1.set_xticklabels(pi_labels)
    ax1.set_yticks(pi_ticks)
    ax1.set_yticklabels(pi_labels)

    # Right panel: Band structure 2
    ax2 = plt.subplot(1, 2, 2)
    im2 = ax2.pcolormesh(QX2, QY2, chi_q2, cmap='viridis', shading='auto')
    cbar2 = plt.colorbar(im2, ax=ax2)
    cbar2.set_label(r'$|\chi(q)|$', fontsize=12, rotation=270, labelpad=20)

    ax2.set_xlabel(r'$q_x$', fontsize=14)
    ax2.set_ylabel(r'$q_y$', fontsize=14)
    ax2.set_title(f't = {t2},  t\' = {tp2},  t\'\' = {tpp2}\nn = {n2:.4f}',
                  fontsize=13, pad=10)
    ax2.set_aspect('equal')

    ax2.axhline(0, color='white', linewidth=1, alpha=0.5, linestyle='--')
    ax2.axvline(0, color='white', linewidth=1, alpha=0.5, linestyle='--')

    ax2.plot(0, 0, 'r*', markersize=12, label=r'$\Gamma$')
    ax2.plot(np.pi, 0, 'g*', markersize=12, label=r'X')
    ax2.plot(np.pi, np.pi, 'b*', markersize=12, label=r'M')
    ax2.legend(loc='upper right', fontsize=9)

    ax2.set_xticks(pi_ticks)
    ax2.set_xticklabels(pi_labels)
    ax2.set_yticks(pi_ticks)
    ax2.set_yticklabels(pi_labels)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('lindhard_function.png', dpi=200, bbox_inches='tight')
    print(f"\nPlot saved to: lindhard_function.png")
    print(f"\nBand structure 1 - Max χ: {chi_q1.max():.6f}, Min χ: {chi_q1.min():.6f}")
    print(f"Band structure 2 - Max χ: {chi_q2.max():.6f}, Min χ: {chi_q2.min():.6f}")
    plt.close()

    return {
        'chi_q1': chi_q1,
        'chi_q2': chi_q2,
        'qx': qx1,
        'qy': qy1,
        'params1': {'t': t1, 'tp': tp1, 'tpp': tpp1, 'mu': mu1, 'T': T, 'n': n1},
        'params2': {'t': t2, 'tp': tp2, 'tpp': tpp2, 'mu': mu2, 'T': T, 'n': n2}
    }


if __name__ == "__main__":
    # Set matplotlib backend
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    # Generate visualization
    print("\n" + "="*70)
    print("LINDHARD FUNCTION COMPARISON")
    print("="*70)
    try:
        results = plot_2d_colormap()
        print("\n" + "="*70)
        print("Computation complete!")
        print("="*70)
    except Exception as e:
        print(f"Plotting failed (non-critical): {e}")
        print("Skipping visualization.")
