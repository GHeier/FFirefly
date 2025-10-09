#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from triqs.gf import *
from triqs.gf.meshes import MeshDLRImFreq, MeshReFreq, MeshImFreq


def iw_to_tau_dlr(Giw_dlr):
    G_dlr_coeff = make_gf_dlr(Giw_dlr)
    Gtau_dlr = make_gf_dlr_imtime(G_dlr_coeff)
    return Gtau_dlr

def tau_to_iw_dlr(Gtau_dlr):
    G_dlr = make_gf_dlr(Gtau_dlr)
    G_iw = make_gf_dlr_imfreq(G_dlr)
    return G_iw

def IPT_iter(G0w, Gw, Sigma_w, U):
    # Transform G0(iω) to G0(τ)
    Gtau = iw_to_tau_dlr(G0w)

    # Calculate self-energy Σ(τ) = U²G0(τ)³
    Sigma_tau = Gtau.copy()
    Sigma_tau << (U**2) * Gtau * Gtau * Gtau

    # Transform Σ(τ) → Σ(iω)
    Sigma_w = tau_to_iw_dlr(Sigma_tau)

    # Calculate new G(iω) using Dyson equation
    Gw << inverse(inverse(G0w) - Sigma_w)

def calculate_dos(U, beta, w_max, eps, n_loops):
    if U == 0:
        n_loops = 0
    # Create DLR mesh
    dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
    # Create initial G(iω) with semicircular DOS
    Giw_dlr = Gf(mesh=dlr_iw_mesh, target_shape=[])
    Giw_dlr << SemiCircular(D)
    Sigma_iw = Giw_dlr.copy()
    G_new = Giw_dlr.copy()
    #G_new << inverse( iOmega_n - (D/2)**2 * Giw_dlr )

    for _ in range(n_loops):
        IPT_iter(Giw_dlr, G_new, Sigma_iw,  U)

    # Convert DLR to standard ImFreq mesh for Pade
    n_iw_standard = 100
    imfreq_mesh = MeshImFreq(beta=beta, statistic='Fermion', n_iw=n_iw_standard)
    Giw_standard = Gf(mesh=imfreq_mesh, target_shape=[])

    # Sample DLR Green's function on standard mesh
    G_dlr_coeff = make_gf_dlr(G_new)
    Giw_temp = make_gf_imfreq(G_dlr_coeff, n_iw=n_iw_standard)

    # Create real frequency mesh
    w_min, w_max_freq, n_w = -6.0, 6.0, 400
    refreq_mesh = MeshReFreq(w_min, w_max_freq, n_w)
    Gw = Gf(mesh=MeshReFreq(window = (-5.0,5.0), n_w=1000), target_shape=[1,1])

    # Perform Pade continuation
    n_pade = min(60, n_iw_standard)
    Gw.set_from_pade(Giw_temp, n_points=n_pade, freq_offset=0.0)

    # Calculate DOS: A(ω) = -Im[G(ω)]/π
    omega = np.array([w.value for w in Gw.mesh])
    dos = -Gw.data.imag / np.pi

    return omega, dos.real


# Main execution
if __name__ == "__main__":
    print("=" * 80)
    print("IPT with DLR: Comparing U=0 and U=1.5")
    print("=" * 80)

    # Parameters
    beta = 20.0
    D = 2
    w_max = 1.2*D
    eps = 1e-14
    n_loops = 25

    print(f"\nParameters: β={beta}, w_max={w_max}, ε={eps}\n")

    # Calculate DOS for U=0 (non-interacting)
    print("Calculating DOS for U=0 (non-interacting)...")
    dos_list = list()
    U_list = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    for U in U_list:
        omega, dos_U0 = calculate_dos(U=U, beta=beta, w_max=w_max, eps=eps, n_loops=n_loops)
        dos_list.append(dos_U0)

    print("\nPlotting results...\n")

    # Plot results
    fig, (ax1) = plt.subplots(1, 1, figsize=(12, 4))

    # Plot DOS - linear scale
    for i in range(len(U_list)):
        ax1.plot(omega, dos_list[i], '-', color='blue', linewidth=2, label=f'U={U_list[i]} ', alpha=0.5)
    ax1.axvline(0, color='blue', linestyle=':', alpha=0.5)
    ax1.set_xlabel('ω', fontsize=12)
    ax1.set_ylabel('DOS(ω)', fontsize=12)
    ax1.set_title('Density of States', fontsize=13, fontweight='bold')
    ax1.grid(alpha=0.3)
    #ax1.legend(fontsize=11)
    ax1.set_xlim(-8, 8)
    ax1.set_ylim(bottom=0)

    ## Plot spectral function - log scale
    #for i in range
    #ax2.semilogy(omega, dos_U0, 'k--', linewidth=2, label='U=0', alpha=0.7)
    #ax2.semilogy(omega, dos_U15, 'r-', linewidth=2, label='U=1.5')
    #ax2.axvline(0, color='gray', linestyle=':', alpha=0.5)
    #ax2.set_xlabel('ω', fontsize=12)
    #ax2.set_ylabel('DOS(ω) [log scale]', fontsize=12)
    #ax2.set_title('Spectral Function (log)', fontsize=13, fontweight='bold')
    #ax2.grid(alpha=0.3, which='both')
    #ax2.legend(fontsize=11)
    #ax2.set_xlim(-3, 3)
    #ax2.set_ylim(1e-2, 200)

    plt.tight_layout()
    plt.savefig('ipt_dos_comparison.png', dpi=150, bbox_inches='tight')
    print(f"✓ Saved DOS comparison plot to ipt_dos_comparison.png")

    # Print results
    print(f"\n{'='*80}")
    print("Results Summary")
    print(f"{'='*80}")
    for i in range(len(U_list)):
        print(f"U={U_list[i]}: Peak DOS = {dos_list[i].max():.4f} at ω = {omega[dos_list[i].argmax()]:.3f}")
    print(f"{'='*80}\n")
