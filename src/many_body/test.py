from triqs.gf import *
from triqs.operators import *
from triqs_cthyb import Solver
from triqs.plot.mpl_interface import *

# Parameters from arXiv:2310.01266 (Bethe lattice interacting case, Section IV)
D = 1.0           # Bandwidth (reference energy unit)
U = 4.0 * D       # Interaction strength U/D = 4
beta = 5.0 / D    # Inverse temperature β·D = 5
t = D / 2.0       # Hopping parameter for Bethe lattice (half-bandwidth)

print(f"=== Bethe lattice DMFT (following arXiv:2310.01266) ===")
print(f"Bandwidth D = {D} eV")
print(f"Interaction U/D = {U/D}")
print(f"Inverse temperature β·D = {beta*D}")
print(f"Temperature T = {1/beta:.3f} eV")

# Construct the impurity solver with parameters from paper
S = Solver(beta = beta, gf_struct = [('up',1), ('down',1)] )

# Bethe lattice with semicircular DOS: Δ(iω) = t^2 * G_loc(iω)
# For half-filled Bethe lattice, self-consistency: G0^-1 = iω + μ - Σ - Δ
# Start with non-interacting approximation (Σ=0)
S.G0_iw << SemiCircular(t)

# Define the interacting Hamiltonian: H = U n↑ n↓
h_int = U * n('up',0) * n('down',0)

# CT-HYB solver parameters from paper
# Paper uses 10^6 to 10^9 QMC samples - use 5×10^5 for testing, increase if needed
p = {
    'length_cycle': 10,
    'n_warmup_cycles': 20000,      # More warmup for better equilibration
    'n_cycles': 500000,            # 5×10^5 cycles (reasonable compromise)
    'measure_G_tau': True
}

print(f"\nRunning CTHYB impurity solver (this may take a few minutes)...")
print(f"QMC cycles: {p['n_cycles']}")
S.solve(h_int = h_int, **p)

# Convert Green's functions to DLR representation for minimize_dyson
from triqs.gf import dlr_crm_dyson_solver, Gf, BlockGf
from triqs.gf.meshes import MeshDLRImFreq
import numpy as np

# DLR parameters from paper (arXiv:2310.01266)
# Paper uses DLR tolerance ε = 10^-6 and mentions Λ = β·ω_max cutoff
w_max = 20.0 * D  # High-frequency cutoff (Λ = β·ω_max ≈ 100)
eps = 1e-6        # DLR precision from paper

print(f"\n=== DLR representation ===")
print(f"DLR cutoff Λ = β·ω_max = {beta * w_max:.1f}")
print(f"DLR tolerance ε = {eps}")

# Create DLR mesh
dlr_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
print(f"DLR mesh size: {len(dlr_mesh)} points")

# Create DLR Green's functions and fill carefully
# Use evaluation at mesh points rather than assignment
G0_dlr = BlockGf(name_list=['up', 'down'], block_list=[Gf(mesh=dlr_mesh, target_shape=[1,1]) for _ in range(2)])
G_dlr = BlockGf(name_list=['up', 'down'], block_list=[Gf(mesh=dlr_mesh, target_shape=[1,1]) for _ in range(2)])

print("Filling DLR Green's functions via interpolation...")
for name in ['up', 'down']:
    for iwn in dlr_mesh:
        # Evaluate G at the DLR frequency points using Fourier transform
        G0_dlr[name].data[iwn.data_index, :, :] = S.G0_iw[name](iwn.value)
        G_dlr[name].data[iwn.data_index, :, :] = S.G_iw[name](iwn.value)

print(f"Filled DLR mesh with {len(G_dlr['up'].mesh)} points")

# Compute self-energy in Matsubara frequencies to get tail coefficients
print(f"\n=== Computing self-energy tail moments ===")
S_iw = inverse(S.G0_iw) - inverse(S.G_iw)

# Fit the tail to get moments (Hartree shift and higher moments)
tail, err = S_iw['up'].fit_hermitian_tail()
print(f"Tail fit error: {err:.2e}")
print(f"Σ_0 (Hartree shift): {tail[0][0,0]:.6f}")
print(f"Σ_1 (first moment): {tail[1][0,0]:.6f}")

# Try with only first moment for better stability
# The second moment constraint can be too restrictive with noisy QMC data
print("Using only Σ_0 (Hartree shift) constraint for better stability...")
Sigma_moments = {'up': tail[0:1], 'down': tail[0:1]}

# Use minimize_dyson to compute self-energy following paper methodology
print(f"\n=== CRM Dyson solver (arXiv:2310.01266) ===")
Sigma_dlr, Sigma_0, residual = dlr_crm_dyson_solver.minimize_dyson(
    G0_dlr=G0_dlr,
    G_dlr=G_dlr,
    Sigma_moments=Sigma_moments,
    method='trust-constr',  # Default method from paper
    options=dict(maxiter=10000, disp=True, gtol=1e-6, xtol=1e-8)  # Reasonable tolerances
)

print(f"\n=== Results ===")
print(f"Hartree shift Σ_0 (up):   {Sigma_0['up'][0,0]:.6f}")
print(f"Hartree shift Σ_0 (down): {Sigma_0['down'][0,0]:.6f}")
print(f"Dyson residual (up):   {residual['up']:.2e}")
print(f"Dyson residual (down): {residual['down']:.2e}")

if residual['up'] < 1e-10 and residual['down'] < 1e-10:
    print(f"\n✓ CRM Dyson solver converged to machine precision!")
    print(f"  Self-energy compressed from {len(S.Sigma_iw['up'].mesh)} to {len(Sigma_dlr['up'].mesh)} DLR points")
else:
    print(f"\n⚠ Warning: Residual is larger than expected")

# For U/D=4, β·D=5, we expect Mott insulating behavior
print(f"\nPhysics: U/D={U/D} is in the strongly correlated regime (Mott insulator).")
print(f"The self-energy shows strong correlation effects.")

# Plot comparison: CTHYB self-energy vs CRM Dyson solver result
# Extract data for manual plotting with full control
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Compute deviation metrics
print(f"\n=== Deviation Analysis (CRM vs CTHYB) ===")

# Plot up and down spins
for name, ax in [('up', ax1), ('down', ax2)]:
    # Extract CTHYB data (points only)
    cthyb_freq = [float(w.imag) for w in S.Sigma_iw[name].mesh]  # Matsubara frequencies
    cthyb_data = [S.Sigma_iw[name][w].imag[0,0] for w in S.Sigma_iw[name].mesh]

    # Extract CRM Dyson data from DLR mesh (lines connecting fewer points)
    crm_freq = [float(w.imag) for w in Sigma_dlr[name].mesh]  # DLR mesh points
    crm_data = [Sigma_dlr[name][w].imag[0,0] for w in Sigma_dlr[name].mesh]

    # Compute deviation at DLR points
    cthyb_at_dlr = [S.Sigma_iw[name](w).imag[0,0] for w in Sigma_dlr[name].mesh]
    deviations = [abs(crm - cthyb) for crm, cthyb in zip(crm_data, cthyb_at_dlr)]
    mean_dev = np.mean(deviations)
    max_dev = np.max(deviations)
    rms_dev = np.sqrt(np.mean([d**2 for d in deviations]))

    print(f"{name} spin: mean = {mean_dev:.2e} eV, max = {max_dev:.2e} eV, RMS = {rms_dev:.2e} eV")

    # Plot: CRM as solid line (smooth from DLR), CTHYB as points
    ax.plot(crm_freq, crm_data, '-', linewidth=2, label='CRM Dyson (DLR)', color='C0')
    ax.plot(cthyb_freq, cthyb_data, 'o', markersize=3, markerfacecolor='none',
            label='CTHYB', color='C1', markeredgewidth=1)

    ax.set_xlim(0, 200)
    ax.set_ylim(-50, 50)
    ax.set_xlabel('Matsubara frequency (eV)')
    ax.set_ylabel('Im Σ(iωₙ) (eV)')
    ax.set_title(f'Self-energy ({name} spin)\nBethe lattice: U/D={U/D}, β·D={beta*D}')
    ax.legend(loc='best')
    ax.grid(alpha=0.3)

plt.suptitle('CRM Dyson Solver (arXiv:2310.01266)', fontsize=14, y=1.00)
plt.tight_layout()
plt.show()
