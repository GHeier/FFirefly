from triqs.gf import *
from triqs.operators import *
from triqs_cthyb import Solver
from h5 import HDFArchive
import triqs.utility.mpi as mpi
from triqs.gf.dlr_crm_dyson_solver import minimize_dyson
from triqs.gf import make_gf_from_fourier
from triqs.gf.gf_factories import fit_gf_dlr, make_gf_imfreq

# Parameters
D, V, U = 1.0, 0.2, 4.0
e_f, beta = -U/2.0, 50

# Construct the impurity solver with the inverse temperature
# and the structure of the Green's functions
S = Solver(beta = beta, gf_struct = [ ('up',1), ('down',1) ], n_l = 100)

# Initialize the non-interacting Green's function S.G0_iw
for name, g0 in S.G0_iw: g0 << inverse(iOmega_n - e_f - V**2 * SemiCircular(D/2))

G0_dlr_blocks = {}
G0_iw_recon_blocks = {}

for name, g0_iw in S.G0_iw:                    # g0_iw: GfImFreq (matrix-valued)
    # 1) ImFreq -> ImTime
    g0_tau = make_gf_from_fourier(g0_iw)

    # 2) DLR fit (choose a physical cutoff!)
    # If you used SemiCircular(D/2), the support is ~[-D/2, D/2] ⇒ pick w_max ≈ 0.5*D (or 1.1× for safety).
    w_max = 0.5 * D                            # spectral width cutoff
    eps   = 1e-10                              # tolerance

    g0_dlr = fit_gf_dlr(g0_tau, w_max, eps)    # matrix-valued gf<dlr>
    G0_dlr_blocks[name] = g0_dlr

    # 3) (optional) back to Matsubara, match original grid size
    n_iw = g0_iw.mesh.size
    g0_iw_recon = make_gf_imfreq(g0_dlr, n_iw=n_iw)
    G0_iw_recon_blocks[name] = g0_iw_recon

# Run the solver. The results will be in S.G_tau, S.G_iw and S.G_l
S.solve(h_int = U * n('up',0) * n('down',0),     # Local Hamiltonian
        n_cycles  = 500000,                      # Number of QMC cycles
        length_cycle = 200,                      # Length of one cycle
        n_warmup_cycles = 10000,                 # Warmup cycles
        measure_G_l = True)                      # Measure G_l

# ====================== FLAT BATH SOLVER PT 1 =================================

w_max = 20.0 * D  # High-frequency cutoff (Λ = β·ω_max ≈ 100)
eps = 1e-6        # DLR precision from paper
# Create DLR mesh
dlr_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
print(f"DLR mesh size: {len(dlr_mesh)} points")

from triqs.gf import MeshProduct, MeshBrillouinZone, Gf, fit_gf_dlr

# Compute self-energy in Matsubara frequencies to get tail coefficients
print(f"\n=== Computing self-energy tail moments ===")
S_iw = inverse(S.G0_iw) - inverse(S.G_iw)

# Fit the tail to get moments (Hartree shift and higher moments)
# For BlockGf, fit each block separately
tail, err = S_iw["up"].fit_hermitian_tail()
print(f"Tail fit error: {err:.2e}")

# Extract scalar moments (just the constant term) as arrays
print(S.G0_iw["up"])
g_up = S.G0_iw["up"].copy()   # <-- this is a Gf (ImFreq), independent copy
print(g_up)
G0_dlr = make_gf_dlr(S.G0_iw["up"])
G_dlr = make_gf_dlr(S.G_iw["up"])

Sigma_moments = {'up': tail[0:1], 'down': tail[0:1]}
S_iw_dlr, Sigma_HF, residual = minimize_dyson(G0_dlr=G0_dlr, G_dlr=G_dlr, Sigma_moments=tail[0:1])
# Compute self-energy in DLR representation
#S_iw_dlr = make_gf_dlr(S_iw)

print(f"\n=== Results ===")

import matplotlib.pyplot as plt
# Plot Matsubara self-energy for both spins
fig, ax = plt.subplots(1,2, figsize=(12,5))

# Plot up spin
ax1 = plt.subplot(1,2,1)
mesh_w = G0_dlr.mesh.components[0]
iw_arr = np.array([complex(iw) for iw in mesh_w], dtype=np.complex128)

# S_iw is matrix-valued
ax1.plot(iw_plot, S_iw['up'].data[iw_indices,0,0].imag, 'o-', label='QMC Σ(iωₙ)', markersize=4)
ax1.plot(iw_plot, S_iw_dlr['up'].data[iw_indices,0,0].imag, 'o-', label='QMC Σ(iωₙ)', linewidth=2)
ax1.set_xlabel('Matsubara frequency iωₙ')
ax1.set_ylabel('Im Σ(iωₙ)')
ax1.set_title(f'Self-energy (up spin)\nDLR fit residual: {residual["up"]:.2e}')
ax1.legend()
ax1.grid()

plt.tight_layout()
plt.show()
