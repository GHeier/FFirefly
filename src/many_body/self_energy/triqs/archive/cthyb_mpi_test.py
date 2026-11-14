"""
Minimal CTHYB test: single-band Anderson impurity model on a Bethe lattice
(semicircular DOS), solved self-consistently for a couple of DMFT iterations.
Sized to run in well under a minute with MPI.
"""
from triqs.gf import Gf, MeshImFreq, iOmega_n, inverse, SemiCircular
from triqs.operators import n
from triqs_cthyb import Solver
import triqs.utility.mpi as mpi

beta = 5.0
U = 2.0
mu = U / 2.0
half_bandwidth = 1.0
n_iw = 200
n_tau = 2001
n_l = 30

S = Solver(beta=beta, gf_struct=[('up', [0]), ('down', [0])],
           n_iw=n_iw, n_tau=n_tau, n_l=n_l)

# Initial guess: non-interacting semicircular bath
for name, g0 in S.G0_iw:
    g0 << inverse(iOmega_n + mu - (half_bandwidth / 2.0) ** 2 * SemiCircular(half_bandwidth))

n_dmft_loops = 2
for it in range(n_dmft_loops):
    if mpi.is_master_node():
        print(f"\n=== DMFT iteration {it} ===")

    S.solve(h_int=U * n('up', 0) * n('down', 0),
            n_cycles=2000,
            length_cycle=50,
            n_warmup_cycles=500,
            measure_G_l=True)

    # Bethe-lattice self-consistency: G0 = (iw + mu - t^2 G)^-1
    # G_iw is already all-reduced across ranks by the solver, so every rank
    # can update its own G0_iw identically without an extra broadcast.
    for name, g0 in S.G0_iw:
        g = S.G_iw[name]
        g0 << inverse(iOmega_n + mu - (half_bandwidth / 2.0) ** 2 * g)

if mpi.is_master_node():
    print("\nDone. Final G_iw (up) at first 5 Matsubara frequencies:")
    print(S.G_iw['up'].data[n_iw:n_iw + 5, 0, 0])
    from h5 import HDFArchive
    with HDFArchive("cthyb_test_output.h5", 'w') as ar:
        ar['G_iw'] = S.G_iw
        ar['G0_iw'] = S.G0_iw
        ar['Sigma_iw'] = S.Sigma_iw
    print("Saved results to cthyb_test_output.h5")
