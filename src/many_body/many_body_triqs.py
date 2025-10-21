from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.gf.meshes import MeshDLRImFreq
from triqs.dos import DOSFromFunction, HilbertTransform
from diagram import *
from load_triqs_H import *
from many_body_solver import ManyBodySolver
from ipt_solver import IPTSolver
import numpy as np
import firefly.config as cfg

def main():
    interaction = cfg.interaction
    if interaction == "FLEX":
        FLEX()
    elif interaction == "DMFT":
        DMFT()
    else:
        print(f"Interaction {interaction} not recognized")

def DMFT():
    H_r, kmesh, e_k = get_energy_mesh()

    # Compute DOS from energy dispersion
    energies = e_k.data.flatten().real
    eps_min, eps_max = energies.min(), energies.max()
    eps_range = eps_max - eps_min
    margin = 0.1 * eps_range

    # Create histogram for DOS
    hist, bin_edges = np.histogram(energies, bins=int(cfg.w_pts), density=True)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    # Create interpolated DOS function
    def dos_func(e):
        if e < eps_min - margin or e > eps_max + margin:
            return 0.0
        idx = np.searchsorted(bin_centers, e)
        if idx >= len(hist):
            idx = len(hist) - 1
        return hist[idx]

    # Create DOS object and Hilbert transform
    dos_obj = DOSFromFunction(function=dos_func, x_min=eps_min-margin, x_max=eps_max+margin, n_pts=int(cfg.w_pts))
    H = HilbertTransform(dos_obj)

    # Run IPT solver
    beta = 1.0 / cfg.Temperature
    S = IPTSolver(beta, H, n_loops=100, mix=0.10, tol=1e-6, w_max=1.2*4, eps=1e-14)
    S.loop(U=cfg.onsite_U)
    print(f"Final Sigma max: {np.max(np.abs(S.Sigma_iw.data)):.4f}")

def FLEX():
    # Get energy mesh (returns tuple: H_r, kmesh, e_k)
    H_r, kmesh, e_k = get_energy_mesh()

    # Set up DLR frequency mesh from config
    beta = 1.0 / cfg.Temperature
    w_max = 10.0  # DLR frequency cutoff
    eps = 1e-14   # DLR precision
    DLRImMesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)

    # Compute non-interacting Green's function
    G0 = lattice_dyson_g0_wk(mu=cfg.fermi_energy, e_k=e_k, mesh=DLRImMesh)

    # Initialize many-body solver with mixing parameter and e_k for mu calculation
    S = ManyBodySolver(G0, U=cfg.onsite_U, mix=0.2, U_maxiter=50, e_k=e_k)

    if cfg.self_consistent:
        S.loop_FLEX(n_loops=50, check_divergence=True)
    else:
        S.solve_FLEX()
    print(f"Final Max Chi: {np.max(np.abs(S.X.obj_wk.data)):.4f}")
    if hasattr(S, 'V'):
        print(f"Final Max Vertex: {np.max(np.abs(S.V.obj_wk.data)):.4f}")
    print(f"Final U*max(Chi): {S.UX:.4f}")
    print(f"Final U: {S.U:.4f}")

    save(S)

def save(S):
    outdir = cfg.outdir
    prefix = cfg.prefix
    pref = outdir + prefix
    S.G.save_as_w(pref + '_G_w.h5')
    S.X.save(pref + '_chi.h5')
    S.V.save(pref + '_vertex.h5')
    S.Sigma.save_as_w(pref + '_sigma_w.h5')

