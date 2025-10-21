from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.gf.meshes import MeshDLRImFreq
from diagram import *
from load_triqs_H import *
from many_body_solver import ManyBodySolver
import numpy as np
import firefly.config as cfg

def main():
    # Get energy mesh (returns tuple: H_r, kmesh, e_k)
    e_k = get_energy_mesh()

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
        S.loop(n_loops=50, check_divergence=True)
    else:
        S.solve()
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
    S.G.save_w(pref + '_G_w.dat')

