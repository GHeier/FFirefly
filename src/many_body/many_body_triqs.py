from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.gf.meshes import MeshDLRImFreq
from diagram import *
from load_triqs_H import *
from many_body_solver import ManyBodySolver
import numpy as np
import firefly.config as cfg

def main():
    # Get energy mesh (returns tuple: H_r, kmesh, e_k, e_k_min, e_k_max)
    H_r, kmesh, e_k, e_k_min, e_k_max = get_energy_mesh()
    print(f"Energy band: [{e_k_min:.4f}, {e_k_max:.4f}]")

    # Set up DLR frequency mesh from config
    beta = 1.0 / cfg.Temperature
    w_max = 10.0  # DLR frequency cutoff
    eps = 1e-14   # DLR precision
    DLRImMesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)

    # Determine chemical potential
    # If num_electrons is specified and different from default, calculate mu from density
    if hasattr(cfg, 'num_electrons') and cfg.num_electrons != 1.0:
        print(f"Target electron density: {cfg.num_electrons:.4f}")
        # For now, use fermi_energy as initial guess
        # The mu calculation will be done within the solver if needed
        mu = cfg.fermi_energy
        use_density_constraint = True
    else:
        mu = cfg.fermi_energy
        use_density_constraint = False

    # Compute non-interacting Green's function
    G0 = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)

    # Initialize many-body solver with mixing parameter
    S = ManyBodySolver(G0, U=cfg.onsite_U, mix=0.2, U_maxiter=50)

    # Store energy band info for mu calculation
    S.e_k_min = e_k_min
    S.e_k_max = e_k_max

    # Optionally calculate mu from target density
    if use_density_constraint:
        print("Calculating chemical potential from target electron density...")
        mu_calculated = S.mu_from_density(cfg.num_electrons, e_k_min, e_k_max)
        print(f"Using mu = {mu_calculated:.4f}")
        # Recompute G0 with new mu
        G0 = lattice_dyson_g0_wk(mu=mu_calculated, e_k=e_k, mesh=DLRImMesh)
        S = ManyBodySolver(G0, U=cfg.onsite_U, mix=0.2, U_maxiter=50)
        S.e_k_min = e_k_min
        S.e_k_max = e_k_max

    S.solve()
    print(f"Initial Max Chi: {np.max(np.abs(S.X.obj_wk.data)):.4f}")
    if hasattr(S, 'V'):
        print(f"Initial Max Vertex: {np.max(np.abs(S.V.obj_wk.data)):.4f}")
    print(f"Initial U*max(Chi): {S.UX:.4f}")

    # Self-consistent loop with divergence checking
    S.loop(n_loops=50, check_divergence=True)
    print(f"Final Max Chi: {np.max(np.abs(S.X.obj_wk.data)):.4f}")
    if hasattr(S, 'V'):
        print(f"Final Max Vertex: {np.max(np.abs(S.V.obj_wk.data)):.4f}")
    print(f"Final U*max(Chi): {S.UX:.4f}")
    print(f"Final U: {S.U:.4f}")
