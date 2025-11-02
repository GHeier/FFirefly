from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.dos import DOSFromFunction, HilbertTransform
from load_triqs_H import get_energy_mesh, create_dlr_meshes
from many_body_solver import ManyBodySolver
from ipt_solver import IPTSolver
from frequency_plots import plot_dmft_results, plot_frequency_data
import numpy as np
import firefly.config as cfg

n = cfg.num_electrons / 2
mu = cfg.fermi_energy

def main():
    interaction = cfg.interaction
    if interaction == "FLEX":
        FLEX()
    elif interaction == "DMFT":
        DMFT()
    elif interaction == "FLEX+DMFT":
        FLEX_DMFT()
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

    # Set up DMFT parameters
    beta = 1.0 / cfg.Temperature
    dlr_iw_mesh = create_dlr_meshes(e_k, beta, statistic='Fermion')
    G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
    G_iw << H(G_iw, mu=cfg.fermi_energy)
    init_n = G_iw.density().real[0][0]
    print(f"Initial chemical potential for n={init_n}: mu = {mu:.6f}")
    calculated_mu = calculate_mu_Gw(G_iw, n_target=n)
    print(f"Calculated chemical potential for n={n}: mu = {calculated_mu:.6f}")

    # Initialize ManyBodySolver in DMFT mode
    S = ManyBodySolver(G_iw, H=H, U=cfg.onsite_U, mix=0.10, n=n, mu=calculated_mu)

    # Run DMFT loop
    S.loop_DMFT(n_loops=100, tol=1e-6)

    # Print results
    print(f"Final Sigma max: {np.max(np.abs(S.Sigma_loc.obj_w.data)):.4f}")
    print(f"Final G max: {np.max(np.abs(S.G_loc.obj_w.data)):.4f}")

    # Save data
    save_DMFT(S)

    # Plot results
    #plot_dmft_results(S, save_dir=cfg.outdir)

def FLEX():
    # Get energy mesh (returns tuple: H_r, kmesh, e_k)
    H_r, kmesh, e_k = get_energy_mesh()

    # Set up DLR frequency mesh from config
    beta = 1.0 / cfg.Temperature
    DLRImMesh = create_dlr_meshes(e_k, beta, statistic='Fermion')

    # Compute non-interacting Green's function
    G0 = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)
    #calculated_mu = calculate_mu_Gk(G0, n_target=n)

    # Initialize many-body solver with mixing parameter and e_k for mu calculation
    S = ManyBodySolver(G0, U=cfg.onsite_U, mix=0.2, U_maxiter=50, n=n, mu=mu)
    new_mu = S.find_mu_for_density(n)

    G0 = lattice_dyson_g0_wk(mu=new_mu, e_k=e_k, mesh=DLRImMesh)
    S = ManyBodySolver(G0, U=cfg.onsite_U, mix=0.2, U_maxiter=100, n=n, mu=new_mu)


    if cfg.self_consistent:
        S.loop_FLEX(n_loops=50, check_divergence=True)
    else:
        S.solve_FLEX()
    if S.diverged:
        print("FLEX calculation diverged due to magnetic instability. Results are unreliable.")
        exit()
    print(f"Final Max Chi: {np.max(np.abs(S.X.obj_wk.data)):.4f}")
    if hasattr(S, 'V'):
        print(f"Final Max Vertex: {np.max(np.abs(S.V.obj_wk.data)):.4f}")
    print(f"Final U*max(Chi): {S.UX:.4f}")
    print(f"Final U: {S.U:.4f}")

    save_FLEX(S)

def FLEX_DMFT():
    # Get energy mesh (returns tuple: H_r, kmesh, e_k)
    H_r, kmesh, e_k = get_energy_mesh()

    # Set up DLR frequency mesh from config
    beta = 1.0 / cfg.Temperature
    DLRImMesh = create_dlr_meshes(e_k, beta, statistic='Fermion')

    # Compute non-interacting Green's function
    G0 = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)

    # Initialize many-body solver with mixing parameter and e_k for mu calculation
    S = ManyBodySolver(G0, U=cfg.onsite_U, mix=0.10, n=n, mu=mu)
    new_mu = S.find_mu_for_density(n)
    G0 = lattice_dyson_g0_wk(mu=new_mu, e_k=e_k, mesh=DLRImMesh)
    S = ManyBodySolver(G0, U=cfg.onsite_U, mix=0.10, n=n, mu=new_mu)

    if cfg.self_consistent:
        S.loop_FLEX_DMFT(n_loops=50, check_divergence=True)
    else:
        S.solve_FLEX_DMFT()
    if S.diverged:
        print("FLEX calculation diverged due to magnetic instability. Results are unreliable.")
        exit()
    print(f"Final Max Chi: {np.max(np.abs(S.X.obj_wk.data)):.4f}")
    if hasattr(S, 'V'):
        print(f"Final Max Vertex: {np.max(np.abs(S.V.obj_wk.data)):.4f}")
    print(f"Final U*max(Chi): {S.UX:.4f}")
    print(f"Final U: {S.U:.4f}")

    save_FLEX(S)

def save_FLEX(S):
    outdir = cfg.outdir
    prefix = cfg.prefix
    pref = outdir + prefix
    S.G.save_as_w(pref + '_G_w.h5')
    S.G.save(pref + '_G.h5')
    S.X.save(pref + '_chi.h5')
    S.V.save(pref + '_vertex.h5')
    S.Sigma.save(pref + '_sigma.h5')
    S.Sigma.save_as_w(pref + '_sigma_w.h5')


def save_DMFT(S):
    outdir = cfg.outdir
    prefix = cfg.prefix
    pref = outdir + prefix
    S.G_loc.save(pref + '_G_iw.h5')
    S.Sigma_loc.save(pref + '_sigma_iw.h5')
    S.G_loc.save_spectral(pref + '_A_w.h5')

def calculate_mu_Gw(G0, n_target):
    diff = 1.0
    G = G0.copy()
    old_mu = mu
    new_mu = mu
    while abs(diff) > 1e-6:
        G << inverse(inverse(G) + new_mu - old_mu)
        old_mu = new_mu
        n = G.density().real[0][0]
        diff = n - n_target
        #print(f"mu: {old_mu:.6f}, n: {n:.6f}, diff: {diff:.6f}")
        new_mu += -diff * 0.5
        
    return old_mu

def calculate_mu_Gk(G0, n_target):
    diff = 1.0
    G = G0.copy()
    print(G)
    n = G.density()
    print(n)
    old_mu = mu
    new_mu = mu
    while abs(diff) > 1e-6:
        G << inverse(inverse(G) + new_mu - old_mu)
        old_mu = new_mu
        n = G.density().real[0][0]
        diff = n - n_target
        print(f"mu: {old_mu:.6f}, n: {n:.6f}, diff: {diff:.6f}")
        new_mu += -diff * 0.5
        
    return mu
