from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.dos import DOSFromFunction, HilbertTransform
from load_triqs_H import get_energy_mesh, create_dlr_meshes
from IPTSolver import *
from FLEXSolver import *
from FLEX_DMFT_Solver import *
#from ipt_solver import IPTSolver
from frequency_plots import plot_dmft_results, plot_frequency_data
import numpy as np
import firefly as fly
import firefly.config as cfg

outdir = cfg.outdir
prefix = cfg.prefix

n = cfg.num_electrons / 2
mu = cfg.fermi_energy
beta = 1.0 / cfg.Temperature
w_pts = cfg.w_pts
mixing = cfg.mixing
U = cfg.onsite_U
max_iters = cfg.max_iters

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
    #H_r, kmesh, e_k = get_energy_mesh()

    ## Compute DOS from energy dispersion
    #energies = e_k.data.flatten().real
    #eps_min, eps_max = energies.min(), energies.max()

    ## Create histogram for DOS
    #hist, bin_edges = np.histogram(energies, bins=int(w_pts), density=True)
    #bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    ## Create interpolated DOS function
    #def dos_func(e):
    #    if e < eps_min - margin or e > eps_max + margin:
    #        return 0.0
    #    idx = np.searchsorted(bin_centers, e)
    #    if idx >= len(hist):
    #        idx = len(hist) - 1
    #    return hist[idx]

    dos_field = fly.Field_R(outdir + prefix + '_DOS.h5')
    data = dos_field.get_data()
    eps_min = np.min(data[:,0])
    eps_max = np.max(data[:,0])
    eps_range = eps_max - eps_min
    margin = 0.1 * eps_range
    def dos_func(e):
        #print("e: ", e, type(e))
        #print("dos: ", dos_field(e))
        return dos_field(float(e))
    # Create DOS object and Hilbert transform
    dos_obj = DOSFromFunction(function=dos_func, x_min=eps_min-margin, x_max=eps_max+margin, n_pts=int(w_pts))
    H = HilbertTransform(dos_obj)

    # Set up DMFT parameters

    # Initialize ManyBodySolver in DMFT mode
    S = IPTSolver(beta, H=H, mix=mixing, mu=mu, n_loops=max_iters)

    # Run DMFT loop
    S.loop(U)

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
    DLRImMesh = create_dlr_meshes(e_k, beta, statistic='Fermion')

    # Compute non-interacting Green's function
    G0 = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)
    print("G0 mesh size: ", G0.data.shape)

    # Initialize many-body solver with mixing parameter and e_k for mu calculation
    S = FLEXSolver(G0, U=U, mix=mixing, U_maxiter=max_iters, n=n, mu=mu)

    if cfg.self_consistent:
        S.loop_FLEX(n_loops=max_iters, check_divergence=True)
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
    DLRImMesh = create_dlr_meshes(e_k, beta, statistic='Fermion')

    # Compute non-interacting Green's function
    G0 = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)

    # Initialize many-body solver with mixing parameter and e_k for mu calculation
    S = FLEX_DMFT_Solver(G0, U=U, mix=mixing, n=n, mu=mu)

    if cfg.self_consistent:
        S.loop_FLEX_DMFT(n_loops=cfg.max_iters, check_divergence=True)
    else:
        S.solve_FLEX_DMFT(S.FLEX, S.IPT)
    if S.diverged:
        print("FLEX calculation diverged due to magnetic instability. Results are unreliable.")
        exit()
    if hasattr(S, 'V'):
        print(f"Final Max Chi: {np.max(np.abs(S.FLEX.X.obj_wk.data)):.4f}")
        print(f"Final Max Vertex: {np.max(np.abs(S.FLEX.V.obj_wk.data)):.4f}")
    print(f"Final U*max(Chi): {S.FLEX.UX:.4f}")
    print(f"Final U: {S.U:.4f}")

    save_FLEX_DMFT(S)

def save_FLEX(S):
    pref = outdir + prefix
    #S.G.save_as_w(pref + '_G_w.h5')
    S.G.save(pref + '_G.h5')
    #S.X.save(pref + '_chi.h5')
    S.V.save(pref + '_vertex.h5')
    S.Sigma.save(pref + '_sigma.h5')
    #S.Sigma.save_as_w(pref + '_sigma_w.h5')

    # Compute and save singlet pairing vertex for superconductivity
    #V_singlet = S.compute_V_singlet()
    #V_singlet.save(pref + '_vertex_singlet.h5')

def save_FLEX_DMFT(S):
    pref = outdir + prefix
    #S.FLEX.G.save_as_w(pref + '_G_w.h5')
    S.FLEX.G.save(pref + '_G.h5')
    #S.X.save(pref + '_chi.h5')
    S.FLEX.V.save(pref + '_vertex.h5')
    S.Sigma_k.save(pref + '_sigma.h5')
    S.Sigma_k.save_as_w(pref + '_sigma_w.h5')
    S.Sigma_loc.save(pref + '_sigma_loc.h5')
    S.Sigma_nonloc.save_as_w(pref + '_sigma_nonloc.h5')
    S.Sigma_imp.save(pref + '_sigma_imp.h5')
    #S.Sigma.save_as_w(pref + '_sigma_w.h5')

def save_DMFT(S):
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

def test():
    pref = outdir + prefix
    dos_field = fly.Field_R(outdir + prefix + '_DOS.h5')
    data = dos_field.get_data()
    eps_min = np.min(data[:,0])
    eps_max = np.max(data[:,0])
    eps_range = eps_max - eps_min
    margin = 0.1 * eps_range
    def dos_func(e):
        #print("e: ", e, type(e))
        #print("dos: ", dos_field(e))
        return dos_field(float(e))
    # Create DOS object and Hilbert transform
    dos_obj = DOSFromFunction(function=dos_func, x_min=eps_min-margin, x_max=eps_max+margin, n_pts=int(w_pts))
    H = HilbertTransform(dos_obj)

    # Set up DMFT parameters

    # Initialize ManyBodySolver in DMFT mode
    S = IPTSolver(beta, H=H, mix=mixing, mu=mu, n_loops=max_iters)
    sigma = S.get_IPT_Sigma(U)
    sigma.save(pref + '_sigma_imp.h5')

    H_r, kmesh, e_k = get_energy_mesh()

    # Set up DLR frequency mesh from config
    DLRImMesh = create_dlr_meshes(e_k, beta, statistic='Fermion')

    # Compute non-interacting Green's function
    G0 = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)

    # Initialize many-body solver with mixing parameter and e_k for mu calculation
    S = FLEX_DMFT_Solver(G0, U=U, mix=mixing, n=n, mu=mu)
    S.FLEX.G = S.G0.copy()
    S.FLEX.get_local_G()
    S.IPT.G_loc = S.FLEX.G_loc.copy()
    S.IPT.set_Weiss()
    sigma = S.IPT.get_IPT_Sigma(U)
    sigma.save(pref + '_sigma_loc.h5')

