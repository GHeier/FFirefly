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
    X0 = S.X.copy()

    if cfg.self_consistent:
        S.loop_FLEX(n_loops=max_iters, check_divergence=True)
    else:
        S.solve_FLEX()
        #S.V = S.FLEX_from_chi(S.X.obj_wk)
        #S.Sigma_from_vertex()
        #S.G.obj_wk << inverse(inverse(S.G0.obj_wk) - S.Sigma.obj_wk)
        S.X = X0
    if S.diverged:
        pref = outdir + prefix
        #S.G.save_as_w(pref + '_G_w.h5')
        S.G0.save(pref + '_G0.h5')
        S.G.save(pref + '_G.h5')
        S.X.save(pref + '_chi.h5')
        print("FLEX calculation diverged due to magnetic instability. Results are unreliable.")
        exit()

    loc_sigma = make_local(S.Sigma.obj_wk.data)
    renorm = get_renorm(loc_sigma, S.Sigma.w_points)
    print(f"Quasiparticle renormalization factor: {renorm:.4f}")

    print(f"Final Max Chi: {np.max(np.abs(S.X.obj_wk.data)):.4f}")
    if hasattr(S, 'V'):
        print(f"Final Max Vertex: {np.max(np.abs(S.V.obj_wk.data)):.4f}")
    print(f"Final U*max(Chi): {S.UX:.4f}")
    print(f"Final U: {S.U:.4f}")
    final_mu = S.find_mu_for_density(n)
    print(f"Target n = {n:.6f}")
    print(f"Final mu = {final_mu:.6f}")

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
    S.G0.save(pref + '_G0.h5')
    S.G.save(pref + '_G.h5')
    S.X.save(pref + '_chi.h5')
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
    S.FLEX.X.save(pref + '_chi.h5')
    S.FLEX.V.save(pref + '_vertex.h5')
    S.Sigma_k.save(pref + '_sigma.h5')
    S.Sigma_k.save_as_w(pref + '_sigma_w.h5')
    S.Sigma_loc.save(pref + '_sigma_loc.h5')
    S.Sigma_nonloc.save_as_w(pref + '_sigma_nonloc.h5')
    S.Sigma_imp.save(pref + '_sigma_imp.h5')
    #S.Sigma.save_as_w(pref + '_sigma_w.h5')

def get_renorm(loc_sigma, w_points):
    loc_sigma = loc_sigma[:, 0, 0]
    signs = np.sign(loc_sigma.imag)
    diff = np.diff(signs)
    zero_crossings = np.where(diff != 0)[0]
    ind = zero_crossings[0]
    w_prev = w_points[ind]
    w_next = w_points[ind+1]
    sigma_prev = loc_sigma[ind]
    sigma_next = loc_sigma[ind+1]
    renorm = 1.0 - (sigma_next.imag - sigma_prev.imag) / (w_next - w_prev)
    print("w prev, next: ", w_prev, w_next)
    print("Sigma prev, next: ", sigma_prev.imag, sigma_next.imag)
    return renorm
