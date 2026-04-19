import firefly as fly
import firefly.config as cfg
from firefly.diagram import Diagram
import numpy as np
from triqs.dos import DOSFromFunction, HilbertTransform
from triqs.gf import Gf, inverse, SemiCircular
from triqs_tprf.lattice import lattice_dyson_g0_wk
from IPTSolver import IPTSolver
from IPTSolver_real import IPTSolver_real
from BubbleSolver import BubbleSolver
from Bubble_DMFTSolver import Bubble_DMFTSolver
from load_triqs_H import create_dlr_meshes, get_energy_mesh
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime

# Load relevant variables from the configuration
outdir = cfg.outdir
prefix = cfg.prefix
debug = cfg.debug

interaction = cfg.interaction
mu = cfg.fermi_energy
n = cfg.num_electrons
print(f"mu = {mu}, n = {n}")
if cfg.mu_from_n and not debug:
    En = fly.Field_R(outdir + prefix + "_E_vs_n.h5")
    mu = En(n)
    print(f"Shifted mu to {mu}")

T = cfg.Temperature
beta = 1 / T
w_pts = cfg.w_pts
mixing = cfg.mixing
U = cfg.U0  # Hubbard U parameter
max_iters = cfg.max_iters

nx, ny, nz = cfg.k_mesh
BZ = np.array(cfg.brillouin_zone)
dim = cfg.dimension
if dim == 2:
    nz = 1

eps = 1e-14

def debug_p(e):
    return 1/ (2 * np.pi) * (4 - e**2)**0.5

def get_H_debug():
    w_pts = 1000
    eps_min = -2.0
    eps_max = 2.0
    eps_range = eps_max - eps_min
    print("Emin, Emax: ", eps_min, eps_max)
    print(f"DOS integrated with {w_pts} w-points")
    margin = 0.1 * eps_range

    def dos_func(e):
        return debug_p(float(e))

    dos_obj = DOSFromFunction(function=dos_func, x_min=eps_min, x_max=eps_max, n_pts=int(w_pts))
    H = HilbertTransform(dos_obj)
    return H, eps_range

def get_H(N):
    data = N.w_points
    eps_min = np.min(data)
    eps_max = np.max(data)
    eps_range = eps_max - eps_min
    print("Emin, Emax: ", eps_min, eps_max)
    print(f"DOS integrated with {w_pts} w-points")
    margin = 0.1 * eps_range

    def dos_func(e):
        return N(float(e))

    dos_obj = DOSFromFunction(function=dos_func, x_min=eps_min, x_max=eps_max, n_pts=int(w_pts))
    H = HilbertTransform(dos_obj)
    return H, eps_range


def run():
    print(f"mu = {mu}, n = {n}")
    print(f"interaction = {interaction}")
    print("mixing = ", mixing)

    if interaction == "DMFT":
        return run_DMFT()
    elif interaction == "Bubble":
        return run_Bubble()
    elif interaction == "Bubble+DMFT":
        return run_Bubble_DMFT()
    else:
        raise ValueError(f"Unknown interaction: {interaction}. Use 'DMFT', 'Bubble', or 'Bubble+DMFT'.")


def get_G0_wk():
    """Create G0(k,iw) from tight-binding model."""
    H_r, kmesh, e_k = get_energy_mesh()
    emax = e_k.data.max().real
    emin = e_k.data.min().real
    D = 1.2 * (emax - emin)
    DLRImMesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=D, eps=eps)
    G0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)
    return G0_wk, D


def run_DMFT():
    """Run DMFT with IPT solver."""
    if debug:
        H, D = get_H_debug()
    else:
        N = fly.Field_R(outdir + prefix + '_DOS.h5')
        H, D = get_H(N)

    if T == 0.0:
        print("Temperature = 0: Using real-axis IPT solver")
        S = IPTSolver_real(H=H, mu=mu, mix=mixing, n_loops=max_iters,
                          w_min=-D, w_max=D, n_w=int(w_pts))
        S.loop(U)
        print(f"Final Sigma max: {np.max(np.abs(S.Sigma_loc.data)):.4f}")
        print(f"Final G max: {np.max(np.abs(S.G_loc.data)):.4f}")
        renorm = get_renorm_real(S.Sigma_loc.data, S.w_points)
        save_DMFT_real(S)
    else:
        print(f"Temperature = {T}: Using Matsubara IPT solver")
        S = IPTSolver(beta, H=H, mix=mixing, mu=mu, n_loops=max_iters, w_max=1.2*D, eps=eps)
        S.loop(U, bethe_lattice=False)
        print(f"Final Sigma max: {np.max(np.abs(S.Sigma_loc.obj_w.data)):.4f}")
        print(f"Final G max: {np.max(np.abs(S.G_loc.obj_w.data)):.4f}")
        renorm = get_renorm(S.Sigma_loc.obj_w.data, S.Sigma_loc.w_points)
        save_DMFT(S, D)

    print(f"Quasiparticle Weight: {1/renorm:.4f}")
    print(f"m*/m: {renorm:.4f}")
    print(f"lambda_z: {renorm - 1:.4f}")
    return 1/renorm


def run_Bubble():
    """Run Bubble (second-order perturbation theory) solver."""
    G0_wk, D = get_G0_wk()
    print(f"Temperature = {T}: Using Bubble solver")

    S = BubbleSolver(G0=G0_wk, U=U, mix=mixing, n=n, mu=mu)
    S.loop_Bubble(n_loops=max_iters)

    print(f"Final Sigma max: {np.max(np.abs(S.Sigma.obj_wk.data)):.4f}")
    print(f"Final G max: {np.max(np.abs(S.G.obj_wk.data)):.4f}")
    # Get local Sigma for renorm calculation
    Sigma_loc = np.einsum('wknm->wnm', S.Sigma.obj_wk.data) / S.G.nk
    w_points = S.G_loc.w_points
    renorm = get_renorm(Sigma_loc, w_points)
    save_Bubble(S, D)

    print(f"Quasiparticle Weight: {1/renorm:.4f}")
    print(f"m*/m: {renorm:.4f}")
    print(f"lambda_z: {renorm - 1:.4f}")
    return 1/renorm


def run_Bubble_DMFT():
    """Run Bubble+DMFT solver."""
    G0_wk, D = get_G0_wk()
    print(f"Temperature = {T}: Using Bubble+DMFT solver")

    S = Bubble_DMFTSolver(G0=G0_wk, U=U, mix=mixing, n=n, mu=mu)
    S.loop_Bubble_DMFT(n_loops=max_iters)

    print(f"Final Sigma max: {np.max(np.abs(S.Bubble.Sigma.obj_wk.data)):.4f}")
    print(f"Final G max: {np.max(np.abs(S.Bubble.G.obj_wk.data)):.4f}")
    # Get local Sigma for renorm calculation
    Sigma_loc = np.einsum('wknm->wnm', S.Sigma_k.obj_wk.data) / S.Bubble.G.nk
    w_points = S.Bubble.G_loc.w_points
    renorm = get_renorm(Sigma_loc, w_points)
    save_Bubble_DMFT(S, D)

    print(f"Quasiparticle Weight: {1/renorm:.4f}")
    print(f"m*/m: {renorm:.4f}")
    print(f"lambda_z: {renorm - 1:.4f}")
    return 1/renorm

def save_DMFT(S, eps_range=0.0):
    pref = outdir + prefix
    S.G_loc.save(pref + '_G_iw.h5')
    S.Sigma_loc.save(pref + '_self_energy.h5')
    S.Sigma_loc.save(pref + '_sigma_iw.h5')
    S.G_loc.save_spectral(pref + '_A_w.h5')
    save_G(S, eps_range)


def save_Bubble(S, eps_range=0.0):
    pref = outdir + prefix
    S.G_loc.save(pref + '_G_iw.h5')
    S.Sigma.save(pref + '_self_energy.h5')
    S.Sigma.save(pref + '_sigma_wk.h5')
    S.G.save(pref + '_G.h5')
    S.X.save(pref + '_chi0.h5')


def save_Bubble_DMFT(S, eps_range=0.0):
    pref = outdir + prefix
    S.Bubble.G_loc.save(pref + '_G_iw.h5')
    S.Bubble.Sigma.save(pref + '_self_energy.h5')
    S.Bubble.Sigma.save(pref + '_sigma_wk.h5')
    S.Bubble.G.save(pref + '_G.h5')
    S.Sigma_imp.save(pref + '_sigma_imp.h5')
    S.Sigma_nonloc.save(pref + '_sigma_nonloc.h5')

def save_G(S, eps_range):
    pref = outdir + prefix
    H_r, kmesh, e_k = get_energy_mesh()
    DLRImMesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=1.2*eps_range, eps=eps)
    G = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)
    G = Diagram(G, 'Fermion')
    G.obj_wk = inverse(G.obj_wk)
    G.obj_wk.data[:] = G.obj_wk.data[:] - S.Sigma_loc.obj_w.data[:, np.newaxis]
    G.obj_wk = inverse(G.obj_wk)
    G.save(pref + '_G.h5')


def save_DMFT_real(S):
    """Save real-axis DMFT results."""
    pref = outdir + prefix
    # For real-axis solver, use built-in save methods
    S.save_results(pref)


def get_renorm(loc_sigma, w_points):
    loc_sigma = loc_sigma[:, 0, 0]
    signs = np.sign(loc_sigma.imag)
    diff = np.diff(signs)
    zero_crossings = np.where(diff != 0)[0]

    if (len(zero_crossings) == 0):
        print("Uncontrolled Self-Energy result, no zero crossing. Returning infinity")
        return float('inf')

    ind = zero_crossings[0]
    w_prev = w_points[ind]
    w_next = w_points[ind+1]
    sigma_prev = loc_sigma[ind]
    sigma_next = loc_sigma[ind+1]
    renorm = 1.0 - (sigma_next.imag - sigma_prev.imag) / (w_next - w_prev)
    print("w prev, next: ", w_prev, w_next)
    print("Sigma prev, next: ", sigma_prev.imag, sigma_next.imag)

    return renorm


def get_renorm_real(loc_sigma, w_points):
    """Calculate quasiparticle renormalization factor for real-axis Green's function."""
    loc_sigma = loc_sigma[:, 0, 0]
    # For real-axis, find zero crossing in imaginary part
    signs = np.sign(loc_sigma.imag)
    diff = np.diff(signs)
    zero_crossings = np.where(diff != 0)[0]
    if len(zero_crossings) == 0:
        print("Warning: No zero crossing found in Im[Sigma]. Using derivative at w=0.")
        # Find index closest to w=0
        ind = np.argmin(np.abs(w_points))
        if ind == 0 or ind == len(w_points) - 1:
            return 1.0  # Can't calculate derivative at boundary
        w_prev = w_points[ind-1]
        w_next = w_points[ind+1]
        sigma_prev = loc_sigma[ind-1]
        sigma_next = loc_sigma[ind+1]
    else:

        ind = zero_crossings[0]
        w_prev = w_points[ind]
        w_next = w_points[ind+1]
        sigma_prev = loc_sigma[ind]
        sigma_next = loc_sigma[ind+1]


    renorm = 1.0 - (sigma_next.imag - sigma_prev.imag) / (w_next - w_prev)
    print("w prev, next: ", w_prev, w_next)
    print("Sigma prev, next: ", sigma_prev.imag, sigma_next.imag)
    return renorm

if __name__ == "__main__": # Runs on file execution
    run()


