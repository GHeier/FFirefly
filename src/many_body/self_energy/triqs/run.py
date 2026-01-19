import firefly as fly
import firefly.config as cfg
import numpy as np
from triqs.dos import DOSFromFunction, HilbertTransform
from IPTSolver import IPTSolver
from IPTSolver_real import IPTSolver_real

# Load relevant variables from the configuration
outdir = cfg.outdir
prefix = cfg.prefix

interaction = cfg.interaction
mu = cfg.fermi_energy
T = cfg.Temperature
w_pts = cfg.w_pts
mixing = cfg.mixing
U = cfg.U0  # Hubbard U parameter
max_iters = cfg.max_iters

def run():
    print("mixing = ", mixing)
    # IPT solver works with any DOS-based approach
    # The interaction type (DMFT, FLEX, etc.) is primarily used in other parts of the many_body code
    if interaction not in ["DMFT", "FLEX"]:
        print(f"Warning: interaction = '{interaction}' may not be fully supported.")
        print("IPT solver will proceed using DOS-based local approximation.")

    N = fly.Field_R(outdir + prefix + '_DOS.h5')
    data = N.w_points
    eps_min = np.min(data)
    eps_max = np.max(data)
    eps_range = eps_max - eps_min
    print("Emin, Emax: ", eps_min, eps_max)
    print(f"DOS integrated with {w_pts} w-points")
    margin = 0.1 * eps_range
    def dos_func(e):
        #print("e: ", e, type(e))
        #print("dos: ", N(e))
        return N(float(e))
    # Create DOS object and Hilbert transform
    dos_obj = DOSFromFunction(function=dos_func, x_min=eps_min, x_max=eps_max, n_pts=int(w_pts))
    H = HilbertTransform(dos_obj)

    # Set up DMFT parameters

    # Check temperature and dispatch to appropriate solver
    if T == 0.0:
        print("Temperature = 0: Using real-axis IPT solver")
        S = IPTSolver_real(H=H, mu=mu, mix=mixing, n_loops=max_iters,
                          w_min=eps_min-margin, w_max=eps_max+margin, n_w=int(w_pts))
        # Run DMFT loop
        S.loop(U)

        # Print results
        print(f"Final Sigma max: {np.max(np.abs(S.Sigma_loc.data)):.4f}")
        print(f"Final G max: {np.max(np.abs(S.G_loc.data)):.4f}")
        renorm = get_renorm_real(S.Sigma_loc.data, S.w_points)
        print(f"Quasiparticle renormalization factor: {renorm:.4f}")

        # Save data
        save_DMFT_real(S)
    else:
        print(f"Temperature = {T}: Using Matsubara IPT solver")
        # Initialize ManyBodySolver in DMFT mode
        beta = 1.0 / cfg.Temperature
        S = IPTSolver(beta, H=H, mix=mixing, mu=mu, n_loops=max_iters)

        # Run DMFT loop
        S.loop(U)

        # Print results
        print(f"Final Sigma max: {np.max(np.abs(S.Sigma_loc.obj_w.data)):.4f}")
        print(f"Final G max: {np.max(np.abs(S.G_loc.obj_w.data)):.4f}")
        renorm = get_renorm(S.Sigma_loc.obj_w.data, S.Sigma_loc.w_points)
        print(f"Quasiparticle renormalization factor: {renorm:.4f}")

        # Save data
        save_DMFT(S)

    return renorm # Return something of any type that can be tested in the test suite.

def save_DMFT(S):
    pref = outdir + prefix
    S.G_loc.save(pref + '_G_iw.h5')
    S.Sigma_loc.save(pref + '_sigma_iw.h5')
    S.G_loc.save_spectral(pref + '_A_w.h5')


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


