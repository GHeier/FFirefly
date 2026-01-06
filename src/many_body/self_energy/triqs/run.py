import firefly as fly
import firefly.config as cfg
from . import IPTSolver

# Load relevant variables from the configuration
outdir = cfg.outdir
prefix = cfg.prefix

interaction = cfg.interaction
mu = cfg.fermi_energy
T = cfg.Temperature
w_pts = cfg.w_pts
mixing = cfg.mixing
U = cfg.onsite_U
max_iters = cfg.max_iters

def run():
    print("mixing = ", mixing)
    if interaction != "DMFT":
        print("Only DMFT interaction is implemented in this example.")
        exit(1)

    N = fly.Field_R(outdir + prefix + '_DOS.h5')
    data = N.get_data()
    eps_min = np.min(data[:,0])
    eps_max = np.max(data[:,0])
    eps_range = eps_max - eps_min
    margin = 0.1 * eps_range
    def dos_func(e):
        #print("e: ", e, type(e))
        #print("dos: ", N(e))
        return N(float(e))
    # Create DOS object and Hilbert transform
    dos_obj = DOSFromFunction(function=dos_func, x_min=eps_min-margin, x_max=eps_max+margin, n_pts=int(w_pts))
    H = HilbertTransform(dos_obj)

    # Set up DMFT parameters

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

if __name__ == "__main__": # Runs on file execution
    run()


