import firefly as fly
import firefly.config as cfg
import numpy as np

# Load relevant variables from the configuration
dim = cfg.dimension
mu = cfg.fermi_energy
w_pts = cfg.w_pts
wc = cfg.cutoff_energy
mu_from_n = cfg.mu_from_n
n = cfg.num_electrons
outdir = cfg.outdir
prefix = cfg.prefix
U = cfg.U0
Z = cfg.qp_weight

if mu_from_n:
    print(f"Initial mu = {mu}")
    En = fly.Field_R(outdir + prefix + "_E_vs_n.h5")
    mu = En(n)
    print(f"Shifted mu = {mu}")


def load_surface():
    """Load Fermi surface k-points and areas using the Bands dispersion."""
    band = fly.Bands()
    # Callback receives a Vec ctypes.Structure with x, y, z fields
    eps_func = lambda k: band(1, [k.x, k.y, k.z])
    surf = fly.Surface(eps_func, mu)
    kpoints, areas = surf.get_faces_and_areas()
    print(f"Found {len(kpoints)} k-points on Fermi surface")
    print(f"Total area: {sum(areas)}")
    return kpoints, areas


def get_surface_data():
    """Generate k-points, DOS weights, and frequency points for Fermi surface calculations."""
    print(f"\nGenerating k-points from Fermi surface at μ = {mu}")
    H = fly.Hamiltonian()

    # Load Fermi surface
    kpoints, areas = load_surface()
    n_k = len(kpoints)
    print(f"Total k-points: {n_k}")

    # Compute Fermi velocities
    print("\nComputing Fermi velocities...")
    velocities = H.get_fermi_velocity(kpoints)  # Shape: (n_k, n_bands, 3)
    v_norms = np.array([np.linalg.norm(velocities[i, 0, :]) for i in range(n_k)])
    print(f"Fermi Velocity Min, Max, Ave: {v_norms.min():.6f}, {v_norms.max():.6f}, {v_norms.mean():.6f}")

    # Compute DOS weights: area / |v_F| / (2π)^dim
    areas_arr = np.array(areas, dtype=np.float32)
    dos_weights = areas_arr / v_norms / (2 * np.pi)**dim
    print(f"DOS Min, Max, Ave, Total: {dos_weights.min():.6f}, {dos_weights.max():.6f}, {dos_weights.mean():.6f}, {dos_weights.sum():.6f}")

    return kpoints, dos_weights


def get_scattering_kpts(kpoints):
    """Compute 2D array of scattering vectors q = k - k' for all pairs (k, k')."""
    n_k = len(kpoints)
    kpts_arr = np.array(kpoints, dtype=np.float32)
    # q[i, j] = k[i] - k[j]
    q_kk = kpts_arr[:, np.newaxis, :] - kpts_arr[np.newaxis, :, :]
    return q_kk


def run():
    """Main entry point for FS_approx renormalization calculation."""
    kpts, dA = get_surface_data()
    dos = np.sum(dA)
    npts = len(kpts)
    kkpts = get_scattering_kpts(kpts).reshape(-1, dim)

    chi = fly.Field_R(outdir + prefix + "_chi.h5")
    V = U**1 * chi(kkpts).reshape(npts, npts)

    val = dA @ V @ dA / dos
    print("lambda_z = ", val)
    print("Quasiparticle Weight: ", 1 / (1+val))


    return val


if __name__ == "__main__":
    run()


