import firefly as fly
import firefly.config as cfg
import numpy as np
from scipy.optimize import linear_sum_assignment

# Load relevant variables from the configuration
outdir = cfg.outdir
prefix = cfg.prefix

kmesh = cfg.k_mesh
nx, ny, nz = kmesh
dim = cfg.dimension
if dim == 2:
    nz = 1

nbnd = cfg.nbnd
BZ = np.array(cfg.brillouin_zone)
mu = cfg.fermi_energy
n = cfg.num_electrons
if cfg.mu_from_n:
    N = fly.Field_R(outdir + prefix + "_E_vs_n.h5")
    mu = N(n)

def get_k_mesh():
    kx, ky, kz = np.meshgrid(np.arange(nx) / nx,
                             np.arange(ny) / ny,
                             np.arange(nz) / nz,
                             indexing='ij')
    kpts = np.stack([kx.ravel(), ky.ravel(), kz.ravel()], axis=-1)
    kpts = kpts @ BZ.T
    return kpts

def save_hamiltonian(H, kpts):
    """Save the raw Hamiltonian matrix on the k-mesh as {prefix}_Hk.h5,
    so fly.Bands() can load and diagonalize it directly (it looks for
    this exact filename and falls back to a single analytic band if
    it isn't found)."""
    Hk = np.stack(H(list(kpts)))  # shape (Nk, n, n), complex
    n = Hk.shape[-1]
    fly.save_data(outdir + prefix + '_Hk.h5', Hk, [nx, ny, nz], BZ, inds=[n, n], centered=False)
    print("Saved to ", outdir + prefix + "_Hk.h5")

def isolate_bands():
    H = fly.Hamiltonian()
    kpts = get_k_mesh()
    nk = kpts.shape[0]

    # Diagonalize the raw Hamiltonian ourselves with numpy.linalg.eig instead of
    # H.get_wavefunctions()/get_bands() (LAPACK cheev with jobz='V'/'N'), which
    # silently re-sorts eigenvalues into ascending order at every k-point and
    # hides band crossings. eig (LAPACK geev) does not sort, so band identity
    # has to be tracked explicitly via eigenvector continuity instead.
    Hk = np.stack(H(list(kpts)))  # (Nk, n, n), complex
    total_bnd = Hk.shape[-1]
    print("Number of bands in Hamiltonian: ", total_bnd)

    ek_raw, evecs_raw = np.linalg.eig(Hk)  # ek_raw: (Nk, n) complex, evecs_raw: (Nk, n, n)
    ek_raw = ek_raw.real

    ek_grid = ek_raw.reshape(nx, ny, nz, total_bnd).copy()
    ev_grid = evecs_raw.reshape(nx, ny, nz, total_bnd, total_bnd).copy()

    # BFS from Gamma=(0,0,0). At each k-point, find the permutation of bands that
    # maximizes overlap with the eigenvectors of an already-assigned neighbor
    # (the standard diabatic/continuity tracking fix for unsorted eigensolvers).
    from collections import deque

    assigned = np.zeros((nx, ny, nz), dtype=bool)
    assigned[0, 0, 0] = True

    queue = deque()
    for dix, diy, diz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
        nix, niy, niz = dix, diy, diz
        if nix < nx and niy < ny and niz < nz:
            queue.append((nix, niy, niz))

    while queue:
        ix, iy, iz = queue.popleft()
        if assigned[ix, iy, iz]:
            continue

        ref_evecs = None
        ref_ek = None
        for dix, diy, diz in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
            nix, niy, niz = ix + dix, iy + diy, iz + diz
            if 0 <= nix < nx and 0 <= niy < ny and 0 <= niz < nz and assigned[nix, niy, niz]:
                ref_evecs = ev_grid[nix, niy, niz]
                ref_ek = ek_grid[nix, niy, niz]
                break
        if ref_evecs is None:
            queue.append((ix, iy, iz))
            continue

        curr_evecs = ev_grid[ix, iy, iz]
        curr_ek = ek_grid[ix, iy, iz]

        # At (near-)exact degeneracies, eig returns arbitrary eigenvectors within the
        # degenerate subspace, so overlap-based matching is unreliable for a single
        # k-point. Detect this case and fall back to matching by closest eigenvalue
        # to the reference instead, since energies stay continuous even there.
        degenerate = np.min(np.abs(curr_ek[:, None] - curr_ek[None, :]) + np.eye(total_bnd) * 1e9) < 0.05
        if degenerate:
            cost = np.abs(ref_ek[:, None] - curr_ek[None, :])
        else:
            cost = -np.abs(ref_evecs.conj().T @ curr_evecs)  # (n, n), rows=ref band, cols=curr band
        row_ind, col_ind = linear_sum_assignment(cost)
        perm = col_ind[np.argsort(row_ind)]  # perm[i] = which current column matches ref band i

        ek_grid[ix, iy, iz] = ek_grid[ix, iy, iz, perm]
        ev_grid[ix, iy, iz] = ev_grid[ix, iy, iz][:, perm]

        assigned[ix, iy, iz] = True
        for dix, diy, diz in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]:
            nix, niy, niz = ix + dix, iy + diy, iz + diz
            if 0 <= nix < nx and 0 <= niy < ny and 0 <= niz < nz and not assigned[nix, niy, niz]:
                queue.append((nix, niy, niz))

    ek = ek_grid.reshape(nk, total_bnd)
    ranks, nk = rank_bands_by_mu_distance(ek)
    ranks = ranks[:nbnd]
    ranks = np.sort(ranks)  # keep ascending energy order among the selected bands
    ek_copy = np.zeros((nk, nbnd))
    for i in range(nbnd):
        ek_copy[:, i] = ek[:, ranks[i]]
    return ek_copy

def rank_bands_by_mu_distance(ek):
    nk, total_bnd = ek.shape
    dists = [np.min(np.abs(ek[:, i] - mu)) for i in range(total_bnd)]
    ranks = np.argsort(dists)
    dists = [float(round(d, 3)) for d in sorted(dists)]
    print("Band minimum distances from fermi level: ", dists)
    if len(dists) >= nbnd and dists[nbnd-1] < 0.05:
        print("Bands not included are within 0.05 eV of Fermi level. Recommend increasing nbnd.")
    return ranks, nk
        

def run():
    ek = isolate_bands()
    print("Bands are tracked for continuity across crossings, then sorted by proximity to Fermi Surface")
    for n in range(1, nbnd+1):
        fly.save_data(outdir + prefix + '_bands_' + str(n) + '.h5', ek[:,n-1], [nx, ny, nz], BZ, inds=[], centered=False)
        print("Saved band #", n, " to ", outdir + prefix + "_bands_" + str(n) + ".h5")

    return np.min(ek[:,0] - mu) # Return something of any type that can be tested in the test suite.

if __name__ == "__main__": # Runs on file execution
    run()


