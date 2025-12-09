import firefly as fly
import firefly.config as cfg
from firefly.diagram import *

from triqs.gf import *
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.lattice import BrillouinZone, BravaisLattice
from triqs.gf.mesh_product import MeshProduct
from triqs_tprf.lattice import *
from triqs_tprf import *

import numpy as np
from numpy.fft import fftn, ifftn
from scipy.sparse.linalg import LinearOperator, eigsh


outdir = cfg.outdir
prefix = cfg.prefix

nstates = cfg.nstates
Nk = cfg.k_mesh[0]
nx, ny, nz = cfg.k_mesh
dim = cfg.dimension
if dim == 2:
    nz = 1
BZ = np.array(cfg.brillouin_zone)

mu = cfg.fermi_energy
beta = 1.0 / cfg.Temperature

def get_k_mesh(BZ, nx, ny, nz):
    # fractional grid (0→1)
    fx = np.linspace(0, 1, nx, endpoint=False)
    fy = np.linspace(0, 1, ny, endpoint=False)
    fz = np.linspace(0, 1, nz, endpoint=False)

    # 3D fractional mesh
    F = np.stack(np.meshgrid(fx, fy, fz, indexing='ij'), axis=-1)  # (nx,ny,nz,3)

    # convert to Cartesian k = f1*b1 + f2*b2 + f3*b3
    K = F @ BZ.T                          # shape (nx,ny,nz,3)

    return K.reshape(-1, 3)           # full grid, flattened grid

def multiply(V_r, f_k, f_r, D_k_flat, mesh, arrsize):
    """
    Perform convolution: result = IFFT[V(r) * FFT[D(k)]]

    Args:
        V_r: Vertex in real space, shape (nkx, nky, nkz, nstates, nstates, nstates, nstates)
        f_k: Square root of form factor in momentum space, shape (nkx, nky, nkz, 1, 1)
        f_r: Square root of form factor in real space, shape (nkx, nky, nkz, 1, 1)
        D_k_flat: Gap function in k-space, flattened, shape (nkx*nky*nkz*nstates*nstates,)
        mesh: (nkx, nky, nkz)
        arrsize: (nstates, nstates)

    Returns:
        Flattened result in k-space
    """
    nkx, nky, nkz = mesh
    nstates_a, nstates_b = arrsize

    # Reshape D_k_flat from (nk*nstates*nstates,) to (nkx, nky, nkz, nstates, nstates)
    D_k = D_k_flat.reshape(nkx, nky, nkz, nstates_a, nstates_b)
    D_r = fftn(D_k, axes=(0, 1, 2))
    #D_r_flip = D_r[::-1, ::-1, ::-1, :, :]
    #D_r_even = (D_r + D_r_flip) / 2.0

    #h_r = np.einsum('xyzabcd,xyzcd->xyzab', V_r, f_r * D_r_even)
    h_r = np.einsum('xyzabcd,xyzcd->xyzab', V_r, D_r)
    h_k = ifftn(h_r, axes=(0, 1, 2))
    return h_k.flatten()
    phi_k = h_k * f_k

    # Flatten back to 1D vector
    return phi_k.flatten()

def multiply_basic(V_r, f_r, D_k_flat, mesh, arrsize):
    nkx, nky, nkz = mesh
    nstates_a, nstates_b = arrsize
    D_k = D_k_flat.reshape(nkx, nky, nkz, nstates_a, nstates_b)
    D_r = fftn(D_k, axes=(0, 1, 2))
    D_r_flip = D_r[::-1, ::-1, ::-1, :, :]
    D_r_even = (D_r + D_r_flip) / 2.0

    #h_r = np.einsum('xyzabcd,xyzcd->xyzab', V_r, f_r * D_r_even)
    h_r = np.einsum('xyzabcd,xyzcd->xyzab', V_r, f_r * D_r)
    h_k = ifftn(h_r, axes=(0, 1, 2))
    return h_k.flatten()

def project_out(v, eigvecs):
    """
    Project out previously found eigenvectors using Gram-Schmidt orthogonalization.

    Args:
        v: Vector to project, shape (n,)
        eigvecs: List of previously found eigenvectors to project out

    Returns:
        Orthogonalized and normalized vector
    """
    for x in eigvecs:
        # Compute projection: proj = (x·v / x·x) * x
        # Use vdot for proper complex conjugation: vdot(a,b) = sum(conj(a) * b)
        proj = (np.vdot(x, v) / np.vdot(x, x)) * x
        v = v - proj

    # Normalize
    nv = np.linalg.norm(v)
    if nv > 1e-14:
        v = v / nv
    else:
        # If deflation results in zero vector, return random orthogonal vector
        print("Warning: Deflation resulted in near-zero vector, using random initialization")
        v = np.random.randn(len(v)) + 1j * np.random.randn(len(v))
        v = project_out(v, eigvecs)  # Recursive call to ensure orthogonality
    return v


def make_power_iteration(V_r, f_r, n_eig=5, max_iter=100, tol=1e-6):
    """
    Solve BCS eigenvalue problem using power iteration with projection.

    This method finds multiple eigenpairs by:
    1. Power iteration to find largest eigenvalue
    2. Project out found eigenvector
    3. Repeat for next eigenvalue

    Args:
        V_r: Vertex in real space, shape (nkx, nky, nkz, nstates, nstates, nstates, nstates)
        f_r: Form factor in real space
        n_eig: Number of eigenpairs to find
        max_iter: Maximum iterations per eigenpair
        tol: Convergence tolerance

    Returns:
        eigenvalues (array of shape (n_eig,)), eigenvectors (array of shape (n, n_eig))
    """
    n = nx * ny * nz * nstates * nstates

    # Define matrix-vector product
    def mv(D_k_flat):
        """Apply the kernel: result = V * D"""
        result = multiply_basic(V_r, f_r, D_k_flat, mesh=(nx, ny, nz), arrsize=(nstates, nstates))
        return result

    eigenvalues = []
    eigenvectors = []
    old_eigvecs = []  # Store previously found eigenvectors for projection

    print(f"\n{'='*70}")
    print(f"Power Iteration Solver (finding {n_eig} eigenpairs)")
    print(f"{'='*70}")

    for ieig in range(n_eig):
        print(f"\nSearching for eigenpair #{ieig+1}...")

        # Random initial guess
        v = np.random.randn(n) + 1j * np.random.randn(n)
        v = v / np.linalg.norm(v)

        # Project out previously found eigenvectors
        if len(old_eigvecs) > 0:
            v = project_out(v, old_eigvecs)

        eig = 0.0
        prev_eig = 0.0

        for it in range(max_iter):
            # Power iteration step: v_new = A * v
            v_new = mv(v)

            # Project out old eigenvectors to avoid convergence to them
            if len(old_eigvecs) > 0:
                v_new = project_out(v_new, old_eigvecs)

            # Compute eigenvalue (Rayleigh quotient): λ = v† A v / v† v
            # Since v is normalized, this simplifies to: λ = v† v_new
            eig = np.vdot(v, v_new).real

            # Normalize
            norm = np.linalg.norm(v_new)
            v_new = v_new / norm

            # Check convergence
            diff = np.abs(eig - prev_eig)
            prev_eig = eig
            v = v_new

            if (it + 1) % 10 == 0:
                print(f"  Iteration {it+1:3d}: λ = {eig:12.6f}, Δλ = {diff:.6e}")

            if diff < tol or np.isnan(diff):
                print(f"  Converged after {it+1} iterations: λ = {eig:.6f}")
                break

        if it == max_iter - 1:
            print(f"  Warning: Did not converge after {max_iter} iterations (Δλ = {diff:.6e})")

        eigenvalues.append(eig)
        eigenvectors.append(v.copy())
        old_eigvecs.append(v.copy())

        # Stop if eigenvalue is too small (hit null space)
        if np.abs(eig) < 1e-3:
            print(f"  Eigenvalue too small (|λ| < 1e-3), stopping search")
            break

    # Convert to arrays and sort by magnitude
    eigenvalues = np.array(eigenvalues)
    eigenvectors = np.column_stack(eigenvectors)

    # Sort by descending |λ|
    idx = np.argsort(np.abs(eigenvalues))[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    print(f"\n{'='*70}")
    print(f"Power Iteration Complete")
    print(f"{'='*70}")
    print(f"Found {len(eigenvalues)} eigenpairs:")
    for i, eig in enumerate(eigenvalues):
        print(f"  λ_{i+1} = {eig:12.6f}")

    return eigenvalues, eigenvectors


def make_lanczos(V_r, f_r, f_k):
    """
    Set up and solve eigenvalue problem on k-grid using Lanczos (ARPACK)

    Args:
        V_r: Vertex in real space, shape (nkx, nky, nkz, nstates, nstates, nstates, nstates)

    Returns:
        eigenvalues, eigenvectors
    """
    n = nx * ny * nz * nstates * nstates

    # Define matrix-vector product for the eigenvalue problem
    def mv(D_k_flat):
        """
        Apply the kernel: result = V * D
        where the multiplication is a convolution in real space
        """
        result = multiply(V_r, f_k, f_r, D_k_flat, mesh=(nx, ny, nz), arrsize=(nstates, nstates))
        return result

    # Create linear operator for ARPACK
    A = LinearOperator((n, n), matvec=mv, dtype=complex)

    # Check both ends of spectrum to avoid null space
    k_check = min(10, n // 2 - 1)

    # Find most positive eigenvalues
    try:
        vals_pos, vecs_pos = eigsh(A, k=k_check, which='LA', tol=1e-10)
    except:
        vals_pos = np.array([])
        vecs_pos = np.zeros((n, 0))

    # Find most negative eigenvalues
    try:
        vals_neg, vecs_neg = eigsh(A, k=k_check, which='SA', tol=1e-10)
    except:
        vals_neg = np.array([])
        vecs_neg = np.zeros((n, 0))

    # Combine and sort by magnitude (largest |λ| first)
    vals = np.concatenate([vals_pos, vals_neg])
    vecs = np.column_stack([vecs_pos, vecs_neg])

    idx = np.argsort(np.abs(vals))[::-1]
    vals = vals[idx]
    vecs = vecs[:, idx]

    # Return top 5 eigenpairs
    n_return = min(5, len(vals))
    return vals[:n_return], vecs[:, :n_return]

def bcs():
    H_r, kmesh, e_k = fly.load_triqs_H.get_energy_mesh()
    kpts = get_k_mesh(BZ, nx, ny, nz)

    vertex_file = outdir + prefix + '_vertex.h5'
    print(f"Loading vertex from {vertex_file}")
    V_vq = fly.Field_CM(vertex_file)
    V_q = np.array(V_vq(kpts))
    #V_q[:, 0, 0] = 1.0
    #V_q[:, 0, 0] = np.cos(kpts[:,0]) + np.cos(kpts[:,1])
    #V_q[:, 0, 0] = np.exp(-kpts[:,0]**2 - kpts[:,1]**2)
    nk = V_q.shape[0]
    V_q = V_q.reshape(nk, nstates, nstates, nstates, nstates)
    V_q = V_q.reshape(nx, ny, nz, nstates, nstates, nstates, nstates)

    # Transform V from k-space to real-space (spatial dimensions only)
    V_r = fftn(V_q, axes=(0, 1, 2))
    print(f"V_r shape: {V_r.shape}")

    sigma_file = outdir + prefix + '_sigma_iw.h5'
    Sigma_w = fly.Field_C(sigma_file)
    dw = 1e-2
    Z = 1.0 - (Sigma_w(dw).imag - Sigma_w(-dw).imag) / (dw)
    print(f"Quasiparticle weight Z: {Z}")

    norb1, norb2 = e_k.data.shape[-2], e_k.data.shape[-1]
    e_k = np.reshape(e_k.data - mu, (nx, ny, nz, norb1, norb2))
    f_k = np.tanh(beta * e_k / Z) / (2.0 * e_k)
    f_r = fftn(f_k.reshape(nx, ny, nz, 1, 1), axes=(0, 1, 2))
    #f_r[:] = 1.0
    f_k_sqrt = np.sqrt(f_k)
    f_r_sqrt = fftn(f_k_sqrt.reshape(nx, ny, nz, 1, 1), axes=(0, 1, 2))

    # Solve BCS eigenvalue problem
    method = cfg.method if hasattr(cfg, 'method') else ''
    print(f"Solving BCS gap equation using {method}...")

    eigs, vecs = make_power_iteration(V_r, f_r, n_eig=5, max_iter=100, tol=1e-6)
        #if method == 'power_iteration':
        #    eigs, vecs = make_power_iteration(V_r, f_r, f_sqrt, n_eig=5, max_iter=100, tol=1e-6)
        #elif method == "lanczos":
        #    eigs, vecs = make_lanczos(V_r, f_r, f_sqrt)
        #else:
        #    raise ValueError(f"Unknown method '{method}' for BCS solver")

    eigs /= (nk * Z)

    for i in range(len(eigs)):
        if abs(eigs[i]) > 1e-3:
            plot_gap(vecs[:,i].reshape(nx, ny))

    print(f"Eigenvalues: {eigs}")
    print(f"Largest eigenvalue: {eigs[0]:.6f}")
    gap_file = outdir + prefix + '_gap.h5'
    gap = vecs[:,0].reshape(nx, ny, nz)
    fly.save_data(gap_file, gap, mesh=[nx, ny, nz], domain=BZ[:2,:2])
    print(f"Saved gap function to {gap_file}")

def plot_gap(gap):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    im0 = ax[0].contourf(gap.real, levels=50, cmap='plasma')
    fig.colorbar(im0, ax=ax[0], label='real Δ(k)')
    im1 = ax[1].contourf(gap.imag, levels=50, cmap='viridis')
    fig.colorbar(im1, ax=ax[1], label='imag Δ(k)')
    plt.show()
