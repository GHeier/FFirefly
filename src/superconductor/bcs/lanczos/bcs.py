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

# Load config variables on file call
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
max_eigs_searched = cfg.num_solutions

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

def f(e, beta, eps=1e-6):
    """
    Elementwise: tanh(beta*e)/(2e), with the e->0 limit = beta/2.
    NOTE: your code used beta/4; the actual limit of tanh(beta e)/(2e) is beta/2.
    If you truly want beta/4, change the constant below.
    """
    e = np.asarray(e)
    small = np.abs(e) < eps
    out = np.empty_like(e, dtype=np.result_type(e, 1.0))
    out[small] = beta / 4.0
    out[~small] = np.tanh(beta * e[~small] / 2) / (2.0 * e[~small])
    return out

def construct_form_factor(e_k, Z):
    """Construct form factor f(k) = tanh(βε_k/Z) / (2ε_k)"""
    norb1, norb2 = e_k.data.shape[-2], e_k.data.shape[-1]
    e_k = np.reshape(e_k.data - mu, (nx, ny, nz, norb1, norb2))
    f_k = f(e_k / Z, beta) / Z
    print("Minum abs value of ek:", np.min(np.abs(e_k)))
    # Keep f_k in momentum space (reshape to match Delta dimensions)
    f_k = f_k.reshape(nx, ny, nz, nstates, nstates)
    return f_k

def load():
    """Load data needed for BCS calculation (FFT-based convolution method)"""
    H_r, kmesh, e_k = fly.load_triqs_H.get_energy_mesh()
    kpts = get_k_mesh(BZ, nx, ny, nz)

    # Load vertex
    vertex_file = outdir + prefix + '_vertex.h5'
    print(f"Loading vertex from {vertex_file}")
    V_vq = fly.Field_CM(vertex_file)
    V_q = np.array(V_vq(kpts))

    nk = V_q.shape[0]
    V_q = V_q.reshape(nk, nstates, nstates, nstates, nstates)
    V_q = V_q.reshape(nx, ny, nz, nstates, nstates, nstates, nstates)

    # Apply k=0 factors
    V_q[0,:, :, :, :, :, :] *= 0.5
    V_q[:,0, :, :, :, :, :] *= 0.5
    if dim > 2:
        V_q[:,:, 0, :, :, :, :] *= 0.5

    # Transform V from k-space to real-space (spatial dimensions only)
    V_r = fftn(V_q, axes=(0, 1, 2))
    print(f"V_r shape: {V_r.shape}")

    # Load self-energy and compute quasiparticle weight
    sigma_file = outdir + prefix + '_sigma_iw.h5'
    Sigma_w = fly.Field_C(sigma_file)
    dw = 1e-2
    Z = 1.0 - (Sigma_w(dw).imag - Sigma_w(-dw).imag) / (dw)
    print(f"Quasiparticle weight Z: {Z}")

    # Compute form factor f(k) = tanh(βε_k/Z) / (2ε_k)
    f_k = construct_form_factor(e_k, Z)

    # Initialize Delta0 (starting gap function)
    Delta0_shape = (nx, ny, nz, nstates, nstates)
    Delta0 = np.zeros(Delta0_shape, dtype=complex)

    return V_r, f_k, Z, Delta0


def load_matrix():
    """Load data and build full V(k-k') matrix for direct matrix multiplication"""
    import time
    start_time = time.time()

    H_r, kmesh, e_k = fly.load_triqs_H.get_energy_mesh()
    kpts = get_k_mesh(BZ, nx, ny, nz)
    nk = len(kpts)

    print(f"Building full V(k-k') matrix for {nk} k-points...")
    print(f"WARNING: This requires O(nk^2) = {nk**2/1e6:.1f}M evaluations and may take several minutes")

    # Load vertex field
    vertex_file = outdir + prefix + '_vertex.h5'
    print(f"Loading vertex from {vertex_file}")
    V_field = fly.Field_CM(vertex_file)

    # Compute all momentum differences k - k'
    # For now, assume single band (nstates=1) for simplicity
    # V_matrix will be shape (nk, nk) for single band
    # or (nk*nstates, nk*nstates) for multi-band

    if nstates == 1:
        # Single band case: V_matrix[k, k'] = V(k-k')[0,0,0,0]
        V_matrix = np.zeros((nk, nk), dtype=complex)
        print("Computing V(k-k') matrix (single band)...")

        dk = kpts[:, None, :] - kpts[None, :, :]
        V_vals = V_field(dk.reshape(-1, 3))
        V_vals = np.array(V_vals).reshape(nk, nk)
        V_matrix[:, :] = V_vals

    else:
        raise NotImplementedError("Multi-band V(k-k') matrix construction not implemented yet")

    # Load self-energy and compute quasiparticle weight
    sigma_file = outdir + prefix + '_sigma_iw.h5'
    Sigma_w = fly.Field_C(sigma_file)
    dw = 1e-2
    Z = 1.0 - (Sigma_w(dw).imag - Sigma_w(-dw).imag) / (dw)
    print(f"Quasiparticle weight Z: {Z}")

    # Compute form factor f(k) = tanh(βε_k/Z) / (2ε_k)
    norb1, norb2 = e_k.data.shape[-2], e_k.data.shape[-1]
    e_k_flat = e_k.data.reshape(nk, norb1, norb2) - mu
    f_k = np.tanh(beta * e_k_flat / Z) / (2.0 * e_k_flat)  # shape (nk, nstates, nstates)

    # Initialize Delta0
    if nstates == 1:
        Delta0 = np.zeros(nk, dtype=complex)
    else:
        Delta0 = np.zeros((nk, nstates, nstates), dtype=complex)

    print(f"V_matrix shape: {V_matrix.shape}")
    print(f"f_k shape: {f_k.shape}")
    print(f"Delta0 shape: {Delta0.shape}")

    return V_matrix, f_k, Z, Delta0


# Main function 1
def run_lanczos():
    V_r, f_r, Z, Delta0 = load()
    eigs, Deltas = solve_bcs_lanczos(V_r, f_r, Z, Delta0)

    # Find max positive eigenvalue
    i = np.where(eigs > 0, eigs, -np.inf).argmax()
    print(f"Max Eig: {eigs[i]:.6f}")

    # Save gap function
    gap_file = outdir + prefix + '_gap.h5'
    gap = Deltas[:, i].reshape(nx, ny, nz)
    fly.save_data(gap_file, gap, mesh=[nx, ny, nz], domain=BZ[:2,:2])
    print(f"Saved gap function to {gap_file}")
    plot_gap(gap.reshape(nx, ny))


# Main function 2
def run_power_iteration():
    V_r, f_r, Z, Delta0 = load()
    eig, Delta = solve_bcs_power_iteration(V_r, f_r, Z, Delta0)
    print(f"Max Eig: {eig:.6f}")

    # Save gap function
    gap_file = outdir + prefix + '_gap.h5'
    gap = Delta.reshape(nx, ny, nz)
    fly.save_data(gap_file, gap, mesh=[nx, ny, nz], domain=BZ[:2,:2])
    print(f"Saved gap function to {gap_file}")
    plot_gap(gap.reshape(nx, ny))


# Main function 3 - Full matrix method
def run_matrix_power_iteration():
    V_matrix, f_k, Z, Delta0 = load_matrix()
    eig, Delta = solve_bcs_matrix_power_iteration(V_matrix, f_k, Z, Delta0)
    print(f"Max Eig: {eig:.6f}")

    # Save gap function
    gap_file = outdir + prefix + '_gap.h5'
    if nstates == 1:
        gap = Delta.reshape(nx, ny, nz)
    else:
        # Multi-band case: save all bands
        gap = Delta.reshape(nx, ny, nz, nstates, nstates)
    fly.save_data(gap_file, gap, mesh=[nx, ny, nz], domain=BZ[:2,:2])
    print(f"Saved gap function to {gap_file}")
    plot_gap(gap.reshape(nx, ny))

def solve_bcs_lanczos(V_r, f_k, Z, Delta0):
    """Returns multiple eigenvalues/vectors for case of competing solutions"""

    # Get shape for flattening
    shape = Delta0.shape
    n = np.prod(shape)
    nk = nx * ny * nz

    # Define matrix-vector product for BCS kernel
    def mv(D_flat):
        # Load in unflattened data, perform the convolution, and flatten result
        D_k = D_flat.reshape(shape)
        Delta_new = BCS_step(V_r, f_k, D_k)
        return Delta_new.flatten()

    # Create linear operator for ARPACK
    A = LinearOperator((n, n), matvec=mv, dtype=complex)

    k_check = max_eigs_searched
    print(f"Searching for {k_check} eigenvalues at each end of spectrum...")

    # Find most positive eigenvalues
    eigs, vecs = eigsh(A, k=k_check, which='LM', tol=1e-8, maxiter=1000)

    # Scale eigenvalues by quasiparticle weight and k-mesh size
    eigs /= (nk * Z)

    for i, eig in enumerate(eigs):
        print(f"     eig{i}: {eig:12.6f}")

    # Return top eigenpairs (vecs is already in column format)
    return eigs, vecs


def solve_bcs_power_iteration(V_r, f_k, Z, Delta0):
    max_iter = 100
    tol = 1e-4
    Delta = Delta0.copy()

    shape = Delta.shape
    nk = nx * ny * nz

    eig = 0.0
    prev_eig = 0.0
    diff = 1.0
    old_Deltas = []

    while eig <= 0.0 and len(old_Deltas) < max_eigs_searched:
        print("shape: ", shape)
        Delta[:] = np.random.rand(*shape) + 1j * np.random.rand(*shape)
        iter = 0

        for it in range(max_iter):
            Delta_new = BCS_step(V_r, f_k, Delta)

            # Compute eigenvalue (Rayleigh quotient): eig = <D|K*D> / <D|D> = <D|D_new> / <D|D>
            norm = np.sum(Delta * np.conj(Delta)).real
            eig = np.sum(np.conj(Delta) * Delta_new).real / norm

            # Project out previously found eigenvectors
            Delta_new_flat = Delta_new.flatten()
            Delta_new_flat = project_out(Delta_new_flat, old_Deltas)
            Delta_new = Delta_new_flat.reshape(shape)

            # Check convergence
            diff = np.abs(eig - prev_eig)
            prev_eig = eig

            # Normalize new Delta
            norm = np.sum(Delta_new * np.conj(Delta_new)).real
            Delta_new = Delta_new / np.sqrt(norm)

            Delta = Delta_new.copy()
            print(f"Eig: {eig} Error = {diff:.6e}")
            iter = it
            if (diff < tol and it > 10) or np.isnan(diff):
                break

        print(f"Iterations: {iter+1}")
        print(f"eig{len(old_Deltas)} = {eig}")
        old_Deltas.append(Delta.flatten().copy())

    # Scale by quasiparticle weight and k-mesh size
    eig /= (nk * Z)

    return eig, Delta.flatten()


def solve_bcs_matrix_power_iteration(V_matrix, f_k, Z, Delta0):
    """
    Solve BCS eigenvalue problem using full matrix V(k-k') with power iteration

    The BCS equation in matrix form:
        Delta(k) = sum_{k'} V(k-k') * f(k') * Delta(k')

    Args:
        V_matrix: Full interaction matrix V(k-k')
                  Single band: shape (nk, nk)
                  Multi-band: shape (nk, nstates, nstates, nk, nstates, nstates)
        f_k: Form factor at each k-point, shape (nk, nstates, nstates)
        Z: Quasiparticle weight
        Delta0: Initial gap function

    Returns:
        eig: Largest eigenvalue
        Delta: Corresponding eigenvector (gap function)
    """
    max_iter = 100
    tol = 1e-4
    Delta = Delta0.copy()

    nk = len(f_k)
    shape = Delta.shape

    eig = 0.0
    prev_eig = 0.0
    diff = 1.0
    old_Deltas = []

    while eig <= 0.0 and len(old_Deltas) < max_eigs_searched:
        print("shape: ", shape)
        Delta[:] = np.random.rand(*shape) + 1j * np.random.rand(*shape)
        iter = 0

        for it in range(max_iter):
            Delta_new = BCS_step_matrix(V_matrix, f_k, Delta)

            # Compute eigenvalue (Rayleigh quotient): eig = <D|K*D> / <D|D> = <D|D_new> / <D|D>
            norm = np.sum(Delta * np.conj(Delta)).real
            eig = np.sum(np.conj(Delta) * Delta_new).real / norm

            # Project out previously found eigenvectors
            Delta_new_flat = Delta_new.flatten()
            Delta_new_flat = project_out(Delta_new_flat, old_Deltas)
            Delta_new = Delta_new_flat.reshape(shape)

            # Check convergence
            diff = np.abs(eig - prev_eig)
            prev_eig = eig

            # Normalize new Delta
            norm = np.sum(Delta_new * np.conj(Delta_new)).real
            Delta_new = Delta_new / np.sqrt(norm)

            Delta = Delta_new.copy()
            print(f"Eig: {eig} Error = {diff:.6e}")
            iter = it
            if (diff < tol and it > 10) or np.isnan(diff):
                break

        print(f"Iterations: {iter+1}")
        print(f"eig{len(old_Deltas)} = {eig}")
        old_Deltas.append(Delta.flatten().copy())

    # Scale by quasiparticle weight and k-mesh size
    eig /= (nk * Z)

    return eig, Delta.flatten()


def BCS_step(V_r, f_k, Delta):
    """
    Perform one step of BCS iteration: Delta_new = V * (f * Delta)
    (FFT-based convolution method)

    IMPORTANT: The correct BCS equation is:
        Delta_new(k) = sum_k' V(k-k') * f(k') * Delta(k')

    In FFT form, f(k) must be applied BEFORE the transform:
        Delta_new = IFFT[V(r) * FFT[f(k) * Delta(k)]]

    Args:
        V_r: Vertex in real space, shape (nx, ny, nz, nstates, nstates, nstates, nstates)
        f_k: Form factor in k-space, shape (nx, ny, nz, 1, 1) or (nx, ny, nz, nstates, nstates)
        Delta: Gap function in k-space, shape (nx, ny, nz, nstates, nstates)

    Returns:
        Delta_new: Updated gap function in k-space
    """
    Delta_reverse = Delta[::-1, ::-1, ::-1, :, :]
    # CRITICAL: Multiply f(k) * Delta(k) in k-space FIRST
    f_Delta_k = f_k * (Delta + Delta_reverse) / 2 # symmetrize Delta

    # Transform to real space
    f_Delta_r = fftn(f_Delta_k, axes=(0, 1, 2))

    # Convolve with vertex in real space
    h_r = -np.einsum('xyzabcd,xyzcd->xyzab', V_r, f_Delta_r)

    # Transform back to k-space
    h_k = ifftn(h_r, axes=(0, 1, 2))

    return h_k


def BCS_step_matrix(V_matrix, f_k, Delta):
    """
    Perform one step of BCS iteration using full matrix multiplication:
        Delta_new(k) = sum_{k'} V(k-k') * f(k') * Delta(k')

    Args:
        V_matrix: Full interaction matrix V(k-k')
                  Single band: shape (nk, nk)
                  Multi-band: shape (nk, nstates, nstates, nk, nstates, nstates)
        f_k: Form factor at each k-point, shape (nk, nstates, nstates)
        Delta: Gap function, shape (nk,) for single band or (nk, nstates, nstates) for multi-band

    Returns:
        Delta_new: Updated gap function, same shape as Delta
    """
    if len(Delta.shape) == 1:
        # Single band case: Delta[k], f_k[k,0,0], V_matrix[k,k']
        # Delta_new[k] = sum_{k'} V[k,k'] * f[k',0,0] * Delta[k']
        f_Delta = f_k[:, 0, 0] * Delta  # shape (nk,)
        Delta_new = V_matrix @ f_Delta  # shape (nk,)
    else:
        # Multi-band case: Delta[k,a,b], f_k[k,c,d], V_matrix[k,a,b,k',c,d]
        # Delta_new[k,a,b] = sum_{k',c,d} V[k,a,b,k',c,d] * f[k',c,d] * Delta[k',c,d]
        f_Delta = f_k * Delta  # shape (nk, nstates, nstates)
        # Reshape for einsum: V[k,a,b,k',c,d] * f_Delta[k',c,d] -> Delta_new[k,a,b]
        Delta_new = np.einsum('kabkcd,kcd->kab', V_matrix, f_Delta)

    return Delta_new


def project_out(v, eigvecs):
    """Project out previously found eigenvectors using Gram-Schmidt orthogonalization."""
    for x in eigvecs:
        # Compute projection: proj = (x*v / x*x) * x
        # Use vdot for proper complex conjugation: vdot(a,b) = sum(conj(a) * b)
        proj = (np.vdot(x, v) / np.vdot(x, x)) * x
        v = v - proj

    # Normalize
    nv = np.linalg.norm(v)
    return v / nv


def plot_gap(gap):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    im0 = ax.contourf(gap.real, levels=50, cmap='bwr')
    fig.colorbar(im0, ax=ax, label='real Δ(k)')
    ax.set_axis_off()

    #im1 = ax[1].contourf(gap.imag, levels=50, cmap='viridis')
    #fig.colorbar(im1, ax=ax[1], label='imag Δ(k)')
    #ax[1].set_axis_off()

    plt.show()
