import firefly as fly
import firefly.config as cfg
from firefly.diagram import *

from triqs.gf import *
from triqs.gf.mesh_product import MeshProduct
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

import load_triqs_H

# Load config variables on file call
outdir = cfg.outdir
prefix = cfg.prefix

nbnd = cfg.nbnd
Nk = cfg.k_mesh[0]
BZ = get_brillouin_zone()

mu = cfg.fermi_energy
beta = 1.0 / cfg.Temperature
max_eigs_searched = cfg.num_eigenvalues_to_save

def load():
    H_r, kmesh, e_k = load_triqs_H.get_energy_mesh()
    emax = e_k.data.max().real
    emin = e_k.data.min().real
    print(f"emax: {emax:.4f}, emin: {emin:.4f}")
    # Build Discrete Lehman Representation (DLR) mesh for imaginary frequencies
    DLRImMesh = load_triqs_H.create_dlr_meshes(e_k, beta, statistic='Fermion')
    k_mesh = MeshBrZone(BZ, n_k=Nk)   # uniform Nk x Nk x Nk (third dim is 1 if 2D)

    # Create TRIQS object for Green's function G(k, iw)
    wk_mesh = MeshProduct(DLRImMesh, k_mesh)
    G = Gf(mesh=wk_mesh, target_shape=[nbnd, nbnd])
    G = fly.diagram.Diagram(G, 'Fermion')

    # Load Green's function using Firefly interface
    G_file = outdir + prefix + '_G.h5'
    print(f"Loading G from {G_file}")
    G_data = fly.Field_CM(G_file)
    G.load(G_data)

    DLRImMesh = load_triqs_H.create_dlr_meshes(e_k, beta, statistic='Boson')
    wk_mesh = MeshProduct(DLRImMesh, k_mesh)
    V = Gf(mesh=wk_mesh, target_shape=[nbnd, nbnd, nbnd, nbnd])
    V = fly.diagram.Diagram(V, 'Boson')

    vertex = fly.Field_CM(outdir + prefix + '_vertex.h5')
    V.load(vertex)

    V.wk_to_tr()

    Delta0 = G.copy()

    return G, V, Delta0


# Main function 1
def run_lanczos():
    G, V, Delta0 = load()
    eigs, Deltas = solve_eliashberg_lanczos(G, V, Delta0)
    # Find max positive eigenvalue
    i = np.where(eigs > 0, eigs, -np.inf).argmax() 
    print(f"Max Eig: {eigs[i]:.6f}")
    Deltas[i].save(outdir + prefix + '_gap.h5')
    for i in range(len(eigs)):
        print(f"Saving eig{i}: {eigs[i]:.6f}")
        Deltas[i].save(outdir + prefix + f'_gap_eig{i}.h5')
    return eigs[i]

# Main function 2
def run_power_iteration():
    G, V, Delta0 = load()
    eig, Delta = solve_eliashberg_power_iteration(G, V, Delta0)
    print(f"Max Eig: {eig:.6f}")
    Delta.save(outdir + prefix + '_gap.h5')


def solve_eliashberg_lanczos(G, V, Delta0):
    """ Returns multiple eigenvalues/vectors for case of competing solutions """

    # Properly flip G(k, iw) to G(-k, -iw) for Cooper pair formation
    G_flip = flip_wk(G)

    # Get shape for flattening
    shape = Delta0.obj_wk.data.shape
    n = np.prod(shape)

    # Define matrix-vector product for Eliashberg kernel
    def mv(D_flat):
        # Load in unflattened data, perform the convolution, and flatten result
        Delta0.obj_wk.data[:] = D_flat.reshape(shape)
        Delta_new = Eliashberg_step(G, G_flip, V, Delta0)
        return Delta_new.obj_wk.data.flatten()

    # Create linear operator for ARPACK
    A = LinearOperator((n, n), matvec=mv, dtype=complex)

    k_check = max_eigs_searched
    print(f"Searching for {k_check} eigenvalues...")

    # Find most positive eigenvalues
    eigs, vecs = eigsh(A, k=k_check, which='LA', tol=1e-8, maxiter=1000)
    for i, eig in enumerate(eigs):
          print(f"     eig{i}: {eig:12.6f}")

    # Convert eigenvectors back to Diagram objects
    vecs_list = []
    for i in range(k_check):
        Delta_eig = Delta0.copy()
        Delta_eig.obj_wk.data[:] = vecs[:, i].reshape(shape)
        vecs_list.append(Delta_eig)

    # Return top 5 eigenpairs
    return eigs, vecs_list

def solve_eliashberg_power_iteration(G, V, Delta0):
    max_iter = 100
    tol = 1e-4
    Delta = Delta0.copy()

    # Properly flip G(k, iw) to G(-k, -iw) for Cooper pair formation
    G_flip = flip_wk(G)

    mesh, BZ = extract_mesh_and_bz(Delta.obj_wk)

    eig = 0.0
    prev_eig = 0.0
    diff = 1.0
    old_Deltas = []

    while eig <= 0.0 and len(old_Deltas) < max_eigs_searched:
        shape = Delta.obj_wk.data.shape
        print("shape: ", shape)
        Delta.obj_wk.data[:] = np.random.rand(*Delta.obj_wk.data.shape) + 1j * np.random.rand(*Delta.obj_wk.data.shape)
        iter = 0
        for it in range(max_iter):
            Delta_new = Eliashberg_step(G, G_flip, V, Delta)

            D_wk = Delta.obj_wk.data
            new_D_wk = Delta_new.obj_wk.data

            # Compute eigenvalue (Rayleigh quotient): eig = <D|K*D> / <D|D> = <D|D_new> / <D|D>
            norm = np.sum(D_wk * np.conj(D_wk)).real
            eig = np.sum(np.conj(D_wk) * new_D_wk).real / norm

            # Project out previously found eigenvectors
            new_D_wk = project_out(new_D_wk, old_Deltas)

            # Check convergence
            diff = np.abs(eig - prev_eig)
            prev_eig = eig

            # Normalize new Delta
            norm = np.sum(new_D_wk * np.conj(new_D_wk)).real
            Delta_new.obj_wk.data[:] = new_D_wk / np.sqrt(norm)

            Delta = Delta_new.copy()
            print(f"Eig: {eig} Error = {diff:.6e}")
            iter = it
            if (diff < tol and it > 10) or np.isnan(diff):
                break

        print(f"Iterations: {iter+1}")
        print(f"eig{len(old_Deltas)} = {eig}")
        old_Deltas.append(Delta.obj_wk.data.copy())

    return eig, Delta

def project_even_odd(Delta, parity='even'):
    flipped = flip_k(Delta)
    if parity == 'even':
        Delta.obj_wk.data[:] = 0.5 * (Delta.obj_wk.data + flipped.obj_wk.data)
    elif parity == 'odd':
        Delta.obj_wk.data[:] = 0.5 * (Delta.obj_wk.data - flipped.obj_wk.data)
    else:
        raise ValueError("parity must be 'even' or 'odd'")
    return Delta

def Eliashberg_step(G, G_flip, V, Delta):
    F = Delta.copy()
    #Delta = project_even_odd(Delta, parity='even')
    # F = -G(k,iw) * G(-k,-iw) * Delta(k,iw)
    F.obj_wk.data[:] = -1.0 * G_flip.obj_wk.data * np.conj(G.obj_wk.data) * Delta.obj_wk.data
    F.wk_to_tr()
    Delta_new = dot_tr(V, F)
    Delta_new.tr_to_wk()
    return Delta_new

def project_out(v, eigvecs):
    """Project out previously found eigenvectors using Gram-Schmidt orthogonalization."""
    for x in eigvecs:
        # Compute projection: proj = (x*v / x*x) * x
        # Use vdot for proper complex conjugation: vdot(a,b) = sum(conj(a) * b)
        v_flat = v.flatten()
        x_flat = x.flatten()
        proj = (np.vdot(x_flat, v_flat) / np.vdot(x_flat, x_flat)) * x
        v = v - proj

    # Normalize
    nv = np.linalg.norm(v)
    return v / nv
