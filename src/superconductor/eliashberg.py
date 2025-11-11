import firefly as fly
import firefly.config as cfg
from firefly.diagram import *
from firefly.diagram import flip_wk

from triqs.gf import *
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs.lattice import BrillouinZone, BravaisLattice
from triqs.gf.mesh_product import MeshProduct
import numpy as np
from triqs_tprf.lattice import *
from triqs_tprf import *
#from interface_triqs import *
#from load_triqs_H import *

outdir = cfg.outdir
prefix = cfg.prefix

nstates = cfg.nstates
Nk = cfg.k_mesh[0]
BZ = get_brillouin_zone()

mu = cfg.fermi_energy
beta = 1.0 / cfg.Temperature

def analyze_gap_symmetry(Delta):
    """Analyze the symmetry of the gap function."""
    mesh, BZ = extract_mesh_and_bz(Delta.obj_wk)
    original_shape = Delta.obj_wk.data.shape
    nk = int(np.sqrt(original_shape[1]))

    # Reshape to (nw, nkx, nky, nkz, orb, orb) and take lowest Matsubara frequency
    data = np.reshape(Delta.obj_wk.data, mesh + original_shape[2:])
    gap_k = data[0, :, :, 0, 0, 0]  # lowest frequency, 2D slice

    print(f"\nGap function analysis:")
    print(f"  Shape: {gap_k.shape}")
    print(f"  Sign changes along kx: {np.sum(np.diff(np.sign(gap_k[:, nk//2].real)) != 0)}")
    print(f"  Sign changes along ky: {np.sum(np.diff(np.sign(gap_k[nk//2, :].real)) != 0)}")
    print(f"  Corner values: Δ(0,0)={gap_k[0,0]:.4f}, Δ(π,0)={gap_k[nk//2,0]:.4f}, Δ(π,π)={gap_k[nk//2,nk//2]:.4f}")
    print(f"  d-wave test: Δ(π,0)/Δ(0,π) = {gap_k[nk//2,0]/gap_k[0,nk//2]:.4f} (expect ~1 for d-wave)")
    print(f"  s-wave test: Δ(0,0)/Δ(π,π) = {gap_k[0,0]/gap_k[nk//2,nk//2]:.4f} (expect >0 for s-wave)")


def print_vertex_structure(V, label):
    """Print diagnostic info about vertex structure."""
    mesh, BZ = extract_mesh_and_bz(V.obj_wk)
    original_shape = V.obj_wk.data.shape
    nw, nk_total = original_shape[0], original_shape[1]
    nk = int(np.sqrt(nk_total))

    # Reshape and average over frequency for visualization
    data = np.reshape(V.obj_wk.data, mesh + original_shape[2:])
    vertex_k = np.mean(np.abs(data[:, :, :, 0, 0, 0, 0, 0]), axis=0)  # average over freq

    max_idx = np.unravel_index(np.argmax(vertex_k), vertex_k.shape)
    print(f"\n{label}:")
    print(f"  Max at k=({max_idx[0]}, {max_idx[1]}) (expect (π,π) at ({nk//2}, {nk//2}))")
    print(f"  Max value: {np.max(vertex_k):.4f}")
    print(f"  Corner values (π,π): V[0,0]={vertex_k[0,0]:.4f}, V[{nk//2},{nk//2}]={vertex_k[nk//2,nk//2]:.4f}")


def symmetrize_vertex(V):
    """
    Symmetrize the vertex by computing V_sym(q) = (V(q) + V(-q)) / 2.

    Args:
        V: Diagram object containing the vertex in wk space
    """
    # Extract mesh dimensions
    mesh, BZ = extract_mesh_and_bz(V.obj_wk)
    original_shape = V.obj_wk.data.shape

    # Reshape to mesh dimensions (nw, nkx, nky, nkz, orb1, orb2, orb3, orb4)
    data = np.reshape(V.obj_wk.data, mesh + original_shape[2:])

    # fftshift k-axes to center k-points
    k_axes = tuple(range(1, len(mesh)))  # axes 1, 2, 3 for kx, ky, kz
    data = np.fft.fftshift(data, axes=k_axes)

    # Flip k-axes to get V(-q) and average
    data = (data + np.flip(data, axis=k_axes)) / 2.0

    # ifftshift back to original k-ordering
    data = np.fft.ifftshift(data, axes=k_axes)

    # Reshape back to TRIQS data shape
    V.obj_wk.data[:] = np.reshape(data, original_shape)

    print("Vertex symmetrized: V_sym(q) = (V(q) + V(-q)) / 2")

def main():
    H_r, kmesh, e_k = fly.load_triqs_H.get_energy_mesh()
    emax = e_k.data.max().real
    emin = e_k.data.min().real
    print(f"emax: {emax}, emin: {emin}")
    DLRImMesh = fly.load_triqs_H.create_dlr_meshes(e_k, beta, statistic='Fermion')
    print(BZ)
    k_mesh = MeshBrZone(BZ, n_k=Nk)   # uniform Nk x Nk x Nk (third dim is 1 if 2D)

    #sigma = fly.Field_C(outdir + prefix + '_sigma.h5')
    wk_mesh = MeshProduct(DLRImMesh, k_mesh) 
    #E = Gf(mesh=wk_mesh, target_shape=[1,1])
    #E.data[:, :, 0, 0] = sigma.get_data()
    #fly.interface_triqs.fill_triqs_from_field(E, sigma)
    #G0 = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)
    #G = inverse(inverse(G0) - E)
    G_data = fly.Field_CM(outdir + prefix + '_G.h5')
    G = Gf(mesh=wk_mesh, target_shape=[nstates, nstates])
    G = fly.diagram.Diagram(G, 'Fermion')
    G.load(G_data)

    # Load singlet pairing vertex (not the FLEX vertex used for self-energy)
    #vertex = fly.Field_C(outdir + prefix + '_vertex_singlet.h5')
    vertex = fly.Field_CM(outdir + prefix + '_vertex.h5')
    DLRImMesh = fly.load_triqs_H.create_dlr_meshes(e_k, beta, statistic='Boson')
    wk_mesh = MeshProduct(DLRImMesh, k_mesh)
    V = Gf(mesh=wk_mesh, target_shape=[nstates, nstates, nstates, nstates])
    V = fly.diagram.Diagram(V, 'Boson')
    #fly.interface_triqs.fill_triqs_from_field(V.obj_wk, vertex)
    V.load(vertex)

    # Check vertex structure before symmetrization
    print_vertex_structure(V, "Before symmetrization")

    # Note: test.py doesn't symmetrize the vertex, so commenting this out
    symmetrize_vertex(V)

    # Check after symmetrization
    print_vertex_structure(V, "After symmetrization")

    V.wk_to_tr()

    Delta0 = G.copy()

    eig, Delta = solve_eliashberg_power_iteration(G, V, Delta0)
    print(f"Max Eig: {eig:.6f}")

    # Analyze gap symmetry
    analyze_gap_symmetry(Delta)

    Delta.save(outdir + prefix + '_gap.h5')


def solve_eliashberg_power_iteration(G, V, Delta0):
    max_iter = 100
    tol = 1e-4
    Delta = Delta0.copy()

    # Properly flip G(k, iω) → G(-k, -iω) for Cooper pair formation
    G_flip = flip_wk(G)

    # Create d-wave initial guess: cos(2πkx) - cos(2πky) (matching test.py)
    # Extract mesh dimensions
    mesh, BZ = extract_mesh_and_bz(Delta.obj_wk)

    eig = 0.0
    prev_eig = 0.0
    old_Deltas = []
    max_eigs_searched = 5

    while eig <= 0.0 and len(old_Deltas) < max_eigs_searched:
        Delta.init_tail()
        Delta.obj_wk.data[:] = np.random.rand(*Delta.obj_wk.data.shape) + 1j * np.random.rand(*Delta.obj_wk.data.shape)
        #Delta.obj_wk.data[:] = 1
        iter = 0
        for it in range(max_iter):
            Delta_new = Eliashberg_step(G, G_flip, V, Delta)
            norm = np.sum(Delta.obj_wk.data * np.conj(Delta.obj_wk.data)).real
            eig = np.sum(Delta_new.obj_wk.data * np.conj(Delta.obj_wk.data)).real / norm
            # Normalize
            Delta_new.obj_wk.data[:] = project_out(Delta_new.obj_wk.data, old_Deltas)

            # Check convergence
            diff = np.abs(eig - prev_eig)
            prev_eig = eig
            norm = np.sum(Delta_new.obj_wk.data * np.conj(Delta_new.obj_wk.data)).real
            Delta_new.obj_wk.data[:] = Delta_new.obj_wk.data / norm
            Delta = Delta_new.copy()
            print(f"Eig: {eig} Error = {diff:.6e}")
            iter = it
            if diff < tol or np.isnan(diff):
                break
        print(f"Iterations: {iter+1}")
        print(f"eig{len(old_Deltas)} = {eig}")
        old_Deltas.append(Delta.obj_wk.data.copy())
    return eig, Delta

def Eliashberg_step(G, G_flip, V, Delta):
    F = Delta.copy()
    # Use conj(G) to match test.py's linearized gap equation formula
    # F = -G(k,iω) * conj(G(k,iω)) * Δ(k,iω) = -|G(k,iω)|² * Δ(k,iω)
    F.obj_wk.data[:] = -1.0 * G.obj_wk.data * np.conj(G.obj_wk.data) * Delta.obj_wk.data
    F.wk_to_tr()
    Delta_new = dot_tr(V, F)
    Delta_new.tr_to_wk()
    return Delta_new

def project_out(v, eigvecs):
    """Project out previously found eigenvectors using Gram-Schmidt orthogonalization."""
    for x in eigvecs:
        # Compute projection: proj = (x·v / x·x) * x
        # Use vdot for proper complex conjugation: vdot(a,b) = sum(conj(a) * b)
        v_flat = v.flatten()
        x_flat = x.flatten()
        proj = (np.vdot(x_flat, v_flat) / np.vdot(x_flat, x_flat)) * x
        v = v - proj

    # Normalize
    nv = np.linalg.norm(v)
    if nv > 0:
        v = v / nv
    else:
        raise ValueError("Deflation resulted in zero vector")
    return v


