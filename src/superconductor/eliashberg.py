import firefly as fly
import firefly.config as cfg
from firefly.diagram import *

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

Nk = cfg.k_mesh[0]
BZ = get_brillouin_zone()

mu = cfg.fermi_energy
beta = 1.0 / cfg.Temperature

def main():
    H_r, kmesh, e_k = fly.load_triqs_H.get_energy_mesh()
    emax = e_k.data.max().real
    emin = e_k.data.min().real
    print(f"emax: {emax}, emin: {emin}")
    wmax = 1.2 * (emax - emin)
    DLRImMesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=wmax, eps=1e-14)
    print(BZ)
    k_mesh = MeshBrZone(BZ, n_k=Nk)   # uniform Nk x Nk x Nk (third dim is 1 if 2D)

    sigma = fly.Field_C(outdir + prefix + '_sigma.h5')
    wk_mesh = MeshProduct(DLRImMesh, k_mesh) 
    E = Gf(mesh=wk_mesh, target_shape=[1,1])
    E.data[:, :, 0, 0] = sigma.get_data()
    #fly.interface_triqs.fill_triqs_from_field(E, sigma)
    G0 = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=DLRImMesh)
    G = inverse(inverse(G0) - E)
    G = fly.diagram.Diagram(G, 'Fermion')

    vertex = fly.Field_C(outdir + prefix + '_vertex.h5')
    DLRImMesh = MeshDLRImFreq(beta=beta, statistic='Boson', w_max=wmax, eps=1e-14)
    wk_mesh = MeshProduct(DLRImMesh, k_mesh) 
    V = Gf(mesh=wk_mesh, target_shape=[1,1,1,1])
    V = fly.diagram.Diagram(V, 'Boson')
    #fly.interface_triqs.fill_triqs_from_field(V.obj_wk, vertex)
    V.obj_wk.data[:, :, 0, 0, 0, 0] = vertex.get_data()
    V.save(outdir + prefix + '_vertex.h5')
    return
    V.obj_wk.data[:, :, 0, 0, 0, 0] = (V.obj_wk.data[:, :, 0, 0, 0, 0] + np.flip(V.obj_wk.data[:, :, 0, 0, 0, 0], axis=1)) / 2.0
    V.wk_to_tr()

    Delta0 = G.copy()

    eig, Delta = solve_eliashberg_power_iteration(G, V, Delta0)
    print(f"Leading eigenvalue: {eig:.6f}")
    Delta.save(outdir + prefix + '_gap.h5')


def solve_eliashberg_power_iteration(G, V, Delta0):
    max_iter = 100
    tol = 1e-4
    max_eigs_searched = 5
    Delta = Delta0.copy()
    G_flip = G.copy()
    # reverse G(r,tau) to G(-r, -tau)
    G_flip.obj_wk.data[:] = np.flip(np.flip(G_flip.obj_wk.data, axis=0), axis=1)
    #G_flip.tr_to_wk()
    eig = 0.0
    prev_eig = 0.0
    old_Deltas = []

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
    F.obj_wk.data[:] = -1.0 * G.obj_wk.data * G_flip.obj_wk.data * Delta.obj_wk.data
    #G2 = G.obj_wk.data * G_flip.obj_wk.data
    #print(f"Max G(iw,k) = {np.max(G2)}")
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

