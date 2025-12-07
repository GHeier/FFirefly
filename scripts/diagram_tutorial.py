
# Check if TRIQS is available
from triqs.gf import Gf
from triqs_tprf.tight_binding import TBLattice
from triqs.gf.meshes import MeshDLRImFreq
from triqs.gf.mesh_product import MeshProduct
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice import chi0_tr_from_grt_PH

import numpy as np
import firefly as fly
Diagram = fly.diagram.Diagram

# Temporary file for saving
file = "temp_sigma.h5"

def TB_2D_square_lattice(t=1.0):
    # Define lattice vectors and hopping
    H_r = TBLattice(
        units=[
            (1, 0, 0),  # basis vector in the x-direction
            (0, 1, 0),  # basis vector in the y-direction
        ],
        hoppings={
            (+1, 0): [[-t]],   # nearest-neighbor hopping in +x
            (-1, 0): [[-t]],   # nearest-neighbor hopping in -x
            (0, +1): [[-t]],   # nearest-neighbor hopping in +y
            (0, -1): [[-t]],   # nearest-neighbor hopping in -y
        }
    )

    # Create k-mesh and compute dispersion
    Nk = 100
    kmesh = H_r.get_kmesh(n_k=Nk)
    e_k = H_r.fourier(kmesh)

    return H_r, kmesh, e_k

if __name__ == "__main__":
    beta = 4.0  # Inverse temperature
    H_r, kmesh, e_k = TB_2D_square_lattice(t=1.0)
    bandwidth = np.max(e_k.data).real - np.min(e_k.data).real
    print("Bandwidth: ", bandwidth)
    # Discrete Lehman Representation Mesh. 
    # Uses functions to dramatically reduce storage size 
    # The cost of the reduction is not being able to store constant values (ie V = const)
    iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=1.2*bandwidth, eps=1e-6)
    iv_mesh = MeshDLRImFreq(beta=beta, statistic='Boson', w_max=1.2*bandwidth, eps=1e-6)

    # Creates non-interacting Green's function using TRIQS library
    G_init = lattice_dyson_g0_wk(mu=0.0, e_k=e_k, mesh=iw_mesh)
    # Converts this to diagram form
    G0 = Diagram(G_init, 'Fermion')
    G0.wk_to_tr() # Transforms G(iw, k) to G(tau, r)

    # Perform simple diagrammatic operation: X0(iv, q) = T ∑_m ∫ dk G0(iw+iv, k+q) * G0(iw, k)
    chi_tr = chi0_tr_from_grt_PH(G0.obj_tr) # X0(tau, r) = G0(tau, r) * G0(-tau, -r)
    # Converts this to diagram form
    X0 = Diagram(chi_tr, 'Boson')
    X0.tr_to_wk() # Transforms X0(tau, r) to X0(iw, q)

    # Make Local Gren's function G(iw)
    G_w = Gf(mesh=iw_mesh, target_shape=[1,1])
    obj_w = np.sum(G0.obj_wk.data, axis=1) / G0.nk  # Sum over k-points
    G_w.data[:] = obj_w
    G_ipt = Diagram(G_w, 'Fermion')

    # Compute IPT self-energy: Sigma(iw) = ∫dtau U^2 * G(tau)^3
    U = 2.0
    Sigma = G_ipt.copy()
    G_ipt.w_to_t()
    Sigma.obj_t << (U**2) * G_ipt.obj_t * G_ipt.obj_t * G_ipt.obj_t
    Sigma.t_to_w()

    print("Max X0(iw, q): ", np.max(np.abs(X0.obj_wk.data)))
    print("Expected: ~0.4")

    print("Max Sigma(iw): ", np.max(np.abs(Sigma.obj_w.data)))
    print("Expected: ~0.1")
    X0.save("X0_iwq.h5")
    Sigma.save("Sigma_ipt.h5")


    # Testing
    Delta = Gf(mesh=kmesh, target_shape=[1,1])
    print(kmesh)
    print("kmesh components: ")
    print(kmesh.dims)
    Delta = Diagram(Delta, 'Fermion')
    Delta.obj_k.data[:] = 1.0
    Delta.k_to_r()
    print("Delta(r) data: ", Delta.obj_r.data)
