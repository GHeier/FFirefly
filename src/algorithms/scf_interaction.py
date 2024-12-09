import sparse_ir_mesh
import sparse_ir
import numpy as np

from ..config.load.python_config import *
from ..hamiltonian.epsilon_module import epsilon
from ..hamiltonian.potential import V

def scf_interaction():
    beta = 1 / Temperature
    IR_tol = 1e-10
    sfc_tol = 1e-4
    nk1, nk2, nk3 = k_mesh
    wmax = 2*get_bandwidth(nk1, nk2, nk3)
    maxiter = 30
    mix = 0.2
    IR_basis_set = sparse_ir.FiniteTempBasisSet(beta, wmax, eps=IR_tol)
    mesh = sparse_ir_mesh.Mesh(IR_basis_set, nk1, nk2, nk3)
    solver = VSolver(mesh, U, n, sigma_init=0, sfc_tol=sfc_tol, maxiter=maxiter, mix=mix)

# perform FLEX loop
    solver.solve()
    save(self.sigma) # save the self-energy
    save(self.ckio) # save the susceptibility

def get_bandwidth(nk1, nk2, nk3):
    x, y, z = np.meshgrid(np.arange(nk1) / nk1, np.arange(nk2) / nk2, np.arange(nk3) / nk3)
    kx, ky, kz = brillouin_zone * np.array([x, y, z])
    ep_grid = epsilon(kx, ky, kz)
    return np.max(ep_grid) - np.min(ep_grid)
