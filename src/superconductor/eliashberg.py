import fcode
import fcode.config as cfg
from fcode import *

import numpy as np
import sparse_ir

BZ = np.array(cfg.brillouin_zone)
U = cfg.onsite_U

T = cfg.Temperature
beta = 1 / T
nk1, nk2, nk3 = cfg.k_mesh
nk = nk1 * nk2 * nk3


def get_bandwidth():
    k1, k2, k3 = np.meshgrid(np.arange(nk1)/nk1, np.arange(nk2)/nk2, np.arange(nk3)/nk3)
    emin, emax = 1000, -1000
    # Flatten the k-points into (N, 3) array
    k_points = np.stack((k1.ravel(), k2.ravel(), k3.ravel()), axis=-1)
    k_points = k_points @ BZ.T
    e_pts = fcode.epsilon(1, k_points)
    e_mesh = e_pts.reshape(nk1, nk2, nk3)
    emin = np.min(e_mesh)
    emax = np.max(e_mesh)
    return emax - emin

def eliashberg():
    W = get_bandwidth()
    W = round(W, 6)
    print("Bandwidth:", W)

    wmax = 2*W
    IR_tol = 1e-5
    IR_basis_set = sparse_ir.FiniteTempBasisSet(beta, wmax, eps=IR_tol)

    result = EliashbergSolver(IR_basis_set, nk1, nk2, nk3)
    result.solve()

class EliashbergSolver:
    def __init__(self, IR_basis_set, nk1, nk2, nk3):
        self.maxiters = 10

        self.mesh = Mesh(IR_basis_set, nk1, nk2, nk3)
        self.phi = np.ones(self.mesh.iwn_f.shape)
        self.Z = np.ones(self.mesh.iwn_f.shape)

        self.file = "/home/g/Research/fcode/" + cfg.input_data_file[1:-1]
        print("Reading susceptibility from", self.file)
        #self.chi0 = Susceptibility(file, 3)
        self.chi = field.load_field_from_file(self.file, 3, False, False)

        print(len(self.mesh.pts))
        print(self.mesh.pts.shape)
        self.e_flat = self.mesh.ek.ravel()

        self.coords = self.mesh.pts[:, None, :] - self.mesh.pts[None, :, :]
        self.chi_vals = self.chi(self.coords)
        print(type(self.chi_vals))
        self.V_mesh = U**2 * self.chi_vals / (1 - U * self.chi_vals) + U**3 * self.chi_vals**2 / (1 - U * self.chi_vals)**2

    def solve(self):
        for i in range(self.maxiters):
            self.denom = (self.w_flat * self.Z)**2 + self.e_flat**2 + self.phi**2
            self.phi = (1 / beta) * np.dot(self.V_mesh, self.phi) / self.denom
            self.Z = 1 + (1 / beta) * np.dot(self.V_mesh, self.Z) / self.denom * (self.w_flat[:, None] / self.w_flat).T

    def save_results(self):
        self.phi = self.phi.reshape(self.mesh.iwn_f.shape)
        self.Z = self.Z.reshape(self.mesh.iwn_f.shape)
        self.phi_tau = self.mesh.wn_to_tau('F', self.phi)
        self.Z_tau = self.mesh.wn_to_tau('F', self.Z)
        self.phi_r = self.mesh.k_to_r(self.phi)
        self.Z_r = self.mesh.k_to_r(self.Z)
        np.savez("results.npz", phi=self.phi, Z=self.Z, phi_tau=self.phi_tau, Z_tau=self.Z_tau, phi_r=self.phi_r, Z_r=self.Z_r)


class Mesh:
    def __init__(self,IR_basis_set,nk1,nk2,nk3):
        self.IR_basis_set = IR_basis_set

        # generate k-mesh and dispersion
        self.nk1, self.nk2, self.nk3, self.nk = nk1, nk2, nk3, nk1*nk2*nk3
        self.k1, self.k2, self.k3 = np.meshgrid(np.arange(self.nk1)/self.nk1, np.arange(self.nk2)/self.nk2, np.arange(self.nk3)/self.nk3)
        self.ek = mesh.energy_mesh([self.nk1, self.nk2, self.nk3])

        # lowest Matsubara frequency index
        self.iw0_f = np.where(self.IR_basis_set.wn_f == 1)[0][0]

        ### Generate a frequency-momentum grid for iwn and ek (in preparation for calculating the Green function)
        # frequency mesh (for Green function)
        self.iwn_f = 1j * self.IR_basis_set.wn_f * np.pi * T
        self.iwn_f_ = np.tensordot(self.iwn_f, np.ones(self.nk), axes=0)

        # ek mesh
        self.ek_ = np.tensordot(np.ones(len(self.iwn_f)), self.ek, axes=0)
        self.pts = np.meshgrid(np.arange(self.nk1), np.arange(self.nk2), np.arange(self.nk3), self.iwn_f)
        self.pts = np.vstack([x.ravel() for x in self.pts]).T

    def smpl_obj(self, statistics):
        smpl_tau = {'F': self.IR_basis_set.smpl_tau_f, 'B': self.IR_basis_set.smpl_tau_b}[statistics]
        smpl_wn  = {'F': self.IR_basis_set.smpl_wn_f,  'B': self.IR_basis_set.smpl_wn_b }[statistics]
        return smpl_tau, smpl_wn

    
    def tau_to_wn(self, statistics, obj_tau):
        smpl_tau, smpl_wn = self.smpl_obj(statistics)

        obj_tau = obj_tau.reshape((smpl_tau.tau.size, self.nk1, self.nk2, self.nk3))
        obj_l   = smpl_tau.fit(obj_tau, axis=0)
        obj_wn  = smpl_wn.evaluate(obj_l, axis=0).reshape((smpl_wn.wn.size, self.nk))
        return obj_wn

    def wn_to_tau(self, statistics, obj_wn):
        smpl_tau, smpl_wn = self.smpl_obj(statistics)

        obj_wn  = obj_wn.reshape((smpl_wn.wn.size, self.nk1, self.nk2, self.nk3))
        obj_l   = smpl_wn.fit(obj_wn, axis=0)
        obj_tau = smpl_tau.evaluate(obj_l, axis=0).reshape((smpl_tau.tau.size, self.nk))
        return obj_tau

    
    def k_to_r(self,obj_k):
        obj_k = obj_k.reshape(-1, self.nk1, self.nk2, self.nk3)
        obj_r = np.fft.fftn(obj_k,axes=(1,2))
        obj_r = obj_r.reshape(-1, self.nk)
        return obj_r

    def r_to_k(self,obj_r):
        obj_r = obj_r.reshape(-1, self.nk1, self.nk2, self.nk3)
        obj_k = np.fft.ifftn(obj_r,axes=(1,2))/self.nk
        obj_k = obj_k.reshape(-1, self.nk)
        return obj_k
eliashberg()
