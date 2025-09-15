import numpy as np
import sparse_ir
from ffirefly import *

[nk1, nk2, nk3] = k_mesh
nk = nk1*nk2*nk3

class Eliashberg:
    def __init__(self, mesh, sfc_tol=1e-4, maxiter=100, mix=0.2, verbose=True):
        self.mesh = mesh
        self.sfc_tol = sfc_tol
        self.maxiter = maxiter
        self.mix = mix
        self.verbose = verbose

        self.phi = np.tensordot(np.ones(len(self.mesh.iwn_f)), self.mesh, axes=0)
        self.Z = np.tensordot(np.ones(len(self.mesh.iwn_f)), self.mesh, axes=0)

        self.V_calc()
        self.phi_Z_fr_calc()
        self.phi_Z_calc()
    
    
    #%%%%%%%%%%% Loop solving instance
    def solve(self):
        # perform loop until convergence is reached:
        for it in range(self.maxiter):
            phi_old = self.phi
            self.loop()

            # check whether solution is converged.
            sfc_check = np.sum(abs(self.phi-phi_old))/np.sum(abs(self.phi))
            if self.verbose:
                print(it, sfc_check)
            if sfc_check < self.sfc_tol:
                print("Eliashberg loop converged at desired accuracy")
                break
    
    def loop(self):
        phi_old = self.phi
        Z_old = self.Z
        
        self.phi_Z_fr_calc()
        self.phi_Z_calc()
        self.phi = self.mix*self.phi + (1-self.mix)*phi_old
        self.Z = self.mix*self.Z + (1-self.mix)*Z_old
        

    def phi_Z_fr_calc(self):
        # Fourier transform
        phi_fr = self.mesh.k_to_r(self.phi)
        self.phi_fr = self.mesh.wn_to_tau('F', phi_fr)
        Z_fr = self.mesh.k_to_r(self.Z)
        self.Z_fr = self.mesh.wn_to_tau('F', Z_fr)

    def phi_Z_calc(self):
        phi_fr = self.V * self.phi_fr[::-1, :] / (self.phi_fr[::-1, :]**2 + self.Z_fr[::-1, :]**2)
        Z_fr = self.V * self.Z_fr[::-1, :] / (self.phi_fr[::-1, :]**2 + self.Z_fr[::-1, :]**2)

        # Fourier transform
        phi = self.mesh.r_to_k(phi_fr)
        self.phi = self.mesh.tau_to_wn('B', phi)
        Z = self.mesh.r_to_k(Z_fr)
        self.Z = self.mesh.tau_to_wn('B', Z)

    def V_calc(self):
        
        V = get_V(self.mesh)

        # Fourier transform
        V = self.mesh.k_to_r(V)
        self.V = self.mesh.wn_to_tau('B', V)

