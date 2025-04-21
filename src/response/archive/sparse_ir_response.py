import numpy as np
import scipy as sc
import scipy.optimize
from warnings import warn
import sparse_ir
import matplotlib.pyplot as plt### System parameters

import ffirefly
import ffirefly.config as cfg
from ffirefly import *


class Mesh:
    """
    Holding class for k-mesh and sparsely sampled imaginary time 'tau' / Matsubara frequency 'iwn' grids.
    Additionally it defines the Fourier transform routines 'r <-> k'  and 'tau <-> l <-> wn'.
    """
    def __init__(self,IR_basis_set,e_mesh):
        T = cfg.Temperature
        self.IR_basis_set = IR_basis_set
        self.nk1, self.nk2, self.nk3 = e_mesh.shape
        self.ek = e_mesh.flatten()
        self.nk = self.nk1 * self.nk2 * self.nk3

        # lowest Matsubara frequency index
        self.iw0_f = np.where(self.IR_basis_set.wn_f == 1)[0][0]
        self.iw0_b = np.where(self.IR_basis_set.wn_b == 0)[0][0]

        self.iwn_f = 1j * self.IR_basis_set.wn_f * np.pi * T
        self.iwn_f_ = np.tensordot(self.iwn_f, np.ones(self.nk), axes=0)

        # ek mesh
        self.ek_ = np.tensordot(np.ones(len(self.iwn_f)), self.ek, axes=0)

    def smpl_obj(self, statistics):
        """ Return sampling object for given statistic """
        smpl_tau = {'F': self.IR_basis_set.smpl_tau_f, 'B': self.IR_basis_set.smpl_tau_b}[statistics]
        smpl_wn  = {'F': self.IR_basis_set.smpl_wn_f,  'B': self.IR_basis_set.smpl_wn_b }[statistics]
        return smpl_tau, smpl_wn

    
    def tau_to_wn(self, statistics, obj_tau):
        """ Fourier transform from tau to in via IR basis """
        smpl_tau, smpl_wn = self.smpl_obj(statistics)

        obj_tau = obj_tau.reshape((smpl_tau.tau.size, self.nk1, self.nk2, self.nk3))
        obj_l   = smpl_tau.fit(obj_tau, axis=0)
        obj_wn  = smpl_wn.evaluate(obj_l, axis=0).reshape((smpl_wn.wn.size, self.nk))
        return obj_wn

    def wn_to_tau(self, statistics, obj_wn):
        """ Fourier transform from iwn to tau via IR basis """
        smpl_tau, smpl_wn = self.smpl_obj(statistics)

        obj_wn  = obj_wn.reshape((smpl_wn.wn.size, self.nk1, self.nk2, self.nk3))
        obj_l   = smpl_wn.fit(obj_wn, axis=0)
        obj_tau = smpl_tau.evaluate(obj_l, axis=0).reshape((smpl_tau.tau.size, self.nk))
        return obj_tau

    
    def k_to_r(self,obj_k):
        """ Fourier transform from k-space to real space """
        obj_k = obj_k.reshape(-1, self.nk1, self.nk2, self.nk3)
        obj_r = np.fft.fftn(obj_k,axes=(1,2))
        obj_r = obj_r.reshape(-1, self.nk)
        return obj_r

    def r_to_k(self,obj_r):
        """ Fourier transform from real space to k-space """
        obj_r = obj_r.reshape(-1, self.nk1, self.nk2, self.nk3)
        obj_k = np.fft.ifftn(obj_r,axes=(1,2))/self.nk
        obj_k = obj_k.reshape(-1, self.nk)
        return obj_k

class FLEXSolver:
    def __init__(self, mesh, U, n, sigma_init=0, sfc_tol=1e-4, 
                 maxiter=100, U_maxiter=10, mix=0.2, verbose=True):
        """
        Solver class to calculate the FLEX loop self-consistently.
        After initializing the Solver by `solver = FLEXSolver(mesh, U, n, **kwargs)` 
        it can be run by `solver.solve()`.
        """
        ## Set internal parameters for the solve 
        self.U = U
        self.n = n
        self.mesh = mesh
        self.sigma = sigma_init
        self.sfc_tol = sfc_tol
        self.maxiter = maxiter
        self.U_maxiter = U_maxiter
        self.mix = mix
        self.verbose = verbose
        
        ## Set initial Green function and irreducible susceptibility
        # NOT running the FLEXSolver.solve instance corresponds to staying on RPA level
        self.mu = 0
        self.mu_calc()
        
        self.gkio_calc(self.mu)
        self.grit_calc()
        self.ckio_calc()
    
    
    #%%%%%%%%%%% Loop solving instance
    def solve(self):
        """ FLEXSolver.solve() executes FLEX loop until convergence """
        # check whether U < U_crit! Otherwise, U needs to be renormalized.
        if np.amax(np.abs(self.ckio))*self.U >= 1:
            self.U_renormalization()
            
        # perform loop until convergence is reached:
        for it in range(self.maxiter):
            sigma_old = self.sigma
            self.loop()

            # check whether solution is converged.
            sfc_check = np.sum(abs(self.sigma-sigma_old))/np.sum(abs(self.sigma))
            if self.verbose:
                print(it, sfc_check)
            if sfc_check < self.sfc_tol:
                print("FLEX loop converged at desired accuracy")
                break
    
    def loop(self):
        """ FLEX loop """
        gkio_old = self.gkio
        
        # calculate interaction and self-energy
        self.V_calc()
        self.sigma_calc()
        
        # set new chemical potential and apply mixing
        self.mu_calc()
        self.gkio_calc(self.mu)
        self.gkio = self.mix*self.gkio + (1-self.mix)*gkio_old
        
        # calculate new irreducible susceptibility
        self.grit_calc()
        self.ckio_calc()


    #%%%%%%%%%%% U renormalization loop instance
    def U_renormalization(self):
        """ Loop for renormalizing U if Stoner enhancement U*max{chi0} >= 1. """
        print('WARNING: U is too large and the spin susceptibility denominator will diverge/turn unphysical!')
        print('Initiate U renormalization loop.')
    
        # save old U for later
        U_old = self.U
        # renormalization loop may run infinitely! Insert break condition after U_it_max steps
        U_it = 0
    
        while U_old*np.amax(np.abs(self.ckio)) >= 1:
            U_it += 1
        
            # remormalize U such that U*chi0 < 1
            self.U = self.U / (np.amax(np.abs(self.ckio))*self.U + 0.01)
            print(U_it, self.U, U_old)
        
            # perform one shot FLEX loop
            self.loop()
        
            # reset U
            self.U = U_old
        
            # break condition for too many steps
            if U_it == self.U_maxiter:
                print('Iteration number of U renormalization reached break condition!')
                break
        print('Leaving U renormalization...')
    
    
    #%%%%%%%%%%% Calculation steps
    def gkio_calc(self, mu):
        """ calculate Green function G(iw,k) """
        self.gkio = (self.mesh.iwn_f_ - (self.mesh.ek_ - mu) - self.sigma)**(-1)

    def grit_calc(self):
        """ Calculate real space Green function G(tau,r) [for calculating chi0 and sigma] """
        # Fourier transform
        grit = self.mesh.k_to_r(self.gkio)
        self.grit = self.mesh.wn_to_tau('F', grit)

    def ckio_calc(self):
        """ Calculate irreducible susciptibility chi0(iv,q) """
        ckio = self.grit * self.grit[::-1, :]

        # Fourier transform
        ckio = self.mesh.r_to_k(ckio)
        self.ckio = self.mesh.tau_to_wn('B', ckio)

    def V_calc(self):
        """ Calculate interaction V(tau,r) from RPA-like spin and charge susceptibility for calculating sigma """
        # check whether U is too large and give warning
        if np.amax(np.abs(self.ckio))*self.U >= 1:
            warn("U*max(chi0) >= 1! Paramagnetic phase is left and calculations will turn unstable!")
        
        # spin and charge susceptibility
        self.chi_spin   = self.ckio / (1 - self.U*self.ckio)
        self.chi_charge = self.ckio / (1 + self.U*self.ckio)

        V = 3/2*self.U**2 * self.chi_spin + 1/2*self.U**2 * self.chi_charge - self.U**2 * self.ckio
        # Constant Hartree Term V ~ U needs to be treated extra, since it cannot be modeled compactly by the IR basis.
        # In the single-band case, the Hartree term can be absorbed into the chemical interaction.

        # Fourier transform
        V = self.mesh.k_to_r(V)
        self.V = self.mesh.wn_to_tau('B', V)

    def sigma_calc(self):
        """ Calculate self-energy Sigma(iw,k) """
        sigma = self.V * self.grit
    
        # Fourier transform
        sigma = self.mesh.r_to_k(sigma)
        self.sigma = self.mesh.tau_to_wn('F', sigma)
    
    
    #%%%%%%%%%%% Setting chemical potential mu
    def calc_electron_density(self, mu):
        """ Calculate electron density from Green function """
        self.gkio_calc(mu)
        gio  = np.sum(self.gkio,axis=1)/self.mesh.nk
        g_l  = self.mesh.IR_basis_set.smpl_wn_f.fit(gio)
        g_tau0 = self.mesh.IR_basis_set.basis_f.u(0)@g_l
    
        n  = 1 + np.real(g_tau0)
        n  = 2*n #for spin
        return n

    def mu_calc(self):
        """ Find chemical potential for a given filling n0 via brent's root finding algorithm """
        n_calc = self.calc_electron_density
        n0 = self.n
        f  = lambda mu : n_calc(mu) - n0

        self.mu = sc.optimize.brentq(f, np.amax(self.mesh.ek)*3, np.amin(self.mesh.ek)*3)

    def save_ckio(self, filename):
        print('Saving chi0(k) to file:', filename)
        x_real = np.real(self.ckio)
        x_imag = np.imag(self.ckio)
        k1, k2, k3, w= np.meshgrid(np.arange(self.mesh.nk1) / self.mesh.nk1, np.arange(self.mesh.nk2) / self.mesh.nk2, np.arange(self.mesh.nk3) / self.mesh.nk3, np.imag(self.mesh.iwn_f[1:]), indexing='ij')
        print(self.mesh.iwn_f)
        print(w.shape)
        print(x_real.shape)
        k1 = k1.ravel()
        k2 = k2.ravel()
        k3 = k3.ravel()
        w = w.ravel()
        x_real = x_real.ravel()
        x_imag = x_imag.ravel()
        print(w.shape)
        print(x_real.shape)
        data = np.stack((k1, k2, k3, w, x_real, x_imag), axis=-1)
        np.savetxt(filename, data, fmt="%.8f")

def sparse_ir_response():
    T    = cfg.Temperature
    beta = 1/T    # inverse temperature
    n    = 0.85   # electron filling, here per spin per lattice site (n=1: half filling)
    U    = 1      # Hubbard interaction

    e_mesh = ffirefly.mesh.energy_mesh(cfg.k_mesh)
### Numerical parameters
    nk1, nk2, nk3 = cfg.k_mesh
    nk        = nk1*nk2*nk3
    wmax = 2 * (np.amax(e_mesh) - np.amin(e_mesh))
    IR_tol    = 1e-10     # desired accuary for l-cutoff of IR basis functions
    sfc_tol   = 1e-4      # desired accuracy for self-consistent iteration
    maxiter   = 30        # maximal number of iterations in self-consistent cycle
    mix       = 0.2       # mixing parameter for new green function
    U_maxiter = 50        # maximal number of iteration steps in U renormalization loop
# initialize calculation
    IR_basis_set = sparse_ir.FiniteTempBasisSet(beta, wmax, eps=IR_tol)
    mesh = Mesh(IR_basis_set, e_mesh)
    solver = FLEXSolver(mesh, U, n, sigma_init=0, sfc_tol=sfc_tol, maxiter=maxiter, U_maxiter=U_maxiter, mix=mix)

    solver.save_ckio('chi_dynamic_mesh.dat')

