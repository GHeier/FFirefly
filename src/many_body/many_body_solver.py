from triqs_tprf.lattice import *
from triqs_tprf import *
from triqs.gf import inverse, Fourier, Gf, make_gf_dlr, make_gf_dlr_imtime, make_gf_dlr_imfreq
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime
import numpy as np
from scipy.optimize import brentq

from triqs.plot.mpl_interface import *
import matplotlib.pyplot as plt
import firefly as fly
from firefly.diagram import Diagram, dot_tr
import firefly.config as cfg

class FLEX_DMFT_Solver:
    def __init__(self, G0, U=0.0, mix=0.2, U_maxiter=50, n=None, mu=None):
        self.U = U
        self.UX = 0.0  # Track U * max(chi)
        self.mix = mix  # Mixing parameter for self-energy
        self.U_maxiter = U_maxiter  # Max iterations for U renormalization
        self.diverged = False  # Track divergence state
        self.n = n   # Target electron density (per spin)
        self.mu = mu

        if G0 is None:
            raise ValueError("Must provide G0 triqs object")

        self.G0 = Diagram(G0, 'Fermion')
        if self.G0.varspace == 'wk':
            stuff = G0.mesh.components[0]
        elif self.G0.varspace == 'w':
            stuff = G0.mesh
        else:
            raise ValueError("G0 must be in 'wk' or 'w' space")
        self.beta = stuff.beta
        self.eps = stuff.eps
        self.w_max = stuff.w_max

        IPT = IPTSolver(self.beta, None, mix=mix, w_max=self.w_max, eps=self.eps)
        FLEX = FLEXSolver(self.G0, U, mix=mix, U_maxiter=U_maxiter, n=n, mu=mu)

        # FLEX+DMFT
        DLR_f = MeshDLRImFreq(beta=self.beta, statistic='Fermion', w_max=self.w_max, eps=self.eps)
        DLR_b = MeshDLRImFreq(beta=self.beta, statistic='Boson', w_max=self.w_max, eps=self.eps)

        sigma = Gf(mesh=DLR_f, target_shape=self.G0.obj_wk.target_shape)
        self.Sigma_loc = Diagram(sigma, 'Fermion')
        self.Sigma_imp = self.Sigma_loc.copy()
        self.Sigma_nonloc = self.G0.copy()
        self.Sigma_nonloc.zero()
        self.Sigma_k = self.G0.copy()
        self.Sigma_k.zero()

        self.X_loc = Diagram(Gf(mesh=DLR_b, target_shape=(self.G0.obj_wk.mesh.nw)), 'Boson')

#    # DMFT methods
#    def get_IPT_sigma(self):
#        """Single DMFT iteration using IPT."""
#        # Transform Weiss field to imaginary time
#        self.G_loc.w_to_t()
#
#        # IPT: Sigma(tau) = U^2 * G_Weiss(tau)^3
#        self.Sigma_loc.obj_t.data[:] = (self.U**2) * self.G_loc.obj_t.data**3
#
#        Sigma_old = self.Sigma_loc.obj_w.copy()
#        # Transform Sigma back to frequency
#        self.Sigma_loc.t_to_w()
#        # Mix self-energy for stability
#        self.Sigma_loc.obj_w.data[:] = self.mix * self.Sigma_loc.obj_w.data + (1.0 - self.mix) * Sigma_old.data
#
#
#    def solve_DMFT(self):
#        self.get_IPT_sigma()
#        # Dyson equation: G_loc = H(Sigma)
#        self.G_loc.obj_w << self.H(Sigma=self.Sigma_loc.obj_w, mu=self.mu)
#
#        # Self-consistency: G_Weiss^-1 = G_loc^-1 + Sigma
#        self.G_loc.obj_w << inverse(inverse(self.G_loc.obj_w) + self.Sigma_loc.obj_w)
#        self.shift_mu_to_target_density(self.n)
#
#    def loop_DMFT(self, n_loops=100, tol=1e-6):
#        """Self-consistent DMFT loop."""
#        print("Beginning DMFT Self-Consistent Loop")
#        for i in range(n_loops):
#            G_old = self.G_loc.obj_w.copy()
#
#            self.solve_DMFT()
#
#            err = np.max(np.abs(self.G_loc.obj_w.data - G_old.data))
#            print(f"DMFT loop {i+1}, err = {err:.3e}")
#
#            if err < tol:
#                print("Convergence achieved.")
#                break
#
    def add_local_to_nonlocal(A, B):
        temp = A.copy()
        temp.obj_wk.data[:] += B.obj_w.data[:, np.newaxis, :, :]
        return temp

    def subtract_local_from_nonlocal(A, B):
        temp = A.copy()
        temp.obj_wk.data[:] -= B.obj_w.data[:, np.newaxis, :, :]
        return temp

    def get_local_chi(self, G_loc):
        G_loc.w_to_t()
        self.X_loc.obj_t.data[:] = G_loc.obj_t.data * G_loc.obj_t.data
        self.X_loc.t_to_w()

    def get_local_sigma(V, G):
        V.w_to_t()
        G.w_to_t()
        self.Sigma_loc = dot_t(V, G)
        self.Sigma_loc.t_to_w()

    def solve_FLEX_DMFT(self, FLEX, IPT):
        """Single iteration of FLEX+DMFT."""
        # DMFT part
        FLEX.dyson_G_from_Sigma()
        FLEX.get_local_G()
        IPT.G0_iw = FLEX.G_loc.copy()
        sigma_imp = IPT.get_IPT_sigma()
        self.Sigma_imp.obj_w = sigma_imp

        # Combine
        self.Sigma_k = add_local_to_nonlocal(self.Sigma_nonloc, self.Sigma_imp)
        FLEX.Sigma = self.Sigma_k

        # Make G
        mu = FLEX.find_mu_for_density(self.n)
        FLEX.make_G(mu)
        FLEX.get_local_G()

        # Local + Nonlocal parts for FLEX
        FLEX.chi0_from_grt_PH()
        self.get_local_chi(FLEX.G_loc)

        V_nonloc = FLEX.FLEX_from_chi(FLEX.X.obj_wk)
        FLEX.V = V_nonloc
        V_loc = FLEX.FLEX_from_chi(self.X_loc.obj_w)
        if FLEX.diverged:
            print("Divergence detected in FLEX part of FLEX+DMFT")
            return

        FLEX.Sigma_from_vertex()
        self.get_local_sigma(FLEX.V, FLEX.G_loc)
        self.sigma_nonloc = subtract_local_from_nonlocal(FLEX.Sigma, self.Sigma_loc)

        # Combine
        self.Sigma_k = add_local_to_nonlocal(self.sigma_nonloc, self.Sigma_imp)
        FLEX.Sigma_k = self.Sigma_k.copy()

    def U_renormalization(self, FLEX, IPT):
        """Renormalize U if U*max(chi) >= 1 to avoid divergence."""
        print("WARNING: U is too large and the spin susceptibility will diverge!")
        print("Initiating U renormalization loop...")

        U_old = self.U
        U_it = 0
        prev_U = 0
        # Check condition: U_old * max(chi0) >= 1
        while U_old * np.max(np.abs(FLEX.X.obj_wk.data)) >= 1.0:
            U_it += 1

            # Reduce U temporarily to bring UX below 1
            max_X = np.max(np.abs(FLEX.X.obj_wk.data))
            self.U = self.U / (max_X * self.U + 0.01)
            FLEX.U = self.U
            IPT.U = self.U
            print(f"{U_it}) U = {self.U}, initial_U = {U_old}")

            # Perform one FLEX loop iteration with reduced U (matching test.py logic)
            G_old = self.G.obj_wk.copy()

            # Calculate V and Sigma with current chi0 and reduced U
            self.solve_FLEX_DMFT(FLEX, IPT)

            # Reset U back to U_old for next iteration
            diff = abs(prev_U - self.U)
            prev_U = self.U
            self.U = U_old

            if U_it >= self.U_maxiter or diff < 1e-3:
                print(f"U_diff = {diff:.4e}, tol = 1e-3")
                print(f"Iteration number {U_it}, max iterations {self.U_maxiter}")
                break

        print("Leaving U renormalization...")
        # Final UX calculation with U_old
        self.UX = self.U * np.max(np.abs(self.X.obj_wk.data))
        FLEX.U = self.U
        IPT.U = self.U

    def loop_FLEX_DMFT(self, n_loops=50, check_divergence=True):
        # Initial check for divergence
        if check_divergence and self.UX >= 1.0:
            print(f"Initial U*max(Chi) = {self.UX:.4f} >= 1")
            self.U_renormalization()

        print("Beginning FLEX+DMFT Self-Consistent Loop")
        for i in range(n_loops):
            #print(f"Starting iteration {i+1}/{n_loops}...")
            G_old = self.G.obj_wk.copy()

            self.solve_FLEX_DMFT()
            if self.diverged:
                print(f"Divergence detected at iteration {i+1}")
                if check_divergence:
                    print("Attempting U renormalization...")
                    self.U_renormalization()
                    continue
                else:
                    print("Code exiting due to divergence.")
                    break

            print(f"Iteration {i+1} completed. U*max(Chi) = {self.UX:.4f}", end=' ')
            err = np.max(np.abs(self.G.obj_wk.data - G_old.data))
            print(f"Max change in G: {err:.3e}")

            if err < 1e-6:
                print("Convergence achieved.")
                break
