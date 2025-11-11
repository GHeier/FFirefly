from triqs_tprf.lattice import *
from triqs_tprf import *
from triqs.gf import inverse, Gf
from triqs.gf.meshes import MeshDLRImFreq
from IPTSolver import *
from FLEXSolver import *
import numpy as np

from firefly.diagram import Diagram, dot_tr, dot_t
import firefly.config as cfg

class FLEX_DMFT_Solver:
    def __init__(self, G0, U=0.0, mix=0.2, U_maxiter=50, n=None, mu=None):
        self.U = U
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

        self.FLEX = FLEXSolver(self.G0.obj_wk, U, mix=mix, U_maxiter=U_maxiter, n=n, mu=mu)
        mu = self.FLEX.find_mu_for_density(n)
        self.IPT = IPTSolver(self.beta, None, mix=mix, w_max=self.w_max, eps=self.eps, mu=mu)

        # FLEX+DMFT
        DLR_f = MeshDLRImFreq(beta=self.beta, statistic='Fermion', w_max=self.w_max, eps=self.eps)
        DLR_b = MeshDLRImFreq(beta=self.beta, statistic='Boson', w_max=self.w_max, eps=self.eps)

        self.Sigma_loc = self.IPT.G_loc.copy()
        self.Sigma_imp = self.Sigma_loc.copy()
        self.Sigma_nonloc = self.G0.copy()
        self.Sigma_nonloc.zero()
        self.Sigma_k = self.Sigma_nonloc.copy()

        self.X_loc = Diagram(Gf(mesh=DLR_b, target_shape=(self.G0.obj_wk.target_shape)), 'Boson')

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

    def get_local_chi(self, G_loc):
        G_loc.w_to_t()
        self.X_loc.obj_t.data[:] = -G_loc.obj_t.data * G_loc.obj_t.data
        self.X_loc.t_to_w()

    def get_local_sigma(self, V, G):
        V.w_to_t()
        G.w_to_t()
        self.Sigma_loc = dot_t(V, G)
        self.Sigma_loc.t_to_w()

    def solve_FLEX_DMFT(self, FLEX, IPT):
        """Single iteration of FLEX+DMFT."""
        # DMFT part
        mu = FLEX.find_mu_for_density(self.n)
        FLEX.make_G(self.mu)
        #FLEX.dyson_G_from_Sigma()
        #FLEX.G = FLEX.G0.copy()
        FLEX.get_local_G()
        IPT.G_loc = FLEX.G_loc.copy()
        IPT.set_Weiss()
        # Update impurity self-energy with mixing
        sigma_ipt_old = self.Sigma_imp.obj_w.copy()
        self.Sigma_imp = IPT.get_IPT_Sigma(self.U)

        # Combine
        new_sigma_k = add_local_to_nonlocal(self.Sigma_nonloc, self.Sigma_imp)
        self.Sigma_k.obj_wk.data[:] = self.mix * new_sigma_k.obj_wk.data + (1.0 - self.mix) * self.Sigma_k.obj_wk.data
        FLEX.Sigma = self.Sigma_k.copy()

        # Make G
        mu = FLEX.find_mu_for_density(self.n)
        FLEX.make_G(mu)
        FLEX.get_local_G()

        # Nonlocal parts for FLEX
        FLEX.chi0_from_grt_PH()
        FLEX.V = FLEX.FLEX_from_chi(FLEX.X.obj_wk)
        if FLEX.diverged or self.U * np.max(np.abs(FLEX.X.obj_wk.data)) >= 1.0:
            self.diverged = True
            print("Divergence detected in FLEX part of FLEX+DMFT with U*max(Chi) = ", FLEX.UX)
            return

        FLEX.Sigma_from_vertex()
        # Local FLEX
        self.Sigma_loc.obj_w.data[:] = np.sum(FLEX.Sigma.obj_wk.data[:, :, :, :], axis=1) / FLEX.Sigma.nk
        # Update non-local self-energy 
        self.Sigma_nonloc = subtract_local_from_nonlocal(FLEX.Sigma, self.Sigma_loc)

        # Combine
        new_sigma_k = add_local_to_nonlocal(self.Sigma_nonloc, self.Sigma_imp)
        self.Sigma_k.obj_wk.data[:] = self.mix * new_sigma_k.obj_wk.data + (1.0 - self.mix) * self.Sigma_k.obj_wk.data
        FLEX.Sigma_k = self.Sigma_k.copy()

    def U_renormalization(self, FLEX, IPT):
        """Renormalize U if U*max(chi) >= 1 to avoid divergence."""
        print("WARNING: U is too large and the spin susceptibility will diverge!")
        print("Initiating U renormalization loop...")

        U_old = self.U
        print("Initial U = ", U_old)
        print("Initial max(Chi) = ", np.max(np.abs(FLEX.X.obj_wk.data)))
        print("Initial U*max(Chi) = ", U_old * np.max(np.abs(FLEX.X.obj_wk.data)))
        U_it = 0
        prev_U = 0
        # Check condition: U_old * max(chi0) >= 1
        while U_old * np.max(np.abs(FLEX.X.obj_wk.data)) >= 1.0:
        # UPDATE while np.max(np.abs(U_old * self.X.obj_wk.data)) >= 1.0:
            U_it += 1

            # Reduce U temporarily to bring UX below 1
            max_X = np.max(np.abs(FLEX.X.obj_wk.data))
            self.U = self.U / (max_X * self.U + 0.01)
            # UPDATE self.U = self.U / (np.max(np.abs(self.U * self.X.obj_wk.data)) + 0.01)
            FLEX.U = self.U
            IPT.U = self.U
            print(f"{U_it}) U = {self.U:.4f}, max_X = {max_X:.4f}, U*max_X = {U_old * max_X:.4f}")
            # UPDATE print(f"{U_it}) U = {self.U:.4f}, U_old*X = {np.max(np.abs(U_old * self.X.obj_wk.data)):.4f}")

            # Perform one FLEX loop iteration with reduced U (matching test.py logic)
            G_old = FLEX.G.obj_wk.copy()

            # Calculate V and Sigma with current chi0 and reduced U
            if U_it == 1:
                self.solve_FLEX_DMFT(FLEX, IPT)
            self.solve_FLEX_DMFT(FLEX, IPT)

            # Reset U back to U_old for next iteration
            diff = abs(prev_U - self.U)
            # UPDATE diff = np.max(np.abs(prev_U - self.U))
            prev_U = self.U
            self.U = U_old

            if U_it >= self.U_maxiter or diff < 1e-3:
                print(f"U_diff = {diff:.4e}, tol = 1e-3")
                print(f"Iteration number {U_it}, max iterations {self.U_maxiter}")
                break

        print("Leaving U renormalization...")
        self.U = prev_U
        print("Final U after renormalization: ", self.U)
        # Final UX calculation with U_old
        # UPDATE FLEX.UX = np.max(np.abs(self.U * self.X.obj_wk.data))
        FLEX.UX = self.U * np.max(np.abs(FLEX.X.obj_wk.data))
        FLEX.U = self.U
        IPT.U = self.U
        if FLEX.UX < 1.0:
            self.diverged = False
            print("U renormalization successful.")
            print("Final U*max(Chi) = ", FLEX.UX)

        if FLEX.UX >= 1.0:
            raise RuntimeError("U renormalization failed to reduce U*max(Chi) below 1.")

    def loop_FLEX_DMFT(self, n_loops=50, check_divergence=True):
        self.FLEX.U = self.U
        self.IPT.U = self.U
        self.solve_FLEX_DMFT(self.FLEX, self.IPT)
        print(f"Initial U*max(Chi) = {self.FLEX.UX:.4f}")
        # Initial check for divergence
        if check_divergence and self.diverged:
            print(f"Initial U*max(Chi) = {self.FLEX.UX:.4f} >= 1")
            self.U_renormalization(self.FLEX, self.IPT)

        print("Beginning FLEX+DMFT Self-Consistent Loop")
        for i in range(n_loops):
            #print(f"Starting iteration {i+1}/{n_loops}...")
            G_old = self.FLEX.G.obj_wk.copy()

            self.solve_FLEX_DMFT(self.FLEX, self.IPT)
            #print("FLEX+DMFT iteration completed.")
            if self.diverged:
                print(f"Divergence detected at iteration {i+1}")
                print(f"U*max(Chi) = {self.FLEX.UX:.4f} >= 1")
                print(f"Max(Chi) = {np.max(np.abs(self.FLEX.X.obj_wk.data)):.4f}")
                print(f"U = {self.U:.4f}")
                if check_divergence:
                    print("Attempting U renormalization...")
                    self.U_renormalization(self.FLEX, self.IPT)
                    continue
                else:
                    print("Code exiting due to divergence.")
                    break

            print(f"Iteration {i+1} completed. U*max(Chi) = {self.FLEX.UX:.4f}", end=' ')
            err = np.max(np.abs(self.FLEX.G.obj_wk.data - G_old.data))
            print(f"Max change in G: {err:.3e}")

            if err < 1e-6:
                print("Convergence achieved.")
                break



def add_local_to_nonlocal(A, B):
    temp = A.copy()
    temp.obj_wk.data[:] += B.obj_w.data[:, np.newaxis, :, :]
    return temp

def subtract_local_from_nonlocal(A, B):
    temp = A.copy()
    temp.obj_wk.data[:] -= B.obj_w.data[:, np.newaxis, :, :]
    return temp
