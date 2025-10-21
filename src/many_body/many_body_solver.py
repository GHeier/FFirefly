from diagram import *
from triqs_tprf.lattice import *
from triqs_tprf import *
from triqs.gf import inverse, Fourier
import numpy as np
from scipy.optimize import brentq

import firefly as fly
import firefly.config as cfg

class ManyBodySolver:
    def __init__(self, G0, U=0.0, mix=0.2, U_maxiter=50, e_k=None):
        self.G0 = Diagram(G0, 'Fermion')
        self.G = copy(self.G0)
        self.U = U
        self.UX = 0.0  # Track U * max(chi)
        self.mix = mix  # Mixing parameter for self-energy
        self.U_maxiter = U_maxiter  # Max iterations for U renormalization
        self.diverged = False  # Track divergence state
        self.e_k = e_k  # Energy dispersion (needed for mu calculation)
        self.n_target = None  # Target electron density
        self.mu = None  # Chemical potential

        # Store Sigma for mu calculation
        if not hasattr(self, 'Sigma'):
            self.Sigma = None

    def chi0_from_grt_PH(self):
        self.G.wk_to_tr()
        chi_tr = chi0_tr_from_grt_PH(self.G.obj_tr)
        self.X = Diagram(chi_tr, "Boson")
        self.X.tr_to_wk()

    def FLEX_from_chi(self):
        # FLEX vertex construction from spin and charge fluctuations
        # V = 3/2 U^2 chi_spin + 1/2 U^2 chi_charge - U^2 chi_0 + U
        U = self.U
        X_data = self.X.obj_wk.data

        # Check for divergence: U * max(chi) should be < 1
        self.UX = U * np.max(np.abs(X_data))
        if self.UX >= 1.0:
            print(f"ERROR: U*max(chi0) = {self.UX:.4f} >= 1! Paramagnetic phase reached - calculations unstable!")
            self.diverged = True
            return

        UX = U * X_data
        U2X = U**2 * X_data
        chi_spin = X_data / (1 - UX)
        chi_charge = X_data / (1 + UX)
        V_wk = 1.5 * U2X * chi_spin + 0.5 * U2X * chi_charge - U2X + U
        V = self.X.obj_wk.copy()
        V.data[:] = V_wk
        self.V = Diagram(V, 'Boson')
        self.diverged = False 

    def Sigma_from_vertex(self):
        self.V.wk_to_tr()
        self.Sigma = dot(self.G, self.V)
        self.Sigma.tr_to_wk()

    def dyson_G_from_Sigma(self):
        G_old = self.G.obj_wk.copy()
        G_new = inverse(inverse(self.G0.obj_wk) - self.Sigma.obj_wk)
        # Mix old and new Green's function for stability
        self.G.obj_wk.data[:] = self.mix * G_new.data + (1 - self.mix) * G_old.data

    def calc_electron_density(self, mu):
        """
        Calculate electron density from Green's function for a given chemical potential.

        Args:
            mu: Chemical potential

        Returns:
            n: Total electron density (including spin degeneracy factor of 2)
        """
        from triqs.gf import Gf, MeshDLRImTime

        # Get G(k,w) with updated mu
        G_wk = self.G0.obj_wk.copy()

        # Manually construct G: G(k,iw) = 1/(iw + mu - e_k - Sigma)
        # The mesh frequencies are already in G_wk.mesh
        for idx_w in range(len(G_wk.mesh.components[0])):
            for idx_k in range(len(G_wk.mesh.components[1])):
                iw = G_wk.mesh.components[0][idx_w]
                k_idx = idx_k

                if self.Sigma is not None:
                    # Interacting case
                    sigma_val = self.Sigma.obj_wk[idx_w, idx_k]
                    ek_val = self.e_k[idx_k]
                    G_wk.data[idx_w, idx_k] = 1.0 / (iw + mu - ek_val - sigma_val)
                else:
                    # Non-interacting case
                    ek_val = self.e_k[idx_k]
                    G_wk.data[idx_w, idx_k] = 1.0 / (iw + mu - ek_val)

        # Sum over k to get local Green's function G(iw)
        nk = len(G_wk.mesh.components[1])
        G_w_data = np.sum(G_wk.data, axis=1) / nk

        # Create a 1D Green's function on the frequency mesh only
        mesh_w = G_wk.mesh.components[0]
        G_w = Gf(mesh=mesh_w, target_shape=[])
        G_w.data[:] = G_w_data

        # Transform to imaginary time
        G_tau = make_gf_dlr_imtime(G_w)

        # Density from G(tau=0-): n = 1 + G(tau=0-) for each spin
        # With DLR, evaluate at tau=0-
        n_per_spin = 1.0 + np.real(G_tau.data[-1, 0, 0])

        # Factor of 2 for spin degeneracy
        n_total = 2.0 * n_per_spin

        return n_total

    def find_mu_for_density(self, n_target):
        """
        Find chemical potential that gives the target electron density using Brent's method.

        Args:
            n_target: Target electron density

        Returns:
            mu: Chemical potential that achieves n_target
        """
        # Get energy range from dispersion
        e_min = np.min(self.e_k.data.real)
        e_max = np.max(self.e_k.data.real)

        # Define function to find root of
        def density_error(mu):
            n = self.calc_electron_density(mu)
            error = n - n_target
            if cfg.verbosity == "high":
                print(f"  mu = {mu:.4f}, n = {n:.4f}, target = {n_target:.4f}, error = {error:.4f}")
            return error

        # Search for mu in expanded energy range
        mu_min = 3 * e_min
        mu_max = 3 * e_max

        print(f"Finding mu for n = {n_target:.4f}...")
        print(f"Energy range: [{e_min:.4f}, {e_max:.4f}]")
        print(f"Search range: [{mu_min:.4f}, {mu_max:.4f}]")

        try:
            mu = brentq(density_error, mu_min, mu_max, xtol=1e-6)
            print(f"Found mu = {mu:.4f} for n = {n_target:.4f}")
            self.mu = mu
            return mu
        except ValueError as e:
            print(f"ERROR: Could not find mu for n = {n_target}")
            print(f"Check that the target density is achievable in the range [{mu_min:.4f}, {mu_max:.4f}]")
            raise

    def solve(self):
        self.chi0_from_grt_PH()
        self.FLEX_from_chi()
        if self.diverged:
            return
        self.Sigma_from_vertex()
        self.dyson_G_from_Sigma()

    def U_renormalization(self):
        """Renormalize U if U*max(chi) >= 1 to avoid divergence."""
        print("WARNING: U is too large and the spin susceptibility will diverge!")
        print("Initiating U renormalization loop...")

        U_old = self.U
        U_it = 0
        U_diff = 1.0
        prev_U = 0.0

        while self.UX >= 1.0:
            U_it += 1
            # Reduce U to bring UX below 1
            self.U = U_old / (self.UX + 0.01)

            print(f"  U renorm iter {U_it}: U = {self.U:.4f}, UX = {self.UX:.4f}")

            # Perform one FLEX iteration with new U
            self.solve()
            if self.diverged:
                # If still diverging, continue reducing
                print(f"  Still diverging, reducing U further...")
                continue

            # Update UX with rescaled U
            self.UX *= U_old / self.U

            # Check convergence
            U_diff = abs(self.U - prev_U)
            prev_U = self.U

            if U_it >= self.U_maxiter or U_diff < 1e-4:
                print(f"U renormalization finished: iterations={U_it}, U_diff={U_diff:.3e}")
                break

        # Recompute UX with final U
        self.UX = self.U * np.max(np.abs(self.X.obj_wk.data))
        print(f"New U = {self.U:.4f}, UX = {self.UX:.4f}")

        # Check if U was reduced too much
        if self.U / U_old < 0.9:
            print("-----------------------------------------------")
            print("WARNING: U reduced by >10%! Paramagnetic phase is unavoidable!")
            print("-----------------------------------------------")
            exit()

    def loop(self, n_loops=50, check_divergence=True):
        # Initial check for divergence
        if check_divergence and self.UX >= 1.0:
            print(f"Initial U*max(Chi) = {self.UX:.4f} >= 1")
            self.U_renormalization()

        print("Beginning Self-Consistent Loop")
        for i in range(n_loops):
            print(f"Starting iteration {i+1}/{n_loops}...")
            G_old = self.G.obj_wk.copy()

            self.solve()
            if self.diverged:
                print(f"Divergence detected at iteration {i+1}")
                if check_divergence:
                    print("Attempting U renormalization...")
                    self.U_renormalization()
                    continue
                else:
                    print("Code exiting due to divergence.")
                    break

            print(f"Iteration {i+1} completed. U*max(Chi) = {self.UX:.4f}")
            err = np.max(np.abs(self.G.obj_wk.data - G_old.data))
            print(f"Max change in G: {err:.3e}")

            if err < 1e-6:
                print("Convergence achieved.")
                break


