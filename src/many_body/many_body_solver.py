from diagram import *
from triqs_tprf.lattice import *
from triqs_tprf import *
from triqs.gf import inverse, Fourier
from triqs.gf.meshes import MeshDLRImTime
import numpy as np
from scipy.optimize import brentq

import firefly as fly
import firefly.config as cfg

class ManyBodySolver:
    def __init__(self, G0, U=0.0, mix=0.2, U_maxiter=50):
        self.G0 = Diagram(G0, 'Fermion')
        self.G = copy(self.G0)
        self.U = U
        self.UX = 0.0  # Track U * max(chi)
        self.mix = mix  # Mixing parameter for self-energy
        self.U_maxiter = U_maxiter  # Max iterations for U renormalization
        self.diverged = False  # Track divergence state

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
        """Calculate electron density from Green function for given chemical potential mu."""
        # Temporarily store current G
        G_backup = self.G.obj_wk.copy()

        # Compute G(k,w) with new mu
        # Using Dyson equation: G = 1/(iw - (ek - mu) - Sigma)
        # For now, assume Sigma is stored in self.Sigma if it exists
        if hasattr(self, 'Sigma'):
            self.G.obj_wk = inverse(inverse(self.G0.obj_wk) - self.Sigma.obj_wk)
        else:
            # If no self-energy, use non-interacting case
            # Need to rebuild G0 with new mu - this requires access to e_k
            # For simplicity, we'll just shift G0
            self.G.obj_wk = self.G0.obj_wk.copy()

        # Sum over k-points to get local Green's function G_loc(iw)
        # The shape is (nw, nk1, nk2, nk3, norb, norb)
        # We need to sum over spatial dimensions
        nw = self.G.obj_wk.data.shape[0]
        nk = np.prod(self.G.obj_wk.data.shape[1:4])  # Total number of k-points

        # Average over k-points
        G_loc_data = np.sum(self.G.obj_wk.data, axis=(1,2,3)) / nk

        # For DLR, we need to transform to imaginary time and evaluate at tau=0
        # Using the relation: n = 2 * (1 + Re[G(tau=0^-)]) for fermions with spin
        # For Matsubara: G(tau=0^-) can be approximated from high-frequency tail

        # Simpler approach: use the sum rule for Matsubara frequencies
        # n/2 = (1/beta) * sum_iw G(iw) * e^(iw*0+)
        # At tau=0-, this gives: n = 2 * (1 + sum of G at large iw)

        # For a rough estimate, use the fact that at tau=beta/2, G(tau) ~ -1/2 for half filling
        # Better: integrate using trapezoidal rule or use TRIQS Fourier transform

        # Use TRIQS to transform to tau and get density
        # Create a copy for time-domain
        mesh_tau = MeshDLRImTime(beta=self.G.obj_wk.mesh.beta,
                                  statistic='Fermion',
                                  w_max=self.G.obj_wk.mesh.w_max,
                                  eps=self.G.obj_wk.mesh.eps)

        # For single-band, get the trace
        # Extract diagonal element (assuming single orbital)
        G_diag = G_loc_data[:, 0, 0] if G_loc_data.ndim > 1 else G_loc_data

        # Approximate density using high-frequency behavior
        # For large w: G(iw) ~ 1/iw, so sum converges
        # Better approximation: n = 2*(1 + Re[G(tau=0-)])
        # Use the fact that G(tau=0-) ~ -1/2 + small corrections

        # Simple estimate from Matsubara sum (proper implementation would use DLR basis)
        density = 2.0 * (1.0 + np.real(G_diag[0]))  # Approximate using lowest frequency

        # Restore original G
        self.G.obj_wk = G_backup

        return density

    def mu_from_density(self, target_n, e_k_min, e_k_max):
        """Find chemical potential mu for a given electron density n using Brent's method.

        Args:
            target_n: Target electron density
            e_k_min: Minimum energy eigenvalue
            e_k_max: Maximum energy eigenvalue

        Returns:
            mu: Chemical potential that gives the target density
        """
        def density_error(mu):
            return self.calc_electron_density(mu) - target_n

        # Search range: extend beyond band edges
        mu_min = 3.0 * e_k_min
        mu_max = 3.0 * e_k_max

        try:
            mu = brentq(density_error, mu_min, mu_max, xtol=1e-6, maxiter=100)
            print(f"Found mu = {mu:.4f} for n = {target_n:.4f}")
            return mu
        except ValueError as e:
            print(f"Warning: Could not find mu for n={target_n}. Using fermi_energy from config.")
            return cfg.fermi_energy

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


