from diagram import *
from triqs_tprf.lattice import *
from triqs_tprf import *
from triqs.gf import inverse, Fourier, Gf, make_gf_dlr, make_gf_dlr_imtime, make_gf_dlr_imfreq
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime
import numpy as np
from scipy.optimize import brentq

import firefly as fly
import firefly.config as cfg

class ManyBodySolver:
    def __init__(self, G0=None, U=0.0, mix=0.2, U_maxiter=50, e_k=None, H=None, beta=None, w_max=None, eps=None):
        """
        Initialize ManyBodySolver for either FLEX or DMFT calculations.

        For FLEX: provide G0 (Green's function on k-mesh)
        For DMFT: provide H (HilbertTransform), beta, w_max, eps
        """
        self.U = U
        self.UX = 0.0  # Track U * max(chi)
        self.mix = mix  # Mixing parameter for self-energy
        self.U_maxiter = U_maxiter  # Max iterations for U renormalization
        self.diverged = False  # Track divergence state
        self.e_k = e_k  # Energy dispersion (needed for mu calculation)
        self.n_target = None  # Target electron density
        self.mu = None  # Chemical potential

        # FLEX mode
        if G0 is not None:
            self.G0 = Diagram(G0, 'Fermion')
            self.G = copy(self.G0)
            self.Sigma = None
            self.mode = 'FLEX'

        # DMFT mode
        elif H is not None:
            self.H = H
            self.beta = beta
            self.mode = 'DMFT'

            # Create DLR meshes
            dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)

            # Initialize Green's functions in frequency space
            G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
            Sigma_iw = G_iw.copy()
            Sigma_iw.zero()
            G_iw << H(Sigma=Sigma_iw, mu=0.0)

            # Create Diagram objects (they handle w <-> t transforms internally)
            self.G_loc = Diagram(G_iw, 'Fermion')        # Local Green's function
            self.Sigma_loc = Diagram(Sigma_iw, 'Fermion') # Local self-energy

        else:
            raise ValueError("Must provide either G0 (for FLEX) or H (for DMFT)")

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
        pass

    def find_mu_for_density(self, n_target):
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

    def solve_FLEX(self):
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
            self.solve_FLEX()
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

    def loop_FLEX(self, n_loops=50, check_divergence=True):
        # Initial check for divergence
        if check_divergence and self.UX >= 1.0:
            print(f"Initial U*max(Chi) = {self.UX:.4f} >= 1")
            self.U_renormalization()

        print("Beginning Self-Consistent Loop")
        for i in range(n_loops):
            print(f"Starting iteration {i+1}/{n_loops}...")
            G_old = self.G.obj_wk.copy()

            self.solve_FLEX()
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

    # DMFT methods
    def solve_DMFT(self):
        """Single DMFT iteration using IPT."""
        # Transform Weiss field to imaginary time
        self.G_loc.w_to_t()

        # IPT: Sigma(tau) = U^2 * G_Weiss(tau)^3
        self.Sigma_loc.obj_t.data[:] = (self.U**2) * self.G_loc.obj_t.data**3

        # Transform Sigma back to frequency
        self.Sigma_loc.t_to_w()

        # Mix self-energy for stability
        Sigma_old = self.Sigma_loc.obj_w.copy()
        self.Sigma_loc.obj_w.data[:] = self.mix * self.Sigma_loc.obj_w.data + (1.0 - self.mix) * Sigma_old.data

        # Dyson equation: G_loc = H(Sigma)
        self.G_loc.obj_w << self.H(Sigma=self.Sigma_loc.obj_w, mu=0.0)

        # Self-consistency: G_Weiss^-1 = G_loc^-1 + Sigma
        self.G_loc.obj_w << inverse(inverse(self.G_loc.obj_w) + self.Sigma_loc.obj_w)

    def loop_DMFT(self, n_loops=100, tol=1e-6):
        """Self-consistent DMFT loop."""
        print("Beginning DMFT Self-Consistent Loop")
        for i in range(n_loops):
            G_old = self.G_loc.obj_w.copy()

            self.solve_DMFT()

            err = np.max(np.abs(self.G_loc.obj_w.data - G_old.data))
            print(f"DMFT loop {i+1}, err = {err:.3e}")

            if err < tol:
                print("Convergence achieved.")
                break

    def solve_FLEX_DMFT(self):
        """Single iteration of FLEX+DMFT."""
        self.solve_DMFT()

        # FLEX steps
        self.solve_FLEX()

        if self.diverged:
            print("Divergence detected in FLEX+DMFT step.")
            return

        self.G_loc.obj_w.data[:] = np.sum(self.G.obj_wk.data, axis=1) / self.G.nk

        # Update local quantities
        self.G_loc = copy(self.G)
        self.Sigma_loc = copy(self.Sigma)
