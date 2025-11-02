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

class ManyBodySolver:
    def __init__(self, G0=None, U=0.0, mix=0.2, U_maxiter=50, H=None, n=None, mu=None):
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
        self.n = n   # Target electron density (per spin)
        self.mu = mu

        # FLEX mode
        if G0 is None:
            raise ValueError("Must provide G0 triqs object")

        self.G0 = Diagram(G0, 'Fermion')
        self.G = self.G0.copy()
        self.Sigma = self.G0.copy()
        self.Sigma.zero()
        if self.G0.varspace == 'wk':
            stuff = G0.mesh.components[0]
        elif self.G0.varspace == 'w':
            stuff = G0.mesh
        else:
            raise ValueError("G0 must be in 'wk' or 'w' space")
        self.beta = stuff.beta
        self.eps = stuff.eps
        self.w_max = stuff.w_max

        dlr_iw_mesh = MeshDLRImFreq(beta=self.beta, statistic='Fermion', w_max=self.w_max, eps=self.eps)

        G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
        Sigma_iw = G_iw.copy()
        Sigma_iw.zero()
        self.G_loc = Diagram(G_iw, 'Fermion')        # Local Green's function
        self.Sigma_loc = Diagram(Sigma_iw, 'Fermion') # Local self-energy
        self.H = H  # HilbertTransform for DMFT

        # Initialize iw and eps(k) arrays for Green's function construction
        self.init_iw_ek()

    def init_iw_ek(self):
        """Initialize iw and eps(k) arrays for Green's function construction."""
        mesh_w = self.G.obj_wk.mesh.components[0]
        self.iw_arr = np.array([complex(iw) for iw in mesh_w], dtype=np.complex128)

        # Extract band energies eps(k) from G0^-1(k,iw) = iw + mu_old - eps(k)
        G0_inv = inverse(self.G0.obj_wk)
        self.eps_k = self.iw_arr[0] + self.mu - G0_inv.data[0, :, 0, 0]  # eps(k) array, shape (nk,)

        # Reshape for broadcasting
        self.iw_broadcast = self.iw_arr[:, np.newaxis, np.newaxis, np.newaxis]
        self.eps_broadcast = self.eps_k[np.newaxis, :, np.newaxis, np.newaxis]


    def make_G(self, mu):
        # Compute G(k,iw) = 1 / (iw + mu - eps(k) - Sigma(k,iw))
        self.G.obj_wk.data[:] = 1.0 / (self.iw_broadcast + mu - self.eps_broadcast - self.Sigma.obj_wk.data)

    def get_local_G(self):
        """Calculate local Green's function by summing over k-points."""
        #G_loc= self.G_loc.obj_w.copy()
        #G_loc.data[:] = np.sum(self.G.obj_wk.data, axis=1) / self.G.nk
        #self.G_loc = Diagram(G_loc_wk, 'Fermion')
        self.G_loc.obj_w.data[:] = np.sum(self.G.obj_wk.data[:, :, :, :], axis=1) / self.G.nk

    def chi0_from_grt_PH(self):
        self.G.wk_to_tr()
        chi_tr = chi0_tr_from_grt_PH(self.G.obj_tr)
        self.X = Diagram(chi_tr, "Boson")
        self.X.tr_to_wk()
        self.UX = self.U * np.max(np.abs(self.X.obj_wk.data))

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
        U2X = U**2 * X_data * X_data
        chi_spin = X_data / (1 - UX)
        chi_charge = X_data / (1 + UX)
        V_wk = 1.5 * U**2 * chi_spin + 0.5 * U**2 * chi_charge - U**2 * X_data #+ U
        #V_wk = U**2 * chi_spin + U**3 * chi_spin * chi_charge
        V = self.X.obj_wk.copy()
        V.data[:] = V_wk
        self.V = Diagram(V, 'Boson')
        self.diverged = False 

    def Sigma_from_vertex(self):
        self.V.wk_to_tr()
        self.Sigma = dot_tr(self.V, self.G)
        #self.Sigma.obj_tr.data[:] = self.G.obj_tr.data * self.V.obj_tr.data[:, :, 0, 0]
        self.Sigma.tr_to_wk()

    def dyson_G_from_Sigma(self):
        G_old = self.G.obj_wk.copy()
        G_new = inverse(inverse(self.G0.obj_wk) - self.Sigma.obj_wk)
        # Mix old and new Green's function for stability
        self.G.obj_wk.data[:] = self.mix * G_new.data + (1 - self.mix) * G_old.data

    def calc_electron_density(self, mu):
        """
        Calculate electron density at given chemical potential.
        Uses direct construction via make_G formula.
        """
        self.make_G(mu)

        # Sum over k-points to get local G at new mu
        G_loc = self.G_loc.obj_w.copy()
        G_loc.data[:, 0, 0] = np.reshape(np.sum(self.G.obj_wk.data, axis=1) / self.G.nk, (self.G.nw))

        # Calculate density
        n = G_loc.density().real[0][0]
        return n

    def find_mu_for_density(self, n_target):
        self.get_local_G()
        n = self.calc_electron_density(self.mu)
        if abs(n - n_target) < 1e-4:
            #print(f"Current mu = {self.mu: .4f} gives n = {n:.4f}, close to target n = {n_target:.4f}")
            return self.mu
        def density_error(mu):
            n = self.calc_electron_density(mu)
            error = n - n_target
            #if cfg.verbosity == "high":
            #    print(f"  mu = {mu: .4f}, n = {n:.4f}, target = {n_target:.4f}, error = {error: .4f}")
            return error

        # Search for mu in expanded energy range
        if abs(n - n_target) < 1e-2:
            mu_min = self.mu - 1.8
            mu_max = self.mu + 1.8
        else:
            e = inverse(self.G.obj_wk).data.real + self.mu
            mu_min = np.min(e)
            mu_max = np.max(e)

        #print(f"Finding mu for n = {n_target:.4f}...")
        #print(f"Energy range: [{e_min:.4f}, {e_max:.4f}]")
        #print(f"Search range: [{mu_min:.4f}, {mu_max:.4f}]")

        try:
            mu = brentq(density_error, mu_min, mu_max, xtol=1e-4)
            #print(f"  Found mu = {mu:.4f} for n = {n_target:.4f}")
            self.mu = mu
            return mu
        except ValueError as e:
            print(f"ERROR: Could not find mu for n = {n_target}")
            print(f"Check that the target density is achievable in the range [{mu_min:.4f}, {mu_max:.4f}]")
            raise

    def shift_mu_to_target_density(self, n_target):
        mu = self.find_mu_for_density(n_target)
        #print(f"  Shifting mu to {mu:.4f} to achieve target density n = {n_target:.4f}")
        self.make_G(mu)

    def solve_FLEX(self):
        # Use existing chi0 to calculate V and Sigma (matching test.py loop order)
        self.FLEX_from_chi()
        if self.diverged:
            return
        self.Sigma_from_vertex()

        # Match test.py: find mu, calculate G, then mix (don't use dyson_G_from_Sigma)
        G_old = self.G.obj_wk.copy()
        self.shift_mu_to_target_density(self.n)  # Updates G at new mu
        # Apply mixing after G update
        self.G.obj_wk.data[:] = self.mix * self.G.obj_wk.data + (1 - self.mix) * G_old.data

        # Calculate new chi0 for next iteration
        self.chi0_from_grt_PH()

    def U_renormalization(self):
        """Renormalize U if U*max(chi) >= 1 to avoid divergence."""
        print("WARNING: U is too large and the spin susceptibility will diverge!")
        print("Initiating U renormalization loop...")

        U_old = self.U
        U_it = 0
        prev_U = 0

        # Check condition: U_old * max(chi0) >= 1
        while U_old * np.max(np.abs(self.X.obj_wk.data)) >= 1.0:
            U_it += 1

            # Reduce U temporarily to bring UX below 1
            max_X = np.max(np.abs(self.X.obj_wk.data))
            self.U = self.U / (max_X * self.U + 0.01)
            print(f"{U_it}) U = {self.U}, initial_U = {U_old}")

            # Perform one FLEX loop iteration with reduced U (matching test.py logic)
            G_old = self.G.obj_wk.copy()

            # Calculate V and Sigma with current chi0 and reduced U
            self.FLEX_from_chi()
            if not self.diverged:
                self.Sigma_from_vertex()

            # Update mu, then update G with mixing
            self.shift_mu_to_target_density(self.n)
            self.G.obj_wk.data[:] = self.mix * self.G.obj_wk.data + (1 - self.mix) * G_old.data

            # Recalculate chi0 from updated G
            self.chi0_from_grt_PH()

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

    def loop_FLEX(self, n_loops=50, check_divergence=True):
        self.chi0_from_grt_PH()
        # Initial check for divergence
        if check_divergence and self.UX >= 1.0:
            print(f"Initial U*max(Chi) = {self.UX:.4f} >= 1")
            self.U_renormalization()

        print("Beginning Self-Consistent Loop")
        for i in range(n_loops):
            G_old = self.G.obj_wk.copy()

            self.solve_FLEX()

            # Print max X for comparison with test.py (from newly calculated chi0)
            max_X = np.max(np.abs(self.X.obj_wk.data))
            # Recalculate UX with new chi0 for accurate reporting
            self.UX = self.U * max_X

            if self.diverged:
                print(f"Divergence detected at iteration {i+1}")
                print("Code exiting due to divergence.")
                break
                    #if check_divergence:
                    #    print("Attempting U renormalization...")
                    #    self.U_renormalization()
                    #    continue
                    #else:
                    #    print("Code exiting due to divergence.")
                    #    break

            err = np.max(np.abs(self.G.obj_wk.data - G_old.data))
            #print(f"max X = {max_X:.6f}")
            print(f"{i}) Max G(iw,k) diff = {err:.4e}, U*max(Chi) = {self.UX:.4f}")

            if err < 1e-6:
                print("Convergence achieved.")
                break

    # DMFT methods
    def get_IPT_sigma(self):
        """Single DMFT iteration using IPT."""
        # Transform Weiss field to imaginary time
        self.G_loc.w_to_t()

        # IPT: Sigma(tau) = U^2 * G_Weiss(tau)^3
        self.Sigma_loc.obj_t.data[:] = (self.U**2) * self.G_loc.obj_t.data**3

        Sigma_old = self.Sigma_loc.obj_w.copy()
        # Transform Sigma back to frequency
        self.Sigma_loc.t_to_w()
        # Mix self-energy for stability
        self.Sigma_loc.obj_w.data[:] = self.mix * self.Sigma_loc.obj_w.data + (1.0 - self.mix) * Sigma_old.data


    def solve_DMFT(self):
        self.get_IPT_sigma()
        # Dyson equation: G_loc = H(Sigma)
        self.G_loc.obj_w << self.H(Sigma=self.Sigma_loc.obj_w, mu=self.mu)

        # Self-consistency: G_Weiss^-1 = G_loc^-1 + Sigma
        self.G_loc.obj_w << inverse(inverse(self.G_loc.obj_w) + self.Sigma_loc.obj_w)
        self.shift_mu_to_target_density(self.n)

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
        # DMFT part
        self.dyson_G_from_Sigma()
        self.get_local_G()
        self.get_IPT_sigma()

        # Make full Sigma(k,iw) from local Sigma(iw)
        # Sigma_loc.obj_w.data has shape (nw, 1, 1), broadcast to (nw, nk, 1, 1)
        self.Sigma.obj_wk.data[:] += self.Sigma_loc.obj_w.data[:, np.newaxis, :, :]

        # FLEX part
        self.dyson_G_from_Sigma()
        self.shift_mu_to_target_density(self.n)
        self.chi0_from_grt_PH()
        self.FLEX_from_chi()
        if self.diverged:
            return
        self.Sigma_from_vertex()
        # Extract local part: sigma_loc has shape (nw, 1, 1)
        sigma_loc = np.sum(self.Sigma.obj_wk.data, axis=1) / self.Sigma.nk
        # Subtract non-local FLEX part and add back local DMFT part
        self.Sigma.obj_wk.data[:] -= sigma_loc[:, np.newaxis, :, :]
        self.Sigma.obj_wk.data[:] += self.Sigma_loc.obj_w.data[:, np.newaxis, :, :]

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
