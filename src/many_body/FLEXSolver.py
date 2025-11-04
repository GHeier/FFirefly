from firefly.diagram import Diagram, dot_tr, dot_t
from triqs.gf import inverse, Gf
from triqs.gf.meshes import MeshDLRImFreq
from triqs_tprf.lattice import chi0_tr_from_grt_PH
from scipy.optimize import brentq
import numpy as np

class FLEXSolver:
    def __init__(self, G0=None, U=0.0, mix=0.2, U_maxiter=50, n=None, mu=None):
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
        else:
            raise ValueError("G0 must be in 'wk' space")
        self.beta = stuff.beta
        self.eps = stuff.eps
        self.w_max = stuff.w_max

        self.init_iw_ek()

        dlr_iw_mesh = MeshDLRImFreq(beta=self.beta, statistic='Fermion', w_max=self.w_max, eps=self.eps)

        G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
        self.G_loc = Diagram(G_iw, 'Fermion')


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

    def chi0_from_grt_PH(self):
        self.G.wk_to_tr()
        chi_tr = chi0_tr_from_grt_PH(self.G.obj_tr)
        self.X = Diagram(chi_tr, "Boson")
        self.X.tr_to_wk()
        self.UX = self.U * np.max(np.abs(self.X.obj_wk.data))

    def FLEX_from_chi(self, X):
        # FLEX vertex construction from spin and charge fluctuations
        # V = 3/2 U^2 chi_spin + 1/2 U^2 chi_charge - U^2 chi_0 + U
        U = self.U
        X_data = X.data

        # Check for divergence: U * max(chi) should be < 1
        self.UX = U * np.max(np.abs(X_data))
        if self.UX >= 1.0:
            #print(f"ERROR: U*max(chi0) = {self.UX:.4f} >= 1! Paramagnetic phase reached - calculations unstable!")
            self.diverged = True

        UX = U * X_data
        chi_spin = X_data / (1 - UX)
        chi_charge = X_data / (1 + UX)

        # Store chi_spin and chi_charge for later use (e.g., singlet vertex)
        self.chi_spin_data = chi_spin
        self.chi_charge_data = chi_charge

        V_wk = 1.5 * U**2 * chi_spin + 0.5 * U**2 * chi_charge - U**2 * X_data #+ U
        #V_wk = U**2 * chi_spin + U**3 * chi_spin * chi_charge
        V = X.copy()
        V.data[:] = V_wk
        self.diverged = False
        return Diagram(V, 'Boson')

    #def compute_V_singlet(self):
    #    """
    #    Compute singlet pairing vertex for superconductivity.
    #    V_singlet = 3/2 U^2 chi_spin - 1/2 U^2 chi_charge
    #    (Note: opposite sign for charge channel compared to FLEX vertex)
    #    """
    #    if not hasattr(self, 'chi_spin_data') or not hasattr(self, 'chi_charge_data'):
    #        raise ValueError("chi_spin and chi_charge not computed yet. Run FLEX_from_chi first.")

    #    U = self.U
    #    V_singlet_wk = 1.5 * U**2 * self.chi_spin_data - 0.5 * U**2 * self.chi_charge_data

    #    V_singlet = self.X.copy()
    #    V_singlet.obj_wk.data[:] = V_singlet_wk
    #    return V_singlet

    def Sigma_from_vertex(self):
        self.V.wk_to_tr()
        self.Sigma = dot_tr(self.V, self.G)
        #self.Sigma.obj_tr.data[:] = self.G.obj_tr.data * self.V.obj_tr.data[:, :, 0, 0]
        self.Sigma.tr_to_wk()

    def dyson_G_from_Sigma(self):
        self.G.obj_wk = inverse(inverse(self.G0.obj_wk) - self.Sigma.obj_wk)

    def get_local_G(self):
        # Calculate local Green's function by summing over k-points
        self.G_loc.obj_w.data[:] = np.sum(self.G.obj_wk.data[:, :, :, :], axis=1) / self.G.nk

    def calc_electron_density(self, mu):
        """
        Calculate electron density at given chemical potential.
        Uses direct construction via make_G formula.
        """
        self.make_G(mu)

        # Sum over k-points to get local G at new mu
        self.G_loc.obj_w.data[:, 0, 0] = np.reshape(np.sum(self.G.obj_wk.data, axis=1) / self.G.nk, (self.G.nw))

        # Calculate density
        n = self.G_loc.obj_w.density().real[0][0]
        return n

    def find_mu_for_density(self, n_target):
        self.get_local_G()
        n = self.calc_electron_density(self.mu)
        if abs(n - n_target) < 1e-4:
            return self.mu
        def density_error(mu):
            n = self.calc_electron_density(mu)
            error = n - n_target
            return error

        # Search for mu in expanded energy range
        if abs(n - n_target) < 1e-2:
            mu_min = self.mu - 1.8
            mu_max = self.mu + 1.8
        else:
            e = inverse(self.G.obj_wk).data.real + self.mu
            mu_min = np.min(e)
            mu_max = np.max(e)

        try:
            mu = brentq(density_error, mu_min, mu_max, xtol=1e-4)
            self.mu = mu
            return mu
        except ValueError as e:
            print(f"ERROR: Could not find mu for n = {n_target}")
            print(f"Check that the target density is achievable in the range [{mu_min:.4f}, {mu_max:.4f}]")
            raise

    def shift_mu_to_target_density(self, n_target):
        mu = self.find_mu_for_density(n_target)
        self.make_G(mu)

    def solve_FLEX(self):
        # Use existing chi0 to calculate V and Sigma (matching test.py loop order)
        self.V = self.FLEX_from_chi(self.X.obj_wk)
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
            self.V = self.FLEX_from_chi(self.X.obj_wk)
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
                if check_divergence:
                    print("Attempting U renormalization...")
                    self.U_renormalization()
                    continue
                else:
                    print("Code exiting due to divergence.")
                    break

            err = np.max(np.abs(self.G.obj_wk.data - G_old.data))
            print(f"{i}) Max G(iw,k) diff = {err:.4e}, U*max(Chi) = {self.UX:.4f}")

            if err < 1e-6:
                print("Convergence achieved.")
                break
