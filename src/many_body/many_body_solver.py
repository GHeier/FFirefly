from diagram import *
from triqs_tprf.lattice import *
from triqs_tprf import *
from triqs.gf import inverse
import numpy as np

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

    def calculate_electron_density_from_G_wk(self):
        G_tk = fourier_wk_to_tk(self.G.obj_wk)
        n_k = np.real(np.trace(G_tk))
        n_total = np.mean([n_k for k in k_grid])    
        return n_total

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


