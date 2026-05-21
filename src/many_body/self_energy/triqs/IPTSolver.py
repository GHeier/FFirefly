from firefly.diagram import Diagram, dot_tr, dot_t
from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime
from triqs.gf import Gf, inverse, iOmega_n
import numpy as np

class IPTSolver:
    def __init__(self, beta, H, mu, n_loops=100, mix=0.10, tol=1e-6, w_max=1.2*6, eps=1e-14):
        self.beta = beta
        self.H = H
        self.max_loops = n_loops
        self.mix = mix
        self.tol = tol
        self.mu = mu 

        # Matsubara frequency Green's functions
        dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
        ws = np.array([float(iw.imag) for iw in dlr_iw_mesh], dtype=np.float32)
        print("Max w: ", max(ws))
        print("Min w: ", min(ws))

        G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
        self.G_weiss = Diagram(G_iw, 'Fermion')
        self.G_weiss_old = self.G_weiss.copy()
        self.Sigma_imp = self.G_weiss.copy()
        self.Sigma_imp.zero()

        if H is not None:
            # mu must be a matrix for matrix-valued Green's functions
            mu_matrix = self.mu * np.eye(1)
            self.G_weiss.obj_w << H(Sigma = self.Sigma_imp.obj_w, mu=mu_matrix)
        self.G_loc = self.G_weiss.copy()

    def get_IPT_Sigma(self, U, G_weiss):
        """Compute IPT self-energy from G_weiss without modifying self.Sigma_imp."""
        G_weiss.w_to_t()
        Sigma = G_weiss.copy()
        Sigma.obj_t << (U**2) * G_weiss.obj_t * G_weiss.obj_t * G_weiss.obj_t
        Sigma.t_to_w()
        return Sigma

    def set_Weiss(self):
        self.G_weiss_old = self.G_weiss.copy()
        self.G_weiss.obj_w << inverse(inverse(self.G_loc.obj_w) + self.Sigma_imp.obj_w)
        self.G_weiss.obj_w << self.mix * self.G_weiss_old.obj_w + (1.0 - self.mix) * self.G_weiss.obj_w

    def solve(self, U):
        self.Sigma_imp = self.get_IPT_Sigma(U, self.G_weiss)

        # Dyson
        mu_matrix = self.mu * np.eye(1)
        self.G_loc.obj_w << self.H(Sigma=self.Sigma_imp.obj_w, mu=mu_matrix)
        self.set_Weiss()

    def solve_bethe_lattice(self, U):
        self.Sigma_imp = self.get_IPT_Sigma(U, self.G_weiss)

        self.G_loc.obj_w = inverse(inverse(self.G_weiss.obj_w) - self.Sigma_imp.obj_w)
        t = 1
        self.G_weiss.obj_w << inverse( iOmega_n - t**2 * self.G_loc.obj_w )

    def loop(self, U, bethe_lattice=False):
        err_history = []
        for i in range(self.max_loops):
            G_iw_prev = self.G_weiss.obj_w.data.copy()  # Save for averaging

            if bethe_lattice:
                self.solve_bethe_lattice(U)
            else:
                self.solve(U)

            err = abs(self.G_weiss.obj_w.data - G_iw_prev).max()
            print("IPT loop %d, err = %.3e" % (i+1, err))

            if err < self.tol:
                break

            # Detect two-cycle oscillation
            #err_history.append(err)
            #if len(err_history) > 5:
            #    err_history.pop(0)
            #    err_std = np.std(err_history)
            #    if err_std < 1e-6:  # Error is constant
            #        print("Detected two-cycle oscillation. Averaging states.")
            #        # Average current and previous states
            #        self.G_weiss.obj_w.data[:] = 0.5 * (self.G_weiss.obj_w.data + G_iw_prev)
            #        self.G_weiss.w_to_t()
            #        err_history = []  # Reset history after intervention

