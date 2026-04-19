from firefly.diagram import Diagram, contract
from triqs.gf import inverse, Gf
from triqs.gf.meshes import MeshDLRImFreq
from triqs_tprf.lattice import chi0_tr_from_grt_PH
from scipy.optimize import brentq
import numpy as np

class BubbleSolver:
    def __init__(self, G0=None, U=0.0, mix=0.2, n=None, mu=None):
        self.U = U
        self.mix = mix
        self.n = n
        self.mu = mu if mu is not None else 0.0

        if G0 is None:
            raise ValueError("Must provide G0 triqs object")

        self.G0 = Diagram(G0, 'Fermion')
        self.G = self.G0.copy()
        self.Sigma = self.G0.copy()
        self.Sigma.zero()

        if self.G0.varspace != 'wk':
            raise ValueError("G0 must be in 'wk' space")

        stuff = G0.mesh.components[0]
        self.beta = stuff.beta
        self.eps = stuff.eps
        self.w_max = stuff.w_max

        self.init_iw_ek()

        dlr_iw_mesh = MeshDLRImFreq(beta=self.beta, statistic='Fermion', w_max=self.w_max, eps=self.eps)
        G_iw = Gf(mesh=dlr_iw_mesh, target_shape=G0.target_shape)
        self.G_loc = Diagram(G_iw, 'Fermion')

        self.mu = self.find_mu_for_density(n)
        self.make_G(self.mu)
        self.chi0_from_grt_PH()

    def init_iw_ek(self):
        mesh_w = self.G.obj_wk.mesh.components[0]
        self.iw_arr = np.array([complex(iw) for iw in mesh_w], dtype=np.complex128)
        G0_inv = inverse(self.G0.obj_wk)
        # Extract H_k from G0^-1 = iw + mu - H_k, assuming G0 was created with some mu
        self.H_k = self.iw_arr[0] + self.mu - G0_inv.data[0, :, :, :]
        self.iw_broadcast = self.iw_arr[:, np.newaxis, np.newaxis, np.newaxis]

    def make_G(self, mu):
        self.G.obj_wk.data[:] = 1.0 / (self.iw_broadcast + mu - self.H_k - self.Sigma.obj_wk.data)

    def chi0_from_grt_PH(self):
        self.G.wk_to_tr()
        chi_tr = chi0_tr_from_grt_PH(self.G.obj_tr)
        self.X = Diagram(chi_tr, "Boson")
        self.X.tr_to_wk()

    def Bubble_from_chi(self, X):
        """Bubble vertex: V = U^2 * chi0"""
        V = X.copy()
        V.data[:] = self.U**2 * X.data
        return Diagram(V, 'Boson')

    def Sigma_from_vertex(self):
        self.X.wk_to_tr()
        self.Sigma.obj_tr = contract(self.X.obj_tr, self.G.obj_tr)
        self.Sigma.tr_to_wk()

    def get_local_G(self):
        self.G_loc.obj_w.data[:] = np.einsum('wknm->wnm', self.G.obj_wk.data) / self.G.nk

    def calc_electron_density(self, mu):
        self.make_G(mu)
        self.get_local_G()
        return self.G_loc.obj_w.density().real[0][0]

    def find_mu_for_density(self, n_target):
        self.get_local_G()
        n = self.calc_electron_density(self.mu)
        if abs(n - n_target) < 1e-4:
            return self.mu

        def density_error(mu):
            return self.calc_electron_density(mu) - n_target

        e = inverse(self.G.obj_wk).data.real + self.mu
        emin, emax = np.min(e), np.max(e)
        bound = 2 if abs(n - n_target) < 1e-2 else 10
        mu_min = max(self.mu - bound, emin)
        mu_max = min(self.mu + bound, emax)

        try:
            self.mu = brentq(density_error, mu_min, mu_max, xtol=1e-4)
        except ValueError:
            # Check achievable density range
            n_at_min = self.calc_electron_density(emin)
            n_at_max = self.calc_electron_density(emax)
            if n_target < n_at_min:
                print(f"Warning: n_target={n_target:.4f} < n_min={n_at_min:.4f}, using emin")
                self.mu = emin
            elif n_target > n_at_max:
                print(f"Warning: n_target={n_target:.4f} > n_max={n_at_max:.4f}, using emax")
                self.mu = emax
            else:
                self.mu = brentq(density_error, emin, emax, xtol=1e-4)
        return self.mu

    def shift_mu_to_target_density(self, n_target):
        mu = self.find_mu_for_density(n_target)
        self.make_G(mu)

    def solve_Bubble(self):
        self.V = self.Bubble_from_chi(self.X.obj_wk)
        self.Sigma_from_vertex()

        G_old = self.G.obj_wk.copy()
        self.shift_mu_to_target_density(self.n)
        self.G.obj_wk.data[:] = self.mix * self.G.obj_wk.data + (1 - self.mix) * G_old.data
        self.chi0_from_grt_PH()

    def loop_Bubble(self, n_loops=50):
        self.chi0_from_grt_PH()
        print("Beginning Bubble Self-Consistent Loop")

        for i in range(n_loops):
            G_old = self.G.obj_wk.copy()
            self.solve_Bubble()

            err = np.max(np.abs(self.G.obj_wk.data - G_old.data))
            print(f"{i}) Max G(iw,k) diff = {err:.4e}")

            if err < 1e-4:
                print("Convergence achieved.")
                break
