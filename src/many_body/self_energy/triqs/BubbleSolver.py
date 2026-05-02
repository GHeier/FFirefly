from firefly.diagram import Diagram, contract, get_renorm, make_local
from triqs.gf import inverse, Gf
from triqs.gf.meshes import MeshDLRImFreq
from triqs_tprf.lattice import chi0_tr_from_grt_PH
from scipy.optimize import brentq, root
import numpy as np

class BubbleSolver:
    def __init__(self, G0=None, U=0.0, mix=0.2, n=None, mu=None):
        self.U = U
        print("U = ", U)
        self.mix = mix
        self.n = n
        self.mu = mu if mu is not None else 0.0

        if G0 is None:
            raise ValueError("Must provide G0 triqs object")

        self.G0 = Diagram(G0, 'Fermion')
        self.G = self.G0.copy()
        self.Sigma = self.G0.copy()
        self.Sigma.zero()

        self.emax = 0.0
        self.emin = 0.0

        if self.G0.varspace != 'wk':
            raise ValueError("G0 must be in 'wk' space")

        stuff = G0.mesh.components[0]
        self.beta = stuff.beta
        self.eps = stuff.eps
        self.w_max = stuff.w_max

        self.ind = self.init_iw_ek()

        dlr_iw_mesh = MeshDLRImFreq(beta=self.beta, statistic='Fermion', w_max=self.w_max, eps=self.eps)
        G_iw = Gf(mesh=dlr_iw_mesh, target_shape=G0.target_shape)
        self.G_loc = Diagram(G_iw, 'Fermion')

        #n = self.calc_electron_density(0.0)
        #print("Initial n for mu=0 is ", n)

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
        self.emax = np.max(G0_inv.data[0, :, :, :].real)
        self.emin = np.min(G0_inv.data[0, :, :, :].real)

        signs = np.sign(self.iw_arr.imag)
        diff = np.diff(signs)
        zero_crossings = np.where(diff != 0)[0]
        ind = zero_crossings[0]
        print(f"ind = {ind}")
        return ind

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
        self.V.wk_to_tr()
        self.Sigma.obj_tr = contract(self.V.obj_tr, self.G.obj_tr)
        self.Sigma.tr_to_wk()

    def get_local_G(self):
        self.G_loc.obj_w.data[:] = np.einsum('wknm->wnm', self.G.obj_wk.data) / self.G.nk

    def calc_electron_density(self, mu):
        self.make_G(mu)
        self.get_local_G()
        # Factor of 2 for spin degeneracy: TRIQS returns density per spin,
        # but num_electrons in config is total density
        return 2 * self.G_loc.obj_w.density().real[0][0]

    def find_mu_for_density(self, n_target):
        # Compute actual band extremes including self-energy at lowest Matsubara
        H_plus_Sigma = (self.H_k + self.Sigma.obj_wk.data[self.ind, :, :, :]).real
        margin = 0.5  # eV margin beyond band edges
        emax = np.max(H_plus_Sigma) + margin
        emin = np.min(H_plus_Sigma) - margin
        #print("Emax = ", emax)
        #print("Emin = ", emin)

        self.get_local_G()
        n = self.calc_electron_density(self.mu)
        if abs(n - n_target) < 1e-4:
            return self.mu

        def density_error(mu):
            return self.calc_electron_density(mu) - n_target

        bound = 2 if abs(n - n_target) < 1e-2 else 10
        mu_min = max(self.mu - bound, emin)
        mu_max = min(self.mu + bound, emax)

        try:
            #self.mu = brentq(density_error, mu_min, mu_max, xtol=1e-4)
            sol = root(density_error, self.mu, method='broyden1')
            #print("sol: ", sol.x)
            self.mu = sol.x
        except ValueError:
            # Check achievable density range
            n_at_min = self.calc_electron_density(emin)
            n_at_max = self.calc_electron_density(emax)
            err_at_min = n_at_min - n_target
            err_at_max = n_at_max - n_target

            # Check if bounds bracket the root (opposite signs)
            if err_at_min * err_at_max < 0:
                self.mu = brentq(density_error, emin, emax, xtol=1e-4)
            elif abs(err_at_min) < abs(err_at_max):
                print(f"Warning: n_target={n_target:.4f}, n_at_emin={n_at_min:.4f}, n_at_emax={n_at_max:.4f}, using emin")
                self.mu = emin
            else:
                print(f"Warning: n_target={n_target:.4f}, n_at_emin={n_at_min:.4f}, n_at_emax={n_at_max:.4f}, using emax")
                self.mu = emax
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

        old_m_star = 0.0
        for i in range(n_loops):
            G_old = self.G.obj_wk.copy()
            self.solve_Bubble()

            err = np.max(np.abs(self.G.obj_wk.data - G_old.data))
            sigma = make_local(self.Sigma.obj_wk.data)
            m_star = get_renorm(sigma, self.Sigma.w_points)
            print(f"Iteration {i+1}, m*/m = {m_star}, Max change in G: {err:.3e}")

            if err < 1e-4 or abs(old_m_star - m_star) < 1e-5:
                print("Convergence achieved.")
                break
            if i == n_loops - 1:
                print("Warning: Convergence not achieved")
            old_m_star = m_star

