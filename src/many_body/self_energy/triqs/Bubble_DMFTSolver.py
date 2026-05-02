from triqs.gf import inverse, Gf
from triqs.gf.meshes import MeshDLRImFreq
from IPTSolver import IPTSolver
from BubbleSolver import BubbleSolver
import numpy as np

from firefly.diagram import Diagram, dot_t, get_renorm

class Bubble_DMFTSolver:
    def __init__(self, G0, U=0.0, mix=0.2, n=None, mu=None):
        self.U = U
        self.mix = mix
        self.n = n
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

        self.Bubble = BubbleSolver(self.G0.obj_wk, U, mix=mix, n=n, mu=mu)
        mu = self.Bubble.find_mu_for_density(n)
        self.IPT = IPTSolver(self.beta, None, mix=mix, w_max=self.w_max, eps=self.eps, mu=mu)

        DLR_b = MeshDLRImFreq(beta=self.beta, statistic='Boson', w_max=self.w_max, eps=self.eps)

        self.Sigma_loc = self.IPT.G_loc.copy()
        self.Sigma_imp = self.Sigma_loc.copy()
        self.Sigma_nonloc = self.G0.copy()
        self.Sigma_nonloc.zero()
        self.Sigma_k = self.Sigma_nonloc.copy()

        self.X_loc = Diagram(Gf(mesh=DLR_b, target_shape=(self.G0.obj_wk.target_shape)), 'Boson')

    def solve_Bubble_DMFT_ave(self, Bubble, IPT):
        mu = Bubble.find_mu_for_density(self.n)
        Bubble.make_G(mu)
        Bubble.get_local_G()
        IPT.G_loc = Bubble.G_loc.copy()
        IPT.set_Weiss()

        self.Sigma_imp = IPT.get_IPT_Sigma(self.U, IPT.G_weiss)

        new_sigma_k = add_local_to_nonlocal(self.Sigma_nonloc, self.Sigma_imp)
        self.Sigma_k.obj_wk.data[:] = self.mix * new_sigma_k.obj_wk.data + (1.0 - self.mix) * self.Sigma_k.obj_wk.data
        Bubble.Sigma = self.Sigma_k.copy()

        mu = Bubble.find_mu_for_density(self.n)
        Bubble.make_G(mu)

        Bubble.chi0_from_grt_PH()
        Bubble.V = Bubble.Bubble_from_chi(Bubble.X.obj_wk)
        Bubble.Sigma_from_vertex()

        self.Sigma_loc.obj_w.data[:] = make_local(Bubble.Sigma.obj_wk.data)
        self.Sigma_nonloc = subtract_local_from_nonlocal(Bubble.Sigma, self.Sigma_loc)

        new_sigma_k = add_local_to_nonlocal(self.Sigma_nonloc, self.Sigma_imp)
        self.Sigma_k.obj_wk.data[:] = self.mix * new_sigma_k.obj_wk.data + (1.0 - self.mix) * self.Sigma_k.obj_wk.data
        Bubble.Sigma = self.Sigma_k.copy()

    def solve_Bubble_DMFT_diag(self, Bubble, IPT):
        # Step 1 - solve DMFT
        IPT.loop(self.U) 

        # Step 2 - construct G(k,iw)
        new_Sigma_k = add_local_to_nonlocal(self.Sigma_nonloc, IPT.Sigma_imp)
        self.Sigma_k.obj_wk.data[:] = self.mix * new_sigma_k.obj_wk.data + (1.0 - self.mix) * self.Sigma_k.obj_wk.data
        Bubble.Sigma = self.Sigma_k.copy()
        Bubble.make_G(self.mu)

        # Step 3 - Find Σ^(2)
        Bubble.chi0_from_grt_PH()
        Bubble.V = Bubble.Bubble_from_chi(Bubble.X.obj_wk)
        Bubble.Sigma_from_vertex()

        # Step 4 - Find Σ^(2)[G]
        Bubble.get_local_G()
        self.Sigma_loc = IPT.get_IPT_Sigma(self.U, Bubble.G_loc)
        self.Sigma_loc.obj_w.data[:] = IPT.Sigma_imp.obj_w.data - self.Sigma_loc.obj_w.data

        # Step 5 - Make Σ(k,iw)
        new_sigma_k = add_local_to_nonlocal(self.Sigma_nonloc, self.Sigma_loc)
        self.Sigma_k.obj_wk.data[:] = self.mix * new_sigma_k.obj_wk.data + (1.0 - self.mix) * self.Sigma_k.obj_wk.data
        Bubble.Sigma = self.Sigma_k.copy()

        # Step 6 - Update mu and G / G_weiss
        self.mu = Bubble.find_mu_for_density(self.n)
        Bubble.make_G(self.mu)
        Bubble.get_local_G()
        IPT.G_loc = Bubble.G_loc.copy()
        IPT.set_Weiss()


    def loop_Bubble_DMFT(self, n_loops=50, dmft_ave=True):
        self.Bubble.U = self.U
        self.IPT.U = self.U
        old_m_star = 0.0

        print("Beginning Bubble+DMFT Self-Consistent Loop")
        for i in range(n_loops):
            G_old = self.Bubble.G.obj_wk.copy()
            if dmft_ave:
                self.solve_Bubble_DMFT_ave(self.Bubble, self.IPT)
            else:
                self.solve_Bubble_DMFT_diag(self.Bubble, self.IPT)

            err = np.max(np.abs(self.Bubble.G.obj_wk.data - G_old.data))
            m_star = get_renorm(self.Sigma_imp.obj_w.data, self.Sigma_imp.w_points)
            print(f"Iteration {i+1}, m*/m = {m_star}, Max change in G: {err:.3e}")

            if err < 1e-4 or abs(old_m_star - m_star) < 1e-5:
                print("Convergence achieved.")
                break
            if i == n_loops - 1:
                print("Warning: Convergence not achieved")
            old_m_star = m_star


def add_local_to_nonlocal(A, B):
    temp = A.copy()
    temp.obj_wk.data[:] += B.obj_w.data[:, np.newaxis, :, :]
    return temp

def subtract_local_from_nonlocal(A, B):
    temp = A.copy()
    temp.obj_wk.data[:] -= B.obj_w.data[:, np.newaxis, :, :]
    return temp

def make_local(data):
    return np.einsum('wknm->wnm', data) / data.shape[1]
