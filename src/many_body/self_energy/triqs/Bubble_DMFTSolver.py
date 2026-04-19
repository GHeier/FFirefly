from triqs.gf import inverse, Gf
from triqs.gf.meshes import MeshDLRImFreq
from IPTSolver import IPTSolver
from BubbleSolver import BubbleSolver
import numpy as np

from firefly.diagram import Diagram, dot_t

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

    def solve_Bubble_DMFT(self, Bubble, IPT):
        mu = Bubble.find_mu_for_density(self.n)
        Bubble.make_G(self.mu)
        Bubble.get_local_G()
        IPT.G_loc = Bubble.G_loc.copy()
        IPT.set_Weiss()

        sigma_ipt_old = self.Sigma_imp.obj_w.copy()
        self.Sigma_imp = IPT.get_IPT_Sigma(self.U)

        new_sigma_k = add_local_to_nonlocal(self.Sigma_nonloc, self.Sigma_imp)
        self.Sigma_k.obj_wk.data[:] = self.mix * new_sigma_k.obj_wk.data + (1.0 - self.mix) * self.Sigma_k.obj_wk.data
        Bubble.Sigma = self.Sigma_k.copy()

        mu = Bubble.find_mu_for_density(self.n)
        Bubble.make_G(mu)
        Bubble.get_local_G()

        Bubble.chi0_from_grt_PH()
        Bubble.V = Bubble.Bubble_from_chi(Bubble.X.obj_wk)
        Bubble.Sigma_from_vertex()

        self.Sigma_loc.obj_w.data[:] = make_local(Bubble.Sigma.obj_wk.data)
        self.Sigma_nonloc = subtract_local_from_nonlocal(Bubble.Sigma, self.Sigma_loc)

        new_sigma_k = add_local_to_nonlocal(self.Sigma_nonloc, self.Sigma_imp)
        self.Sigma_k.obj_wk.data[:] = self.mix * new_sigma_k.obj_wk.data + (1.0 - self.mix) * self.Sigma_k.obj_wk.data
        Bubble.Sigma_k = self.Sigma_k.copy()

    def loop_Bubble_DMFT(self, n_loops=50):
        self.Bubble.U = self.U
        self.IPT.U = self.U
        self.solve_Bubble_DMFT(self.Bubble, self.IPT)

        print("Beginning Bubble+DMFT Self-Consistent Loop")
        for i in range(n_loops):
            G_old = self.Bubble.G.obj_wk.copy()
            self.solve_Bubble_DMFT(self.Bubble, self.IPT)

            err = np.max(np.abs(self.Bubble.G.obj_wk.data - G_old.data))
            print(f"Iteration {i+1}, Max change in G: {err:.3e}")

            if err < 1e-4:
                print("Convergence achieved.")
                break


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
