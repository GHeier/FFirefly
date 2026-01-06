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
        G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
        self.G_weiss = Diagram(G_iw, 'Fermion')
        self.Sigma_loc = self.G_weiss.copy()
        self.Sigma_loc.zero()

        if H is not None:
            self.G_weiss.obj_w << H(Sigma = self.Sigma_loc.obj_w, mu=self.mu)
        self.G_loc = self.G_weiss.copy()

    def get_IPT_Sigma(self, U):
        self.G_weiss.w_to_t()
        self.Sigma_loc.obj_t << (U**2) * self.G_weiss.obj_t * self.G_weiss.obj_t * self.G_weiss.obj_t
        self.Sigma_loc.t_to_w()
        return self.Sigma_loc

    def set_Weiss(self):
        self.G_weiss.obj_w << inverse(inverse(self.G_loc.obj_w) + self.Sigma_loc.obj_w)

    def solve(self, U):
        Sigma_iw = self.get_IPT_Sigma(U)
        self.Sigma_loc.obj_w = self.mix * Sigma_iw.obj_w + (1.0 - self.mix) * self.Sigma_loc.obj_w

        # Dyson
        #self.G << inverse(inverse(self.G0) - self.Sigma_iw)
        self.G_loc.obj_w << self.H(Sigma=self.Sigma_loc.obj_w, mu=self.mu)
        #self.G0 << inverse( iOmega_n - t**2 * self.G )
        self.set_Weiss()
        #self.G_iw = self.G0_iw * self.mix + self.G_iw * (1.0 - self.mix)

    def solve_bethe_lattice(self, U):
        self.G_weiss.w_to_t()
        self.Sigma_loc.obj_t << (U**2) * self.G_weiss.obj_t * self.G_weiss.obj_t * self.G_weiss.obj_t
        self.Sigma_loc.t_to_w()

        self.G_loc.obj_w = inverse(inverse(self.G_weiss.obj_w) - self.Sigma_loc.obj_w)
        t = 1
        self.G_weiss.obj_w << inverse( iOmega_n - t**2 * self.G_loc.obj_w )
        #self.G_iw = self.G0_iw * self.mix + self.G_iw * (1.0 - self.mix)

    def loop(self, U, bethe_lattice=False):
        for i in range(self.max_loops):
            G_iw_old = self.G_loc.obj_w.copy()
            if bethe_lattice:
                self.solve_bethe_lattice(U)
            else:
                self.solve(U)
            err = abs(self.G_loc.obj_w.data - G_iw_old.data).max()
            print("IPT loop %d, err = %.3e" % (i+1, err))
            if err < self.tol:
                break

