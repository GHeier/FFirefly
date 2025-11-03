class IPTSolver:
    def __init__(self, beta, H, n_loops=100, mix=0.10, tol=1e-6, w_max=1.2*6, eps=1e-14):
        self.beta = beta
        self.H = H
        self.max_loops = n_loops
        self.mix = mix
        self.tol = tol
        self.mu = cfg.fermi_energy

        # Matsubara frequency Green's functions
        dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
        self.G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1,1])
        self.Sigma_iw = self.G_iw.copy()
        self.Sigma_iw.zero()
        self.G_iw << H(Sigma = self.Sigma_iw, mu=self.mu)
        self.G0_iw = self.G_iw.copy() # self.G0 will be set by the user after initialization

        # Imaginary time
        tau_mesh = MeshDLRImTime(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
        self.G_tau = Gf(mesh=tau_mesh, target_shape=[1,1])
        self.Sigma_tau = self.G_tau.copy()

    def get_IPT_Sigma(self, U):
        G_tau = iw_to_tau_dlr(self.G0_iw)
        Sigma_tau = G_tau * G_tau * G_tau * (U**2)
        Sigma_iw = tau_to_iw_dlr(Sigma_tau)
        return Sigma_iw

    def solve(self, U):
        Sigma_iw = self.get_IPT_Sigma(U)
        self.Sigma_iw = self.mix * Sigma_iw + (1.0 - self.mix) * self.Sigma_iw

        # Dyson
        #self.G_iw << inverse(inverse(self.G0_iw) - self.Sigma_iw)
        self.G_iw << self.H(Sigma=self.Sigma_iw, mu=self.mu)
        #self.G0_iw << inverse( iOmega_n - t**2 * self.G_iw )
        self.G0_iw << inverse(inverse(self.G_iw) + self.Sigma_iw)
        #self.G_iw = self.G0_iw * self.mix + self.G_iw * (1.0 - self.mix)

    def solve_bethe_lattice(self, U):
        self.G_tau = iw_to_tau_dlr(self.G0_iw)
        self.Sigma_tau << (U**2) * self.G_tau * self.G_tau * self.G_tau
        self.Sigma_iw = tau_to_iw_dlr(self.Sigma_tau)

        self.G_iw = inverse(inverse(self.G0_iw) - self.Sigma_iw)
        self.G0_iw << inverse( iOmega_n - t**2 * self.G_iw )
        #self.G_iw = self.G0_iw * self.mix + self.G_iw * (1.0 - self.mix)

    def loop(self, U, bethe_lattice=False):
        for i in range(self.max_loops):
            G_iw_old = self.G_iw.copy()
            if bethe_lattice:
                self.solve_bethe_lattice(U)
            else:
                self.solve(U)
            err = abs(self.G_iw.data - G_iw_old.data).max()
            print("IPT loop %d, err = %.3e" % (i+1, err))
            if err < self.tol:
                break

