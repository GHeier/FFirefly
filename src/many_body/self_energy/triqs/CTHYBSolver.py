from triqs.gf import inverse, iOmega_n
from triqs.operators import n as n_op
from triqs_cthyb import Solver
from scipy.optimize import brentq
import numpy as np


def effective_mass(Sigma_block):
    """m*/m = Z^-1, estimated from the slope of Im Sigma(iwn) at the two lowest
    positive Matsubara frequencies (same formula used in dmft2.py, for comparability)."""
    data = Sigma_block.data[:, 0, 0]
    n_iw = data.shape[0] // 2
    iw = np.array([complex(w).imag for w in Sigma_block.mesh])
    dImSigma = (data[n_iw + 1].imag - data[n_iw].imag) / (iw[n_iw + 1] - iw[n_iw])
    return 1.0 - dImSigma


class ParamagBlockGf:
    """Thin wrapper around a paramagnetic (up==down) BlockGf on CTHYB's plain
    MeshImFreq. Supports both the dict-style block access (obj['up'], iteration)
    that CTHYBSolver/run.py use internally, AND the .obj_w/.w_points/.copy()/.zero()
    interface that IPTSolver's Diagram objects expose -- so Bubble_DMFTSolver (or
    anything else written against IPTSolver) can treat CTHYBSolver's G_weiss/
    Sigma_imp/G_loc identically to IPTSolver's, without caring which is active."""

    def __init__(self, block_gf):
        self.block_gf = block_gf

    def __getitem__(self, block):
        return self.block_gf[block]

    def __iter__(self):
        return iter(self.block_gf)

    def copy(self):
        return ParamagBlockGf(self.block_gf.copy())

    def zero(self):
        self.block_gf.zero()

    @property
    def obj_w(self):
        # Paramagnetic: up == down, so either block represents the physical object.
        return self.block_gf['up']

    @property
    def w_points(self):
        return np.array([complex(w).imag for w in self.block_gf['up'].mesh], dtype=np.float32)


class CTHYBSolver:
    def __init__(self, beta, H, mu, n_loops=100, mix=0.10, tol=1e-2, w_max=1.2*6, eps=1e-14,
                 n_iw=1025, n_tau=10001,
                 n_cycles=200000, length_cycle=100, n_warmup_cycles=10000, n=None):
        self.beta = beta
        self.H = H
        self.max_loops = n_loops
        self.mix = mix
        self.tol = tol
        self.mu = mu
        self.n = n   # target electron density per spin; None disables mu-search (fixed mu)

        # CTHYB QMC run parameters, passed straight through to S.solve() each iteration
        self.n_cycles = n_cycles
        self.length_cycle = length_cycle
        self.n_warmup_cycles = n_warmup_cycles

        # CTHYB impurity solver: paramagnetic single-orbital Hubbard model (spin up/down blocks).
        # Everything below lives directly on the solver's own equally-spaced MeshImFreq --
        # no DLR mesh, no conversion back and forth.
        self.S = Solver(beta=beta, gf_struct=[('up', 1), ('down', 1)], n_iw=n_iw, n_tau=n_tau)

        self.G_weiss = ParamagBlockGf(self.S.G0_iw.copy())
        for block, g in self.G_weiss:
            g << inverse(iOmega_n + mu)
        self.G_weiss_old = self.G_weiss.copy()

        self.Sigma_imp = self.G_weiss.copy()
        self.Sigma_imp.zero()
        self.G_loc = self.G_weiss.copy()

        if H is not None:
            if self.n is not None:
                # Find the correct mu for the target density FIRST (at Sigma=0), same
                # ordering as FLEXSolver.find_mu_for_density used in __init__.
                self.mu = self.find_mu_for_density(self.n)
                print(f"Initial mu set to {self.mu:.4f} for n = {self.n:.6f}")
            mu_matrix = self.mu * np.eye(1)   # mu must be a matrix for matrix-valued Gfs
            for block, g in self.G_loc:
                g << H(Sigma=self.Sigma_imp[block], mu=mu_matrix)
            self.G_weiss = self.G_loc.copy()

    def calc_electron_density(self, mu):
        """Electron density per spin (paramagnetic) at the given mu, for the CURRENT
        (fixed) Sigma_imp -- mirrors FLEXSolver.calc_electron_density, but reuses the
        already-solved Sigma_imp instead of needing a fresh QMC solve per trial mu."""
        mu_matrix = mu * np.eye(1)
        G_trial = self.G_loc['up'].copy()
        G_trial << self.H(Sigma=self.Sigma_imp['up'], mu=mu_matrix)
        return G_trial.density()[0, 0].real

    def find_mu_for_density(self, n_target):
        """Root-find mu (via bisection) so that calc_electron_density(mu) == n_target,
        holding Sigma_imp fixed. Same bracketing strategy as FLEXSolver.find_mu_for_density."""
        n = self.calc_electron_density(self.mu)
        if abs(n - n_target) < 1e-4:
            return self.mu

        def density_error(mu):
            return self.calc_electron_density(mu) - n_target

        bound = 10.0
        if abs(n - n_target) < 1e-2:
            bound = 2.0
        mu_min = self.mu - bound
        mu_max = self.mu + bound

        try:
            mu = brentq(density_error, mu_min, mu_max, xtol=1e-4)
            self.mu = mu
            return mu
        except ValueError:
            print(f"ERROR: Could not find mu for n = {n_target}")
            print(f"Check that the target density is achievable in the range [{mu_min:.4f}, {mu_max:.4f}]")
            raise

    def shift_mu_to_target_density(self, n_target):
        mu = self.find_mu_for_density(n_target)
        mu_matrix = mu * np.eye(1)
        for block, g in self.G_loc:
            g << self.H(Sigma=self.Sigma_imp[block], mu=mu_matrix)

    def _solve_cthyb(self, U):
        """Run CTHYB on whatever is currently loaded into self.S.G0_iw and return the
        paramagnetically symmetrized Sigma_iw, wrapped for Diagram-style access."""
        self.S.solve(
            h_int=U * n_op('up', 0) * n_op('down', 0),
            n_cycles=self.n_cycles,
            length_cycle=self.length_cycle,
            n_warmup_cycles=self.n_warmup_cycles,
        )

        Sigma = self.S.Sigma_iw.copy()
        Sigma_sym = 0.5 * (Sigma['up'] + Sigma['down'])
        for block, sg in Sigma:
            sg << Sigma_sym
        return ParamagBlockGf(Sigma)

    def get_CTHYB_Sigma(self, U):
        """Solve the impurity problem with CTHYB using self.G_weiss, directly on
        CTHYB's own MeshImFreq -- same as dmft2.py."""
        self.S.G0_iw << self.G_weiss.block_gf
        return self._solve_cthyb(U)

    def get_IPT_Sigma(self, U, G_weiss):
        """Same name/signature as IPTSolver.get_IPT_Sigma, so Bubble_DMFTSolver can
        call either solver identically. Unlike IPT's cheap analytic G_weiss^3
        formula, this runs a full CTHYB solve using the GIVEN G_weiss (not
        self.G_weiss, and without touching self.Sigma_imp) -- much more expensive,
        but functionally equivalent from the caller's point of view."""
        g_weiss_gf = G_weiss.obj_w
        for block in self.S.G0_iw.indices:
            self.S.G0_iw[block] << g_weiss_gf
        return self._solve_cthyb(U)

    def set_G_loc(self, G_loc):
        """Same purpose as IPTSolver.set_G_loc: overwrite G_loc from an externally
        computed (paramagnetic) source, e.g. a combined Bubble+DMFT lattice update.
        Broadcasts into both spin blocks since CTHYB tracks up/down explicitly."""
        g = G_loc.obj_w
        for block, gb in self.G_loc:
            gb << g

    def set_Weiss(self):
        self.G_weiss_old = self.G_weiss.copy()
        G_weiss_new = self.G_weiss.copy()
        for block, g in G_weiss_new:
            g << inverse(inverse(self.G_loc[block]) + self.Sigma_imp[block])
        for block, g in self.G_weiss:
            g << self.mix * G_weiss_new[block] + (1.0 - self.mix) * self.G_weiss_old[block]

    def solve(self, U):
        self.Sigma_imp = self.get_CTHYB_Sigma(U)

        # Dyson (re-tune mu to hit the target density if one was given, else fixed mu)
        if self.n is not None:
            self.shift_mu_to_target_density(self.n)
        else:
            mu_matrix = self.mu * np.eye(1)
            for block, g in self.G_loc:
                g << self.H(Sigma=self.Sigma_imp[block], mu=mu_matrix)
        self.set_Weiss()

    def solve_bethe_lattice(self, U):
        self.Sigma_imp = self.get_CTHYB_Sigma(U)

        t = 1
        for block, g in self.G_loc:
            g << inverse(inverse(self.G_weiss[block]) - self.Sigma_imp[block])
        for block, g in self.G_weiss:
            g << inverse(iOmega_n - t**2 * self.G_loc[block])

    def loop(self, U, bethe_lattice=False):
        old_m = {'up': 0.0, 'down': 0.0}
        break_cond = False
        for i in range(self.max_loops):
            G_iw_prev = {block: g.data.copy() for block, g in self.G_weiss}

            if bethe_lattice:
                self.solve_bethe_lattice(U)
            else:
                self.solve(U)

            err = max(abs(g.data - G_iw_prev[block]).max() for block, g in self.G_weiss)
            print("CTHYB loop %d, err = %.3e" % (i+1, err))

            for block, _ in self.Sigma_imp:
                m = effective_mass(self.Sigma_imp[block])
                print(f"    m*/m [{block}] = {m:.4f}")
                m_err = abs(old_m[block] - m)
                print(f"    m*/m err = {m_err:.4f}")
                old_m[block] = m
                if m_err < 1e-4:
                    break_cond = True

            if self.n is None:
                # mu is fixed (no mu_from_n target) -- read off the density this mu
                # actually gives at the current Sigma_imp, purely for visibility.
                n_actual = self.calc_electron_density(self.mu)
                print(f"    n [mu={self.mu:.4f}] = {n_actual:.6f}")
            else:
                # mu_from_n target given -- read off the mu that was found for it.
                print(f"    mu [n={self.n:.6f}] = {self.mu:.4f}")

            if err < self.tol or break_cond:
                break

    def save_results(self, prefix, renorm=None):
        """Save G_loc/Sigma_imp/spectral function directly via fly.save_data, same
        pattern as IPTSolver_real.save_results -- Diagram expects a DLR mesh, which
        CTHYB's plain MeshImFreq isn't, so we bypass Diagram here too."""
        import firefly as fly
        from triqs.gf import Gf
        from triqs.gf.meshes import MeshReFreq

        G_data = self.G_loc['up'].data[:, 0, 0]
        fly.save_data(f"{prefix}_G_iw.h5", G_data, mesh=None, domain=None,
                       w_points=self.G_loc.w_points)

        Sigma_data = self.Sigma_imp['up'].data[:, 0, 0]
        fly.save_data(f"{prefix}_self_energy.h5", Sigma_data, mesh=None, domain=None,
                       w_points=self.Sigma_imp.w_points)
        fly.save_data(f"{prefix}_sigma_iw.h5", Sigma_data, mesh=None, domain=None,
                       w_points=self.Sigma_imp.w_points)

        # Spectral function via Pade continuation directly from the native MeshImFreq
        # G_loc -- no DLR round-trip needed since it's already on a standard imfreq mesh.
        w_min, w_max, n_w = -5.0, 5.0, 500
        Gw = Gf(mesh=MeshReFreq(window=(w_min, w_max), n_w=n_w), target_shape=[1, 1])
        Gw.set_from_pade(self.G_loc['up'])
        A_w = -Gw.data[:, 0, 0].imag / np.pi
        w_points_real = np.linspace(w_min, w_max, n_w)
        fly.save_data(f"{prefix}_A_w.h5", A_w, mesh=None, domain=None, w_points=w_points_real)

        # Same format as save_DMFT's IPT renormalization output: a bare scalar, no mesh.
        if renorm is not None:
            fly.save_data(f"{prefix}_renormalization.h5", np.array([renorm]),
                           mesh=None, domain=None, w_points=None)

        print(f"CTHYB results saved with prefix: {prefix}")


# ---------------------------------------------------------------------------
# Old DLR-mesh implementation (superseded -- CTHYB works natively on a plain
# MeshImFreq, so the DLR <-> MeshImFreq round-trip below is unnecessary).
# Kept for reference.
# ---------------------------------------------------------------------------
#
# from firefly.diagram import Diagram, dot_tr, dot_t
# from triqs.gf.meshes import MeshDLRImFreq, MeshDLRImTime
# from triqs.gf import Gf, make_gf_dlr, make_gf_imfreq, make_gf_dlr_imfreq, fit_gf_dlr
#
# class CTHYBSolver:
#     def __init__(self, beta, H, mu, n_loops=100, mix=0.10, tol=1e-6, w_max=1.2*6, eps=1e-14,
#                  n_iw=1025, n_tau=10001, n_l=80,
#                  n_cycles=200000, length_cycle=100, n_warmup_cycles=10000, measure_G_l=True):
#         self.beta = beta
#         self.H = H
#         self.max_loops = n_loops
#         self.mix = mix
#         self.tol = tol
#         self.mu = mu
#         self.w_max = w_max
#         self.eps = eps
#
#         self.n_cycles = n_cycles
#         self.length_cycle = length_cycle
#         self.n_warmup_cycles = n_warmup_cycles
#         self.measure_G_l = measure_G_l
#
#         dlr_iw_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
#         ws = np.array([float(iw.imag) for iw in dlr_iw_mesh], dtype=np.float32)
#         print("Max w: ", max(ws))
#         print("Min w: ", min(ws))
#
#         G_iw = Gf(mesh=dlr_iw_mesh, target_shape=[1, 1])
#         self.G_weiss = Diagram(G_iw, 'Fermion')
#         self.G_weiss_old = self.G_weiss.copy()
#         self.Sigma_imp = self.G_weiss.copy()
#         self.Sigma_imp.zero()
#
#         if H is not None:
#             mu_matrix = self.mu * np.eye(1)
#             self.G_weiss.obj_w << H(Sigma=self.Sigma_imp.obj_w, mu=mu_matrix)
#         self.G_loc = self.G_weiss.copy()
#
#         self.S = Solver(beta=beta, gf_struct=[('up', 1), ('down', 1)], n_iw=n_iw, n_tau=n_tau, n_l=n_l)
#
#     def _push_Weiss_to_cthyb(self):
#         n_iw_cthyb = self.S.G0_iw['up'].data.shape[0] // 2
#         G0_dlr_coeffs = make_gf_dlr(self.G_weiss.obj_w)
#         G0_imfreq = make_gf_imfreq(G0_dlr_coeffs, n_iw_cthyb)
#         for block, g0 in self.S.G0_iw:
#             g0.data[:] = G0_imfreq.data[:]
#
#     def get_CTHYB_Sigma(self, U):
#         self._push_Weiss_to_cthyb()
#
#         self.S.solve(
#             h_int=U * n_op('up', 0) * n_op('down', 0),
#             n_cycles=self.n_cycles,
#             length_cycle=self.length_cycle,
#             n_warmup_cycles=self.n_warmup_cycles,
#             measure_G_l=self.measure_G_l,
#         )
#
#         G_tau_avg = 0.5 * (self.S.G_tau['up'] + self.S.G_tau['down'])
#         G_dlr_coeffs = fit_gf_dlr(G_tau_avg, self.w_max, self.eps)
#         G_imp = Diagram(make_gf_dlr_imfreq(G_dlr_coeffs), 'Fermion')
#
#         Sigma = self.G_weiss.copy()
#         Sigma.obj_w << inverse(self.G_weiss.obj_w) - inverse(G_imp.obj_w)
#         return Sigma
#
#     def set_Weiss(self):
#         self.G_weiss_old = self.G_weiss.copy()
#         self.G_weiss.obj_w << inverse(inverse(self.G_loc.obj_w) + self.Sigma_imp.obj_w)
#         self.G_weiss.obj_w << self.mix * self.G_weiss_old.obj_w + (1.0 - self.mix) * self.G_weiss.obj_w
#
#     def solve(self, U):
#         self.Sigma_imp = self.get_CTHYB_Sigma(U)
#         mu_matrix = self.mu * np.eye(1)
#         self.G_loc.obj_w << self.H(Sigma=self.Sigma_imp.obj_w, mu=mu_matrix)
#         self.set_Weiss()
#
#     def solve_bethe_lattice(self, U):
#         self.Sigma_imp = self.get_CTHYB_Sigma(U)
#         self.G_loc.obj_w = inverse(inverse(self.G_weiss.obj_w) - self.Sigma_imp.obj_w)
#         t = 1
#         self.G_weiss.obj_w << inverse( iOmega_n - t**2 * self.G_loc.obj_w )
#
#     def loop(self, U, bethe_lattice=False):
#         for i in range(self.max_loops):
#             G_iw_prev = self.G_weiss.obj_w.data.copy()
#             if bethe_lattice:
#                 self.solve_bethe_lattice(U)
#             else:
#                 self.solve(U)
#             err = abs(self.G_weiss.obj_w.data - G_iw_prev).max()
#             print("CTHYB loop %d, err = %.3e" % (i+1, err))
#             if err < self.tol:
#                 break
