from triqs.gf import inverse, iOmega_n
from triqs.operators import n as n_op
from triqs_cthyb import Solver
from triqs.utility import mpi
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
    that CTHYBSolver uses internally, AND the .obj_w/.w_points/.copy()/.zero()
    interface that Diagram exposes -- so the internal native-mesh bookkeeping can
    reuse the same access patterns as the DLR-mesh public attributes."""

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

    def save(self, filename, force_real=False):
        """Same format as Diagram.save's 'w'-space branch, so callers written
        against IPTSolver's Diagram-wrapped G_loc/Sigma_imp/G_weiss can call
        .save() on CTHYBSolver's equivalents without caring which is active."""
        # Every rank reaches this call (CTHYB's own solve() already collectively
        # reduces the QMC result onto all ranks) but only rank 0 should write --
        # concurrent ranks writing the same HDF5 file race and throw H5::FileIException.
        if not mpi.is_master_node():
            return
        import firefly as fly
        data = self.obj_w.data.real if force_real else self.obj_w.data
        obj = np.reshape(data, (data.shape[0],))
        fly.save_data(filename, obj, mesh=None, domain=None, w_points=self.w_points,
                       is_complex=not force_real)
        print(f"Diagram saved to {filename}")


class CTHYBSolver:
    def __init__(self, beta, H, mu, n_loops=100, mix=0.10, tol=1e-2, w_max=1.2*6, eps=1e-14,
                 n_iw=1025, n_tau=10001,
                 n_cycles=200000, length_cycle=100, n_warmup_cycles=10000, n=None, damp=True):
        self.beta = beta
        self.H = H
        self.max_loops = n_loops
        self.mix = mix
        self.tol = tol
        self.mu = mu
        self.n = n   # target electron density per spin; None disables mu-search (fixed mu)
        self.w_max = w_max
        self.eps = eps
        self.damp = damp

        # CTHYB QMC run parameters, passed straight through to S.solve() each iteration
        self.n_cycles = n_cycles
        self.length_cycle = length_cycle
        self.n_warmup_cycles = n_warmup_cycles

        # S.auto_corr_time (from S.solve()) collected here, one entry per CTHYB solve
        # across the whole loop() run, for the end-of-run mean/std report.
        self.auto_corr_times = []

        # CTHYB impurity solver: paramagnetic single-orbital Hubbard model (spin up/down blocks).
        # The inner G_weiss <-> Sigma_imp <-> G_loc self-consistency (solve()/set_Weiss(),
        # driven by loop()) lives ENTIRELY on this native mesh, in the private _G_weiss/
        # _Sigma_imp/_G_loc attributes below -- that's a hard requirement of the QMC backend,
        # and mixing/damping on this mesh is what keeps the loop numerically stable. Round-
        # tripping this inner loop through a DLR mesh every iteration was tried and made the
        # hybridization function drift non-causal over repeated solves until CTHYB's
        # determinant went singular -- see git history. So the native mesh here is never
        # touched by a DLR conversion.
        self.S = Solver(beta=beta, gf_struct=[('up', 1), ('down', 1)], n_iw=n_iw, n_tau=n_tau)

        # If a converged Weiss field from a previous run is sitting on disk, use it
        # (interpolated onto CTHYB's native Matsubara mesh) directly as the initial
        # self._G_weiss instead of starting cold from the Sigma=0 bath -- lets a
        # restarted run skip re-thermalizing CTHYB from a non-interacting start.
        #
        # This used to warm-start from a saved local SELF-ENERGY instead (see git
        # history). Sigma(iwn) = G0^-1(iwn) - G^-1(iwn) is a difference of two
        # near-equal, near-iwn quantities at large iwn, so its absolute QMC noise --
        # small near w=0 -- grows into outright causality violations far out in the
        # tail (confirmed directly: 62 violations of Im[G_weiss(iwn>0)]<=0 from one
        # such file, which is what was making CTHYB's determinant go singular almost
        # immediately on a from-file restart, even after capping the trusted window).
        # G_weiss itself needs none of that: it's a bounded spectral function
        # integrated against 1/(iwn-eps), so it inherits a gentle, universal ~1/iwn
        # tail directly with no subtraction-driven noise blowup, and it's exactly
        # what CTHYB's own G0_iw gets set from every iteration anyway (get_CTHYB_Sigma)
        # -- so loading it needs no reconstruction through Sigma or the lattice H at
        # all, just a direct copy onto this run's native mesh.
        import os
        import firefly as fly
        import firefly.config as cfg
        g_weiss_guess = None
        g_weiss_field = None
        # outdir+prefix+"_self_energy.h5" is the exact filename
        # Bubble_DMFTSolver.save_Bubble_DMFT() writes for the full k-dependent
        # lattice self-energy -- a completely different, incompatible field. Use our
        # own dedicated filename instead (see save_full_g_weiss below).
        g_weiss_file = cfg.outdir + cfg.prefix + "_G_weiss_full.h5"
        if os.path.isfile(g_weiss_file):
            g_weiss_field = fly.Field_C(g_weiss_file)
            if len(g_weiss_field.mesh) != 0:
                # Defensive: nothing currently writes a k-dependent field under this
                # name, but silently feeding a k-point array into a field that expects
                # a bare frequency there returns one bogus (mesh-mismatched) constant,
                # stamped identically onto every in-range Matsubara point regardless of
                # sign -- since a real causal G(-iwn) must be G(iwn)'s complex
                # conjugate, not an identical value, that alone breaks G0_iw's
                # G(iw)=G(-iw)* symmetry outright (the "WARNING: S.G0_iw violates
                # fundamental Green Function property" symmetrization message
                # downstream is what catches it).
                print(f"Found {g_weiss_file}, but it has k-dependence (not a local "
                      f"impurity Green's function) -- skipping as initial guess.")
                g_weiss_field = None

        if g_weiss_field is not None:
            print(f"Found existing G_weiss at {g_weiss_file}, using as initial guess")
            # iwn can reach further on this run's native mesh than on the one that
            # saved the file (different beta/n_iw), so evaluating outside the file's
            # own range still needs guarding against "w_val out of bounds". Unlike
            # the old Sigma-based guess, zero is the WRONG fallback here -- G_weiss's
            # own tail decays like 1/iwn, not to zero -- so fall back to the analytic
            # Sigma=0 bath G0(iwn)=1/(iwn+mu) instead, matching the known-good
            # cold-start behavior out there. A plain MeshImFreq is exactly
            # mirror-symmetric by construction (unlike a DLR mesh), so a plain
            # symmetric min/max bound is exact here, not a hack.
            w_points = g_weiss_field.w_points
            w_bound = min(abs(w_points.min()), abs(w_points.max())) if len(w_points) else 0.0
            w_lo, w_hi = -w_bound, w_bound
            g_weiss_guess = self.S.G0_iw.copy()
            for block, g in g_weiss_guess:
                iw_native = np.array([complex(iw).imag for iw in g.mesh], dtype=np.float32)
                in_range = (iw_native >= w_lo) & (iw_native <= w_hi)
                g.data[~in_range, 0, 0] = 1.0 / (1j * iw_native[~in_range] + mu)
                g.data[in_range, 0, 0] = g_weiss_field(iw_native[in_range])
            g_weiss_guess = ParamagBlockGf(g_weiss_guess)

        self._G_weiss = ParamagBlockGf(self.S.G0_iw.copy())
        if g_weiss_guess is not None:
            self._G_weiss = g_weiss_guess
        else:
            for block, g in self._G_weiss:
                g << inverse(iOmega_n + mu)
        self._G_weiss_old = self._G_weiss.copy()

        # No warm start for Sigma_imp -- the loaded G_weiss above is the only thing
        # bootstrapping this run; the first real CTHYB solve produces a fresh Sigma.
        self._Sigma_imp = self._G_weiss.copy()
        self._Sigma_imp.zero()
        self._G_loc = self._G_weiss.copy()

        if H is not None:
            if self.n is not None:
                # Find the correct mu for the target density FIRST (at Sigma=0), same
                # ordering as FLEXSolver.find_mu_for_density used in __init__.
                self.mu = self.find_mu_for_density(self.n)
                print(f"Initial mu set to {self.mu:.4f} for n = {self.n:.6f}")
            mu_matrix = self.mu * np.eye(1)   # mu must be a matrix for matrix-valued Gfs
            for block, g in self._G_loc:
                g << H(Sigma=self._Sigma_imp[block], mu=mu_matrix)
            # Only replace G_weiss with the lattice-derived G_loc when we didn't just
            # warm-start it from file -- otherwise this would throw the warm start away.
            if g_weiss_guess is None:
                self._G_weiss = self._G_loc.copy()

        # Public exposure, on the shared DLR mesh (beta/w_max/eps) the rest of the
        # Bubble+DMFT pipeline uses -- same convention as IPTSolver.Sigma_imp/G_loc, so
        # Bubble_DMFTSolver/run.py can treat either solver identically. Zero-initialized
        # here; filled in for real by _project_to_dlr() once loop() converges.
        from firefly.diagram import Diagram
        from triqs.gf import Gf
        from triqs.gf.meshes import MeshDLRImFreq
        dlr_mesh = MeshDLRImFreq(beta=beta, statistic='Fermion', w_max=w_max, eps=eps)
        self.Sigma_imp = Diagram(Gf(mesh=dlr_mesh, target_shape=[1, 1]), 'Fermion')
        self.Sigma_imp.zero()
        self.G_loc = self.Sigma_imp.copy()
        if H is not None:
            mu_matrix = self.mu * np.eye(1)
            self.G_loc.obj_w << H(Sigma=self.Sigma_imp.obj_w, mu=mu_matrix)

    def calc_electron_density(self, mu):
        """Electron density per spin (paramagnetic) at the given mu, for the CURRENT
        (fixed) Sigma_imp -- mirrors FLEXSolver.calc_electron_density, but reuses the
        already-solved Sigma_imp instead of needing a fresh QMC solve per trial mu."""
        mu_matrix = mu * np.eye(1)
        G_trial = self._G_loc['up'].copy()
        G_trial << self.H(Sigma=self._Sigma_imp['up'], mu=mu_matrix)
        return G_trial.density()[0, 0].real

    def find_mu_for_density(self, n_target, damp=False):
        """Root-find mu (via bisection) so that calc_electron_density(mu) == n_target,
        holding Sigma_imp fixed. Same bracketing strategy as FLEXSolver.find_mu_for_density.

        damp=True blends the bisection result with the previous self.mu using the same
        self.mix as set_Weiss() uses for G_weiss. Sigma_imp is a fresh, noisy QMC estimate
        each outer iteration, and calc_electron_density(mu) is evaluated directly against
        that noisy Sigma_imp (not a damped/averaged one) -- so an undamped bisection result
        fully inherits that iteration's QMC shot noise, which shows up as mu/m*/m still
        oscillating even once err (the damped G_weiss's own convergence metric) is well
        under tol. damp=False (the default) preserves the original undamped behavior, used
        for the one-off Sigma=0 search in __init__ where there's no previous mu to blend with."""
        n = self.calc_electron_density(self.mu)
        if abs(n - n_target) < 1e-4:
            return self.mu

        def density_error(mu):
            return self.calc_electron_density(mu) - n_target

        bound = 2.0 if abs(n - n_target) < 1e-2 else 10.0
        # A nonzero Sigma_imp (from CTHYB, or from an initial-guess file loaded in
        # __init__) can flatten n(mu) into a Mott-like plateau near n_target, so the
        # density already being close at self.mu (the bound=2.0 case above) doesn't
        # guarantee a sign change within +-2 of it -- widen the bracket geometrically
        # until density_error actually changes sign, instead of handing brentq a
        # window that can't bracket a root.
        max_bound = 40.0
        mu_min, mu_max = self.mu - bound, self.mu + bound
        e_min, e_max = density_error(mu_min), density_error(mu_max)
        while e_min * e_max >= 0 and bound < max_bound:
            bound *= 2
            mu_min, mu_max = self.mu - bound, self.mu + bound
            e_min, e_max = density_error(mu_min), density_error(mu_max)
            print(f"    mu bracket widened to +-{bound:.1f} (n({mu_min:.3f})={e_min + n_target:.4f}, "
                  f"n({mu_max:.3f})={e_max + n_target:.4f})")

        try:
            mu = brentq(density_error, mu_min, mu_max, xtol=1e-4)
            self.mu = self.mix * mu + (1.0 - self.mix) * self.mu if damp else mu
            return mu
        except ValueError:
            print(f"ERROR: Could not find mu for n = {n_target}")
            print(f"Check that the target density is achievable in the range [{mu_min:.4f}, {mu_max:.4f}]")
            raise

    def shift_mu_to_target_density(self, n_target):
        mu = self.find_mu_for_density(n_target, damp=self.damp)
        mu_matrix = mu * np.eye(1)
        for block, g in self._G_loc:
            g << self.H(Sigma=self._Sigma_imp[block], mu=mu_matrix)

    def _solve_cthyb(self, U):
        """Run CTHYB on whatever is currently loaded into self.S.G0_iw and return the
        paramagnetically symmetrized Sigma_iw, wrapped for Diagram-style access."""
        self.S.solve(
            h_int=U * n_op('up', 0) * n_op('down', 0),
            n_cycles=self.n_cycles,
            length_cycle=self.length_cycle,
            n_warmup_cycles=self.n_warmup_cycles,
        )
        self.auto_corr_times.append(self.S.auto_corr_time)

        Sigma = self.S.Sigma_iw.copy()
        Sigma_sym = 0.5 * (Sigma['up'] + Sigma['down'])
        for block, sg in Sigma:
            sg << Sigma_sym
        return ParamagBlockGf(Sigma)

    def get_CTHYB_Sigma(self, U):
        """Solve the impurity problem with CTHYB using self._G_weiss, directly on
        CTHYB's own MeshImFreq -- same as dmft2.py."""
        self.S.G0_iw << self._G_weiss.block_gf
        return self._solve_cthyb(U)

    def get_IPT_Sigma(self, U, G_weiss):
        """Same name/signature as IPTSolver.get_IPT_Sigma, so Bubble_DMFTSolver can
        call either solver identically for a one-off external solve. NOTE: G_weiss
        here must already be on CTHYB's native MeshImFreq (e.g. from another
        CTHYBSolver's _G_weiss) -- unlike self.Sigma_imp/self.G_loc, this does not
        accept a DLR-mesh Diagram (no DLR<->native conversion is done here; see
        class docstring for why that conversion is kept out of any per-solve path)."""
        g_weiss_gf = G_weiss.obj_w
        for block in self.S.G0_iw.indices:
            self.S.G0_iw[block] << g_weiss_gf
        return self._solve_cthyb(U)

    def set_G_loc(self, G_loc):
        """Same purpose as IPTSolver.set_G_loc: overwrite G_loc from an externally
        computed (paramagnetic) source, e.g. a combined Bubble+DMFT lattice update.
        G_loc must be on CTHYB's native MeshImFreq (see get_IPT_Sigma note above);
        broadcasts into both spin blocks since CTHYB tracks up/down explicitly."""
        g = G_loc.obj_w
        for block, gb in self._G_loc:
            gb << g

    def set_Weiss(self):
        self._G_weiss_old = self._G_weiss.copy()
        G_weiss_new = self._G_weiss.copy()
        for block, g in G_weiss_new:
            g << inverse(inverse(self._G_loc[block]) + self._Sigma_imp[block])
        for block, g in self._G_weiss:
            g << self.mix * G_weiss_new[block] + (1.0 - self.mix) * self._G_weiss_old[block]

    def solve(self, U):
        self._Sigma_imp = self.get_CTHYB_Sigma(U)

        # Dyson (re-tune mu to hit the target density if one was given, else fixed mu)
        if self.n is not None:
            self.shift_mu_to_target_density(self.n)
        else:
            mu_matrix = self.mu * np.eye(1)
            for block, g in self._G_loc:
                g << self.H(Sigma=self._Sigma_imp[block], mu=mu_matrix)
        self.set_Weiss()

    def solve_bethe_lattice(self, U):
        self._Sigma_imp = self.get_CTHYB_Sigma(U)

        t = 1
        for block, g in self._G_loc:
            g << inverse(inverse(self._G_weiss[block]) - self._Sigma_imp[block])
        for block, g in self._G_weiss:
            g << inverse(iOmega_n - t**2 * self._G_loc[block])

    def loop(self, U, bethe_lattice=False):
        m_vals, mu_vals, errs = [], [], []
        old_m = {'up': 0.0, 'down': 0.0}
        stable_iters = 0
        for i in range(self.max_loops):
            print("Iteration ", i)
            G_iw_prev = {block: g.data.copy() for block, g in self._G_weiss}

            if bethe_lattice:
                self.solve_bethe_lattice(U)
            else:
                self.solve(U)

            err = max(abs(g.data - G_iw_prev[block]).max() for block, g in self._G_weiss)
            print("CTHYB loop %d, err = %.3e" % (i+1, err))
            errs.append(err)

            all_blocks_stable = True
            for block, _ in self._Sigma_imp:
                m = effective_mass(self._Sigma_imp[block])
                #print(f"    m*/m [{block}] = {m:.4f}")
                m_err = abs(old_m[block] - m)
                print(f"    m*/m err = {m_err:.4f}")
                old_m[block] = m
                if m_err >= 1e-4:
                    all_blocks_stable = False
            # Require sustained agreement across blocks, not a single lucky
            # pair of noisy CTHYB estimates landing close by chance.
            stable_iters = stable_iters + 1 if all_blocks_stable else 0
            m_vals.append(old_m['up'])

            if self.n is None:
                # mu is fixed (no mu_from_n target) -- read off the density this mu
                # actually gives at the current Sigma_imp, purely for visibility.
                n_actual = self.calc_electron_density(self.mu)
                print(f"    n [mu={self.mu:.4f}] = {n_actual:.6f}")
            else:
                # mu_from_n target given -- read off the mu that was found for it.
                print(f"    mu [n={self.n:.6f}] = {self.mu:.4f}")
            mu_vals.append(self.mu)

            print("Last 10 Snapshot")
            print(" m*     |   mu")
            l = len(m_vals)
            for j in range(max(0, l - 11), l):
                print(f"{m_vals[j]:0.3f} | {mu_vals[j]:0.3f}")

            window_stable = (err < self.tol and len(m_vals) >= 2
                              and abs(np.max(m_vals[-10:]) - np.min(m_vals[-10:])) < self.tol
                              and abs(mu_vals[-1] - mu_vals[-2]) < self.tol)
            if window_stable or stable_iters >= 10:
                break
        print("Iterations: ", i)
        self.print_cthyb_stats(m_vals, mu_vals, errs)
        print(f"Auto-correlation time: mean = {np.mean(self.auto_corr_times):.4f}, "
              f"std = {np.std(self.auto_corr_times):.4f}")
        self._project_to_dlr()

    def _project_to_dlr(self):
        """One-time projection of the converged (native-mesh) impurity solution onto
        the shared DLR mesh (beta/w_max/eps), refreshing the public self.Sigma_imp/
        self.G_loc. Called once, after loop()'s inner self-consistency has already
        converged -- NOT part of that inner loop, so there's no per-iteration DLR
        round-trip noise to compound (that's what caused the singular-determinant
        instability previously; see class docstring).

        Fits G_imp from G_tau (dense, QMC-averaged over 20000 measures -- well suited
        to a tau-domain least-squares DLR fit) rather than fitting Sigma_iw directly:
        Sigma has sharper structure than G and fitting it directly would amplify QMC
        noise. G_weiss is Fourier-transformed to tau (it's the analytic mixing result
        of the loop above, not a QMC measurement, so there's no noise to fit around)
        purely so Sigma_imp can be extracted via the impurity Dyson equation
        G_weiss^-1 - G_imp^-1, same relation set_Weiss() uses on the native mesh.
        """
        from firefly.diagram import Diagram
        from triqs.gf import fit_gf_dlr, make_gf_dlr_imfreq, make_gf_from_fourier

        G_tau_avg = 0.5 * (self.S.G_tau['up'] + self.S.G_tau['down'])
        G_imp_coeffs = fit_gf_dlr(G_tau_avg, self.w_max, self.eps)
        G_imp_dlr = Diagram(make_gf_dlr_imfreq(G_imp_coeffs), 'Fermion')

        G_weiss_tau = make_gf_from_fourier(self._G_weiss['up'], len(G_tau_avg.mesh))
        G_weiss_coeffs = fit_gf_dlr(G_weiss_tau, self.w_max, self.eps)
        G_weiss_dlr = Diagram(make_gf_dlr_imfreq(G_weiss_coeffs), 'Fermion')

        self.Sigma_imp.obj_w << inverse(G_weiss_dlr.obj_w) - inverse(G_imp_dlr.obj_w)

        mu_matrix = self.mu * np.eye(1)
        self.G_loc.obj_w << self.H(Sigma=self.Sigma_imp.obj_w, mu=mu_matrix)

    def print_cthyb_stats(self, m_vals, mu_vals, errs):
        print("\n-------- CTHYB STATS --------")
        print(f" m*/m's  |  mu's  |  errors")
        for i in range(len(m_vals)):
            print(f"{i} | {m_vals[i]:^9.4f} | {mu_vals[i]:^8.4f} | {errs[i]:^8.4e}")
        print()

    def save_full_self_energy(self, filename):
        """Save the local impurity self-energy at full native-mesh resolution
        (thousands of raw Matsubara points straight from the QMC solve, on CTHYB's
        own MeshImFreq -- see self._Sigma_imp). This is the file __init__ looks for
        as an initial-guess warm start on a later run, so anything that runs a CTHYB
        solve without going through save_results (e.g. Bubble_DMFTSolver, which
        writes its own combined k-dependent lattice self-energy to
        outdir+prefix+"_self_energy.h5" -- a different, incompatible file this
        solver's own init guard now refuses to load) should call this too."""
        self._Sigma_imp.save(filename)

    def save_full_g_weiss(self, filename):
        """Save the Weiss/bath field at full native-mesh resolution (self._G_weiss --
        exactly what get_CTHYB_Sigma() feeds into CTHYB's own G0_iw every iteration).
        This is the file __init__ looks for as an initial-guess warm start on a later
        run (see the note there for why G_weiss, not Sigma_imp, is what's reused).
        Anything that runs a CTHYB solve without going through save_results (e.g.
        Bubble_DMFTSolver) should call this too."""
        self._G_weiss.save(filename)

    def save_results(self, prefix, renorm=None):
        """Save G_loc/Sigma_imp/spectral function via Diagram's own save/save_spectral,
        now that both are DLR-mesh Diagram objects (see _project_to_dlr) -- same pattern
        as IPTSolver's DMFT save path in run.py's save_DMFT."""
        self.G_loc.save(f"{prefix}_G_iw.h5")
        self.Sigma_imp.save(f"{prefix}_self_energy.h5")
        self.Sigma_imp.save(f"{prefix}_sigma_iw.h5")
        self.G_loc.save_spectral(f"{prefix}_A_w.h5")

        # Also save the full-resolution self-energy (see save_full_self_energy above)
        # alongside the coarse DLR fit -- the DLR fit's outermost nodes are the
        # least-constrained coefficients of that fit and can come out as outliers,
        # so anything reusing this self-energy at full resolution should prefer this.
        self.save_full_self_energy(f"{prefix}_self_energy_full.h5")

        # The actual initial-guess warm-start file __init__ looks for (see there).
        self.save_full_g_weiss(f"{prefix}_G_weiss_full.h5")

        if renorm is not None and mpi.is_master_node():
            import firefly as fly
            fly.save_data(f"{prefix}_renormalization.h5", np.array([renorm]),
                           mesh=None, domain=None, w_points=None)

        if mpi.is_master_node():
            print(f"CTHYB results saved with prefix: {prefix}")
