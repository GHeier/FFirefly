from firefly.diagram import Diagram
from triqs.gf.meshes import MeshReFreq
from triqs.gf import Gf, inverse
import numpy as np

class IPTSolver_real:
    def __init__(self, H, mu, n_loops=100, mix=0.10, tol=1e-6,
                 w_min=-6.0, w_max=6.0, n_w=1000, eta=0.01):
        """
        Initialize real-axis IPT solver.

        Args:
            H: Hamiltonian callable that takes Sigma and mu
            mu: Chemical potential
            n_loops: Maximum number of SCF iterations
            mix: Mixing parameter for self-consistency (0 < mix <= 1)
            tol: Convergence tolerance
            w_min: Minimum real frequency
            w_max: Maximum real frequency
            n_w: Number of frequency points
            eta: Small imaginary broadening (positive)
        """
        self.H = H
        self.max_loops = n_loops
        self.mix = mix
        self.tol = tol
        self.mu = mu
        self.eta = eta

        # Real frequency mesh with imaginary offset
        # For retarded Green's function: G^R(ω) = G(ω + i*eta)
        refreq_mesh = MeshReFreq(window=(w_min, w_max), n_w=n_w)

        # Initialize Green's functions on real axis
        G_w = Gf(mesh=refreq_mesh, target_shape=[1,1])

        # We don't use Diagram wrapper for real-axis since Diagram expects DLR meshes
        # Work directly with TRIQS Gf objects
        self.G_weiss = G_w.copy()
        self.Sigma_loc = G_w.copy()
        self.Sigma_loc.zero()

        if H is not None:
            # Initialize Weiss field from Hamiltonian
            self.G_weiss << H(Sigma=self.Sigma_loc, mu=self.mu)

        self.G_loc = self.G_weiss.copy()

        # Store mesh points for later use
        self.w_points = np.array([float(w.real) for w in refreq_mesh])

    def get_IPT_Sigma(self, U):
        """
        Calculate IPT self-energy on real axis.

        At T=0, the IPT self-energy simplifies to:
        Σ(ω) = U² G(ω)³

        This is exact for the infinite-U atomic limit and gives reasonable
        results for moderate U in the Bethe lattice.
        """
        G_data = self.G_weiss.data[:, 0, 0]

        # IPT formula: Σ(ω) = U² G(ω)³
        # For real frequencies with small eta, this captures the main physics
        self.Sigma_loc.data[:, 0, 0] = (U**2) * G_data * G_data * G_data

        return self.Sigma_loc

    def set_Weiss(self):
        """Update Weiss field from Dyson equation."""
        self.G_weiss << inverse(inverse(self.G_loc) + self.Sigma_loc)

    def solve(self, U):
        """
        Perform one IPT iteration step.

        1. Calculate self-energy from current Weiss field
        2. Mix with previous self-energy
        3. Update local Green's function via Dyson equation
        4. Update Weiss field for next iteration
        """
        # Get new self-energy
        Sigma_new = self.get_IPT_Sigma(U)

        # Mix self-energies
        self.Sigma_loc << self.mix * Sigma_new + (1.0 - self.mix) * self.Sigma_loc

        # Update local Green's function via Dyson equation
        if self.H is not None:
            self.G_loc << self.H(Sigma=self.Sigma_loc, mu=self.mu)

        # Update Weiss field
        self.set_Weiss()

    def solve_bethe_lattice(self, U):
        """
        Solve IPT for Bethe lattice on real axis.

        For Bethe lattice: G_weiss(ω) = 1/(ω + iη - t² G_loc(ω))
        """
        # Calculate self-energy with mixing for stability
        G_data = self.G_weiss.data[:, 0, 0]
        Sigma_new = (U**2) * G_data * G_data * G_data

        # Mix self-energy
        self.Sigma_loc.data[:, 0, 0] = (
            self.mix * Sigma_new + (1.0 - self.mix) * self.Sigma_loc.data[:, 0, 0]
        )

        # Update local Green's function
        self.G_loc << inverse(inverse(self.G_weiss) - self.Sigma_loc)

        # Update Weiss field for Bethe lattice
        # G_0^{-1}(ω) = ω + iη - t² G(ω)
        t = 1.0
        w_mesh = self.G_weiss.mesh

        # Build frequency mesh with small imaginary part
        for idx, w in enumerate(w_mesh):
            omega = complex(w.real, self.eta)
            G_val = self.G_loc.data[idx, 0, 0]
            # Avoid division by small numbers
            denom = omega - t**2 * G_val
            if abs(denom) < 1e-10:
                denom = 1e-10 * (1.0 + 0.0j)
            self.G_weiss.data[idx, 0, 0] = 1.0 / denom

    def loop(self, U, bethe_lattice=False):
        """
        Run self-consistency loop until convergence.

        Args:
            U: Hubbard interaction strength
            bethe_lattice: If True, use Bethe lattice self-consistency
        """
        for i in range(self.max_loops):
            G_old = self.G_loc.data.copy()

            if bethe_lattice:
                self.solve_bethe_lattice(U)
            else:
                self.solve(U)

            # Check convergence
            err = np.abs(self.G_loc.data - G_old).max()
            print("IPT loop %d, err = %.3e" % (i+1, err))

            if err < self.tol:
                print("Converged after %d iterations" % (i+1))
                break
        else:
            print("Warning: Did not converge after %d iterations" % self.max_loops)

    def get_spectral_function(self):
        """
        Extract spectral function A(ω) = -1/π Im[G(ω + iη)].

        Returns:
            w_points: Array of frequency points
            A_w: Spectral function values
        """
        A_w = -self.G_loc.data[:, 0, 0].imag / np.pi
        return self.w_points, A_w

    def save_results(self, prefix='ipt_real'):
        """
        Save Green's function and self-energy to files.

        Args:
            prefix: Filename prefix
        """
        import firefly as fly

        # Save local Green's function
        G_data = self.G_loc.data[:, 0, 0]
        fly.save_data(f"{prefix}_G_loc.h5", G_data,
                     mesh=None, domain=None, w_points=self.w_points)

        # Save self-energy
        Sigma_data = self.Sigma_loc.data[:, 0, 0]
        fly.save_data(f"{prefix}_Sigma_loc.h5", Sigma_data,
                     mesh=None, domain=None, w_points=self.w_points)

        # Save spectral function
        w, A = self.get_spectral_function()
        fly.save_data(f"{prefix}_spectral.h5", A,
                     mesh=None, domain=None, w_points=w)

        print(f"Results saved with prefix: {prefix}")
