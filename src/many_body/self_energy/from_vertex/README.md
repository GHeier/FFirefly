# many_body / self_energy / from_vertex

## Overview

Computes the real-axis FLEX/spin-fluctuation self-energy Sigma(q,w) from a
previously-computed complex vertex V(k,w), by convolving V with the band
structure over the Brillouin zone, at any temperature (Temperature<=0 is the
exact T=0 limit).

## Quick Description

Computes the complex self-energy Sigma(q,w), the retarded Green's function G(q,w), and the spectral function A(q,w) on the real frequency axis, at any temperature, from the vertex V(k,w)

## Dependencies
- None

## Install Instructions

```bash
No special installation required - built automatically by fly-build.sh
```

### Parameters

outdir
prefix
filetype
fermi_energy - chemical potential mu used in the Fermi/Bose occupation factors
Temperature - physical temperature; Temperature<=0 (or unset, default 0.0) is the exact T=0 limit
q_mesh - output grid Sigma(q,w)/G(q,w) are evaluated and saved on
k_mesh - Brillouin zone mesh used for the internal k-sum (independent of q_mesh); V is evaluated at these points

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}_{prefix}_self_energy.{ext}` - Complex Sigma(q,w), saved on the same
  w-q grid convention as the input vertex/chi files
- `{outdir}_{prefix}_G.{ext}` - Complex retarded Green's function G(q,w) =
  1/(w - e(q) - Sigma(q,w)), same w-q grid convention
- `{outdir}_{prefix}_A_w.{ext}` - Real spectral function A(q,w) = -2 Im G(q,w),
  same w-q grid convention

## Testing

None

## Calculation Details

### Algorithm

The real-axis self-energy integral

    Sigma(q,w) = (1/N) sum_k integral (dv/pi) ImV(k,v) [n_B(v)+f(e(q-k))] / (w-v-e(q-k)+i*eta)

collapses analytically onto v = w - e(q-k) via the delta-function identity
Im[1/(x+i*eta)] -> -pi*delta(x) (no broadening or tetrahedron resolution
needed, since V(k,v) is evaluated directly via its own frequency
interpolation, at any temperature):

1. For each output q-point and frequency w, sum over k:
   `Im Sigma(q,w) = -(1/N) sum_k [n_B(w-e(q-k)) + f(e(q-k))] * ImV(k, w-e(q-k))`
   At `Temperature<=0`, `n_B(v) -> -Theta(-v)` and `f(e) -> Theta(mu-e)`
   exactly, reproducing the T=0 step-function formula bit-for-bit.
2. Terms where `w - e(q-k)` falls outside the tabulated frequency window of
   V(k,w) are dropped (Im V is assumed negligible there).
3. `n_B(v)` has a removable singularity at `v=0` (since `ImV(k,0)=0` for all
   k, by causality - Im V is odd in frequency); handled via the L'Hopital
   limit `-Temperature * d(ImV(k,v))/dv|_0` whenever `w-e(q-k)` lands within
   a small tolerance of 0.
4. `Re Sigma(q,w)` is obtained from `Im Sigma(q,w)` via
   `kk_transform_fft()` - an O(n log n) FFT-based Hilbert transform (the
   "analytic signal" sign-filter method, `Re Sigma = IFFT[i*sgn(xi)*FFT(Im Sigma)]`),
   requiring a uniformly-spaced w-grid and falling back to the exact O(n^2)
   `kk_transform()` otherwise. Validated against `kk_transform()` on
   synthetic data: agreement improves with grid resolution and converges to
   exact agreement (RMS difference and best-fit scale factor both -> 
   0%/1.0 respectively) by n~4096; expect a few-percent difference at the
   coarse w-grids (~10-20 points) typical of real configs, from ordinary
   discretization error in the FFT/trapezoidal quadratures, not a bug -
   `kk_transform()` is kept in the file as the correctness reference.
5. The retarded Green's function is built directly from the complex
   Sigma(q,w) via the Dyson equation:
   `G(q,w) = 1 / (w - e(q) - Sigma(q,w))`
   No extra `+i*eta` broadening is added - causality is already carried by
   `Im Sigma(q,w) <= 0`.
6. The spectral function `A(q,w) = -2 Im G(q,w)` is computed directly from G
   and saved alongside it. It's real and should be >= 0 everywhere (the same
   causality condition as step 5, `Im Sigma(q,w) <= 0`); a negative value
   printed as a diagnostic flags either a genuine causality violation
   upstream, or - if this file is later evaluated off its native mesh rather
   than at the saved grid points - a linear-interpolation artifact, since
   interpolating a near-singular complex function componentwise doesn't
   preserve its causal structure near a pole.

### Implementation Notes

- Band index is fixed to 0 (single/diagonal-band approximation); extend with
  a band loop for nbnd > 1 if interband contributions matter.
- Output q-points are a fresh mesh built from `q_mesh`, decoupled from the
  `k_mesh` used for the internal BZ sum - V(k,w) is evaluated via its own
  interpolation, so the output grid doesn't need to match whatever mesh V/chi
  were computed on. Falls back to V's own stored points only when V isn't on
  a regular mesh (as_mesh=false).
- Cost: `O(N_q*N_w*N_k)` for the main sum, `O(N_q*N_w*log(N_w))` for the
  Kramers-Kronig transform via `kk_transform_fft()`; parallelized over
  q-points via OpenMP.
- `kk_transform_fft()`'s FFT calls are wrapped in an OpenMP critical section
  (`fftw_planner`), since it's called once per q-point inside the
  parallelized q-loop and `fft_1d_f()` creates/destroys an FFTW plan on
  every call - FFTW's planner is not thread-safe by default in this repo's
  build (no `fftw3f_threads`/`fftw3f_omp` linked, which is what would let
  `fftwf_make_planner_thread_safe()` be used instead), and calling it
  concurrently from multiple threads without this corrupts FFTW's internal
  state (observed directly during implementation: a `planner.c` assertion
  failure followed by a segfault). Only plan creation/execution is
  serialized, not the surrounding O(n log n) computation, so this doesn't
  undermine the complexity win.
- An FFT-accelerated variant, `call_flex_self_energy_fft()`, is also
  available (selected via a flag in `run.cpp`) - it uses a momentum-space
  cyclic convolution instead of the direct sum, trading exact analytic
  delta-function resolution for a mesh-adaptive Lorentzian broadening. It
  requires `k_mesh == q_mesh` and is only faster when `N_w` is small relative
  to `N_k` (its cost scales as `O(N_w^2 * N_k log N_k)`, unaffected by this
  method's own O(n log n) Kramers-Kronig improvement, since it never used
  `kk_transform()` in the first place - Re/Im Sigma come out together from
  the same Lorentzian convolution); for a densely sampled vertex (hundreds
  of w-points), the direct method is usually faster.
- Only saves as hdf5 file at the moment.

## References

1. Bickers, Scalapino, White, "Conserving Approximations for Strongly Correlated Electron Systems: Bethe-Salpeter Equation and Dynamics for the Two-Dimensional Hubbard Model", Phys. Rev. Lett. 62, 961 (1989).
2. Scalapino et al., "d-wave pairing near a spin-density-wave instability", Phys. Rev. B. 34, (1986).
