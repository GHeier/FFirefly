# many_body / self_energy / triqs

## Overview

Calculates the self-energy using Iterated Perturbation Theory (IPT) using the TRIQS library. Calculations and output are on imaginary axis, and it depends on the Density of States. Outputs Self-Energy and Green's function on imaginary axis, as well as spectral function on real axis (found using Pade approximants on Green's function).

## Quick Description

Calculates the self-energy using Iterated Perturbation Theory (IPT) on the imaginary axis.

## Dependencies
- Density of States calculation is taken as input
- TRIQS library

## Install Instructions

```bash
Install TRIQS and tprf libraries as per their documentation.
- https://triqs.github.io/triqs/latest/install.html#packaged-versions-of-triqs
- https://triqs.github.io/tprf/latest/install.html#
```

### Parameters

w_points: Number of Matsubara frequency points to use.
Temperature: Temperature for the calculation.
U: On-site interaction strength.
mu: Chemical potential.
mixing: Mixing parameter for self-energy convergence.
max_iters: Maximum number of iterations for self-consistency.

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir_}{prefix}_G_iw.h5` - Local Green's function on imaginary axis.
- `{outdir_}{prefix}_Sigma_iw.h5` - Self-energy on imaginary axis.
- `{outdir_}{prefix}_A_w.h5` - Spectral function (new DOS) on real axis.

File format details:
- All files are f(w) vs w datasets in HDF5 format.

## Testing

- No reference tests created, but it qualitatively matches known IPT results.

## Calculation Details

Performs $\Sigma(i\omega_n) = \int d\tau U^2 G(\tau)^3$ where $G$ is the interacting Green's function calculated self-consistently.

### Algorithm

Description of the algorithm:
1. Using DOS, calculate HilbertTransform for fast energy integrals.
2. Calculate local non-interacting Green's function
3. Calculate self-energy using IPT formula
4. Update Green's function using Dyson equation
5. Loop steps 3-4 until convergence or max_iters reached
6. Use Pade approximants to analytically continue Green's function to real axis and compute spectral function.

### Implementation Notes

- Only works for 1-band at the moment
- Playing with mixing can be important for convergence

## References

1. https://spm-lab.github.io/sparse-ir-tutorial/src/DMFT_IPT_py.html
