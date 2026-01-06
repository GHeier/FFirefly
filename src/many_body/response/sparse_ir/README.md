# many_body / response / sparse_ir

## Overview

Calculates bare susceptibility χ(iν,k) from band structure using sparse_ir for efficient frequency representation.

## Quick Description

Computes the non-interacting response function by constructing the Green's function from the band structure and performing the convolution χ(r,τ) = G(r,τ) · G(r,-τ).

## Dependencies
- SparseIR julia library
- FFTW

## Install Instructions

```bash
julia -Pkg; Pkg.add("SparseIR"); Pkg.add("FFTW")
```

## Parameters

Configuration parameters from `input.cfg`:

- `k_mesh` - k-space grid for calculation
- `dimension` - 2D or 3D calculation
- `Temperature` - Sets inverse temperature β = 1/T for Matsubara frequencies
- `fermi_energy` - Chemical potential μ
- `brillouin_zone` - BZ matrix for coordinate transformation
- `outdir` - Output directory
- `prefix` - File prefix for output
- `filetype` - Output format (h5, dat, etc.)

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{prefix}_chi.{ext}` - Bare susceptibility χ(iν,k) on bosonic frequency and momentum grid

## Testing

Returns maximum value of χ. For 3D TB at μ=0, T=0.25, typical max(χ) ~ 0.4

## Calculation Details

### Algorithm

1. Load band structure ε(k) via `Bands()`
2. Construct bare Green's function: G(iω,k) = 1/(iω - ε(k) + μ)
3. Transform to real space and imaginary time: G(r,τ)
4. Compute convolution: χ(r,τ) = G(r,τ) · G(r,-τ)
5. Transform back to momentum and bosonic frequencies: χ(iν,k)
6. Save result

### Implementation Notes

- Uses intermediate representation (IR) for efficient τ ↔ ω transforms
- Convolution performed in (r,τ) space for efficiency
- Output k-points centered at (-π,π) instead of (0,2π)

## References

1. https://spm-lab.github.io/sparse-ir-tutorial/
