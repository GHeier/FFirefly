# superconductor / eliashberg / power_iteration

## Overview

Solves the linearized Eliashberg equation using power iteration to find the dominant superconducting gap eigenvalue and eigenvector on Fermi surface patches.

## Quick Description

Iteratively applies the Eliashberg kernel to find the leading eigenvalue λ and gap function Δ(k,iω), using Fermi surface patching and sparse_ir for efficient frequency representation.

## Dependencies

### Required
- SparseIR (Julia library)
- FastGaussQuadrature (Julia)
- Firefly (for band structure and vertex interaction)

### Optional
- MPI (for parallelization, currently unused)

## Install Instructions

```bash
julia -e 'using Pkg; Pkg.add(["SparseIR", "FastGaussQuadrature", "FFTW"])'
```

## Results Saved

Output files created by this calculation (using `prefix` from config):

- Returns maximum gap value `max(|Δ|)` as convergence indicator

File format:
- Currently returns scalar value for testing
- Can be extended to save full gap function Δ(k,iω) as HDF5

## Testing

Run the test suite:
```bash
fly.x  # Runs all tests including this one
```

Expected test behavior:
- Constructs Fermi surface patches using Gaussian quadrature
- Loads vertex interaction V(k,iν)
- Iterates until convergence or max iterations
- Returns max gap magnitude
- Typical values: |Δ| ~ 0.001-0.1 depending on interaction strength

## Calculation Details

### Algorithm

1. Construct Fermi surface patches at energies ε = ωc × (Gauss-Legendre points)
2. Load band structure and pairing vertex V from Firefly
3. Create IR mesh for frequency sampling
4. Initialize gap Δ, renormalization Z, and shift χ on each patch
5. Power iteration loop:
   - Compute F(ε,iω) = Δ / [(Z·iω)² + Δ² + (ε - μ + χ)²]
   - Compute G(ε,iω) = -(Z·iω + ε + χ - μ) / [-(Z·iω)² + Δ² + (ε + χ - μ)²]
   - Convolve with vertex: Δ_new = Σ V(k-k',iω-iω') F(k',iω')
   - Convolve for self-energy: Σ = Σ V(k-k',iω-iω') G(k',iω')
   - Extract Z and χ from Σ
   - Check convergence: |Δ_new - Δ| < tol
6. Return maximum gap value

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | [8,8,8] | k-space mesh (for future extensions) |
| `nbnd` | int | 1 | Number of bands |
| `dimension` | int | 3 | Spatial dimension (2 or 3) |
| `Temperature` | float | 0.01 | Temperature in energy units (β = 1/T) |
| `fermi_energy` | float | 0.0 | Chemical potential μ |
| `cutoff_energy` | float | 1.0 | Energy cutoff ωc for Fermi surface integration |
| `onsite_U` | float | 1.0 | Interaction strength (if using model vertex) |
| `brillouin_zone` | matrix | I | BZ matrix for coordinate transformation |
| `prefix` | string | "sample" | Prefix for output files |
| `outdir` | string | "./" | Output directory |

### Implementation Notes

- Uses Fermi surface patching to reduce dimensionality from 3D k-space to 2D surfaces
- Gaussian-Legendre quadrature provides optimal energy sampling points
- Sparse_ir efficiently handles Matsubara frequency sums via IR basis
- Power iteration finds only the dominant eigenvalue (fastest method for single eigenvalue)
- Convolution performed in (k,ω) space using element-wise operations
- Currently uses 7 Fermi surface patches (hardcoded in `num_surfaces`)
- Self-consistency in Z and χ included for full Eliashberg treatment

## References

1. Eliashberg, G. M., "Interactions between electrons and lattice vibrations in a superconductor", Sov. Phys. JETP 11, 696 (1960)
2. Scalapino, D. J., "A common thread: The pairing interaction for unconventional superconductors", Rev. Mod. Phys. 84, 1383 (2012)
3. SparseIR documentation: https://spm-lab.github.io/sparse-ir-tutorial/
