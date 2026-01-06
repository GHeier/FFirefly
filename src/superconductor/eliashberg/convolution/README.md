# superconductor / eliashberg / convolution

## Overview

Solves the linearized Eliashberg equation using Lanczos algorithm to find multiple superconducting gap eigenvalues and eigenvectors, useful for identifying competing pairing channels.

## Quick Description

Uses ARPACK's Lanczos eigensolver to find the leading eigenvalues of the Eliashberg kernel K, where Δ = λK[Δ], returning multiple eigenpairs to identify dominant and subdominant pairing symmetries.

## Dependencies

### Required
- TRIQS (The Toolbox for Research on Interacting Quantum Systems)
- NumPy
- SciPy (sparse.linalg.eigsh)
- HDF5

### Optional
- None

## Install Instructions

```bash
# Install TRIQS following instructions at triqs.github.io
# Python dependencies installed automatically via pip
```

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{prefix}_gap.h5` - Superconducting gap function Δ(k,iω) for the leading eigenvalue (HDF5 format, Field_CM object)

File format:
- HDF5 file with mesh structure compatible with Firefly Field_CM
- Contains gap function on k-mesh × DLR frequency mesh
- Matrix-valued for multi-band systems (shape: [nω, nk, nbnd, nbnd])

## Testing

Run the test suite:
```bash
fly.x  # Runs all tests including this one
```

Expected test behavior:
- Loads pre-computed G(k,iω) and V(k,iν) from HDF5 files
- Solves for leading eigenvalues
- Returns maximum eigenvalue λ_max
- λ > 1 indicates superconducting instability

## Calculation Details

### Algorithm

1. Load Green's function G(k,iω) from `{prefix}_G.h5`
2. Load pairing interaction vertex V(k,iν) from `{prefix}_vertex.h5`
3. Transform V to real space and imaginary time: V(r,τ)
4. Define Eliashberg kernel: K[Δ] = -∫ V(k-k',iω-iω') G(k',iω') G(-k',-iω') Δ(k',iω') dk'dω'
5. Set up LinearOperator for matrix-free Lanczos iteration
6. Call scipy.sparse.linalg.eigsh to find top `num_solutions` eigenvalues
7. Save gap function for leading eigenvalue

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | [8,8,8] | k-space mesh for Green's function |
| `nstates` | int | 1 | Number of bands/orbitals |
| `Temperature` | float | 0.01 | Temperature in energy units (sets β = 1/T) |
| `fermi_energy` | float | 0.0 | Chemical potential μ |
| `num_solutions` | int | 5 | Number of eigenvalues to compute |
| `prefix` | string | "sample" | Prefix for input/output files |
| `outdir` | string | "./" | Output directory |

### Implementation Notes

- Uses TRIQS DLR (Discrete Lehmann Representation) for efficient frequency sampling
- Matrix-free implementation: only stores G, V, Δ (not full Eliashberg matrix)
- Lanczos converges to extreme eigenvalues, ideal for finding instabilities
- Returns multiple eigenpairs to identify competing superconducting channels
- Convolution performed using Diagram class for efficient k,ω sums

## References

1. Eliashberg, G. M., "Interactions between electrons and lattice vibrations in a superconductor", Sov. Phys. JETP 11, 696 (1960)
2. Scalapino, D. J., "A common thread: The pairing interaction for unconventional superconductors", Rev. Mod. Phys. 84, 1383 (2012)
3. TRIQS documentation: https://triqs.github.io/
