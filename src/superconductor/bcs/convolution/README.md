# superconductor / bcs / convolution

## Overview

Solves the linearized BCS gap equation using FFT-based convolution in real space. This approach exploits momentum-space translation symmetry to achieve O(N log N) scaling, making it efficient for large k-meshes where direct matrix methods become prohibitive. Uses Lanczos iteration to find multiple eigenvalues and identify competing pairing symmetries.

## Quick Description

Solves the linearized BCS gap equation across the Brillouin Zone 

## Dependencies

- NumPy and SciPy for FFT and sparse eigensolvers
- ARPACK (via scipy.sparse.linalg) for Lanczos iteration
- Firefly Field classes for loading vertex/self-energy
- HDF5 for I/O

## Install Instructions

```bash
pip install numpy scipy h5py
```

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | - | k-point mesh dimensions [nx, ny, nz] |
| `dimension` | int | 3 | Spatial dimensionality |
| `Temperature` | float | 0.01 | Temperature T in eV |
| `fermi_energy` | float | 0.0 | Chemical potential μ in eV |
| `num_solutions` | int | 1 | Number of eigenvalues to find |
| `qp_weight` | float | 1.0 | Quasiparticle weight Z |

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}/{prefix}_gap.h5` - Dominant gap function Δ(k)
- `{outdir}/{prefix}_gap_eig{i}.h5` - Additional gap solutions if num_solutions > 1

File format details:
- Complex scalar fields on k-mesh in HDF5 format
- Eigenvalues stored as metadata

## Testing

Expected test behavior:
- For attractive Hubbard model, dominant eigenvalue should correspond to s-wave pairing
- For repulsive model near half-filling, d-wave should dominate
- Eigenvalue ordering reveals pairing symmetry hierarchy

## Calculation Details

### Algorithm

1. Load vertex V(k) and self-energy Σ(k) from files
2. Compute BCS form factor: f(k) = tanh(β × ε(k) / Z) / (2 × ε(k))
3. Transform vertex to real space: V(r) = FFT[V(k)]
4. Define kernel action: K[Δ](k) = IFFT[V(r) × FFT[f(k) × Δ(k)]]
5. Use Lanczos/ARPACK to find largest eigenvalues of K
6. Scale eigenvalues by Z and mesh normalization
7. Save gap functions for each eigenvalue

### Implementation Notes

- FFT convolution reduces O(N²) matrix-vector product to O(N log N)
- Lanczos iteration finds extremal eigenvalues without forming full matrix
- Multiple eigenvalues identify s-wave, d-wave, p-wave, etc.
- Memory efficient: only stores O(N) vectors, not O(N²) matrix

## References

1. D. J. Scalapino, E. Loh, and J. E. Hirsch, "d-wave pairing near a spin-density-wave instability", Phys. Rev. B 34, 8190 (1986).
2. T. A. Maier et al., "Systematic study of d-wave superconductivity in the 2D repulsive Hubbard model", Phys. Rev. Lett. 95, 237001 (2005).
