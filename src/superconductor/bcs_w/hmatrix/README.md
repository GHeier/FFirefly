# superconductor / bcs+w / hmatrix

## Overview

Solves the frequency-dependent BCS gap equation using Hierarchical Matrix (H-Matrix) compression combined with Lanczos iteration. This method includes Matsubara frequency dependence in the pairing kernel while maintaining computational efficiency through low-rank approximations of the interaction matrix blocks.

## Quick Description

Solves the linearized BCS gap equation across the Fermi Surface using compressed Hierarchical Matrices and a lanczos matrix solver.

## Dependencies

- HMatrices.jl for hierarchical matrix operations
- LinearAlgebra for eigenvalue computations
- Firefly Field classes for data I/O
- HDF5 for file storage

## Install Instructions

```julia
using Pkg
Pkg.add("HMatrices")
```

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | - | k-point mesh dimensions |
| `w_pts` | int | 100 | Number of Matsubara frequencies |
| `Temperature` | float | 0.01 | Temperature T in eV |
| `fermi_energy` | float | 0.0 | Chemical potential μ |
| `compression_tol` | float | 1e-6 | H-matrix compression tolerance |

## Results Saved

- `{outdir}/{prefix}_gap.h5` - Gap function Δ(k,iω_n)
- `{outdir}/{prefix}_eigenvalue.dat` - Leading eigenvalue and T_c estimate

## Testing

Expected test behavior:
- For constant vertex V(k,iω) = V₀, should recover standard BCS result
- Eigenvalue should converge with increasing w_pts
- H-matrix compression should achieve significant memory savings for large systems

## Calculation Details

### Algorithm

1. Generate Fermi surface points and Matsubara frequency grid
2. Construct pairing kernel K(k,iω; k',iω') including frequency dependence
3. Compress kernel into H-matrix format using adaptive cross approximation
4. Apply Lanczos iteration to find largest eigenvalue of compressed kernel
5. Extract gap function from dominant eigenvector
6. Estimate T_c from eigenvalue crossing λ = 1

### Implementation Notes

- H-matrix compression reduces O(N²) storage to O(N log N)
- Particularly effective when interaction is smooth in momentum space
- Compression tolerance controls accuracy vs memory tradeoff
- Frequency dependence important for retardation effects

## References

1. W. Hackbusch, "Hierarchical Matrices: Algorithms and Analysis", Springer (2015).
2. S. Borm, L. Grasedyck, and W. Hackbusch, "Introduction to hierarchical matrices with applications", Eng. Anal. Bound. Elem. 27, 405 (2003).
