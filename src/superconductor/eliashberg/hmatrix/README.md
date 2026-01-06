# superconductor / eliashberg / hmatrix

## Overview

Solves Eliashberg equation on real axis using HMatrix compression and Lanczos solver.

## Quick Description

Solves Eliashberg equation on real axis using HMatrix compression and Lanczos solver.

## Dependencies
- HMatrices
- KrylovKit
- LinearAlgebra
- Printf
- StaticArrays

## Install Instructions

Download packages from Julia REPL

### Parameters

Configuration parameters from `input.cfg`:

k_mesh - used implicitly in c++ call of Surface creation (uses tetrahedra)
w_pts - number of real axis frequency points
cutoff_energy - cutoff energy for real frequencies

## Results Saved

- `{outdir}_{prefix}_gap.{ext}` - Dominant gap function on real axis

## Testing

With constant vertex (V(k,k')=-1), see if it converges to Tc = 4 (derived analytically).

## Calculation Details

### Algorithm

Description of the algorithm:
1. Define Fermi Surface
2. Load Vertex from file
3. Construct HMatrix representation of interaction kernel
4. Use Lanczos solver to solve Eliashberg equation iteratively
5. Iterate until Tc converges and lambda=1

### Implementation Notes

- HMatrix improves memory and computational efficiency
- Lanczos solver accelerates convergence by finding dominant eigenvalues

## References

