# many_body / response / tetrahedra

## Overview

Calculates chi0(w,q) using recursive tetrahedron method for BZ integration. Takes cutoff_energy and w_pts to determine frequency grid, and q_mesh to determine momentum grid.

## Quick Description

Calculates non-interacting response function chi0(w,q) using recursive tetrahedron integration

## Dependencies
- BZIntegral
- LinearAlgebra
- Printf
- Interpolations
- Base.Threads

## Install Instructions

```bash
Install all in Julia REPL
```

### Parameters

kmesh - input for e(k) creation
qmesh - grid of q-points for chi0(q,w) calculation
cutoff_energy - energy cutoff for frequency grid
w_pts - number of frequency points
dim - system dimension (2 or 3)
nbnds - number of bands to consider
prefix - prefix for output files
outdir - output directory
fermi_energy - Fermi energy
brillouin_zone - BZ matrix


## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}_{prefix}_chi.h5` - chi0(w,q) on w-k grid. 

## Testing

For 3D TB, mu=0.0, t'=0.0, T=0.25, max_chi~=0.4

## Calculation Details

### Algorithm

Description of the algorithm:
1. Tetrahedrizes the Brillouin Zone
2. Recursively subdivides each tetrahedra
3. Integrates over recursive weights and maps back to k-grid
4. Sums over weights to get value at w,q point
4. Repeats for every w-q

### Implementation Notes

- k and q meshes have to be odd to avoid the Gamma point divergence, so the code enforces this, adding 1 if even values are given.
- Parallelized over q-points using Julia's `Threads.@threads`
- Uses `Interpolations.jl` for fast H(k) interpolation used in tetrahedron integration

## References

1. https://github.com/SelimLin/BZIntegral.jl
