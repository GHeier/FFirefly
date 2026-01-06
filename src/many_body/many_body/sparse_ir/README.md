# many_body / many_body / sparse_ir

## Overview

Performs FLEX calculations using DLR sparse_ir code, with the option of self-consistency

## Quick Description

Performs FLEX calculations using DLR sparse_ir code, with the option of self-consistency

## Dependencies
- SparseIR julia library

## Install Instructions

```bash
Install SparseIR for julia
```

### Parameters

k_mesh - handles input mesh for H(k). Output is the same
dimension - 2 and 3D calculations allowed
Temperature - Temperature calculations are performed at (matsubara frequencies)
max_iters - upper limit for self-consistent loop
outdir
prefix
filetype
self_consistent
nbnd
fermi_energy
onsite_U
brillouin_zone

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}_{prefix}_self_energy.{ext}` - self-energy on w-k grid
- `{outdir}_{prefix}_chi.{ext}` - response on w-k grid
- `{outdir}_{prefix}_vertex.{ext}` - vertex on w-k grid

## Testing

Returns maximum value of chi. For 3D TB t'=0 at mu=0 T=0.25, max_chi~0.4

## Calculation Details

### Algorithm

1. Calculate X(w,k) as integral over G(w,k) * G(-w,-k)
2. Calculate V(w,k) using FLEX formula
3. Calculate Sigma(w,k) as integral over G(w,k) * V(w,k)
4. Calculate new G(w,k) using Dyson equation
5. Loop back to step 1 until consistency

### Implementation Notes

- Error is determined by difference in Sigma(w,k) between iterations

## References

1. https://spm-lab.github.io/sparse-ir-tutorial/src/FLEX_jl.html
