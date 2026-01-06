# many_body / many_body / triqs

## Overview

Calculates FLEX or FLEX+DMFT, taking the electron density as a conserved quantity, calculating the new chemical potential, vertex, response, green's function, and self-energies at every iteration in the self-consistent loop. FLEX+DMFT substitutes the local part of the FLEX self-energy with the DMFT self-energy calculated using the IPT solver.

## Quick Description

Solves FLEX or FLEX+DMFT self-consistently with DLR calculations using matsubara frequencies at finite Temperature

## Dependencies
- triqs
- triqs_tprf

## Install Instructions

```bash
Check TRIQS install instructions at https://triqs.github.io/triqs/latest/install.html#
Check TRIQS/TPRF install instructions at https://triqs.github.io/tprf/latest/install.html#
```

### Parameters

interaction - vertex type, options are "FLEX" or "FLEX+DMFT"
num_electrons - used as a conserved quantity to determine the chemical potential
fermi_level - initial guess for the chemical potential
Temperature - chosen temperature for calculations
mixing - mixing parameter for self-consistency
max_iters 
outdir
prefix
onsite_U 

## Results Saved

- `{outdir}_{prefix}_vertex.{ext}` - Vertex on w-k grid
- `{outdir}_{prefix}_chi.{ext}` - Susceptibility on w-k grid
- `{outdir}_{prefix}_sigma.{ext}` - Self-energy on w-k grid
- `{outdir}_{prefix}_G.{ext}` - Green's function on w-k grid
- `{outdir}_{prefix}_G0.{ext}` - Non-interacting Green's function on w-k grid

## Testing

Returns maximum chi value, for 3D TB, t'=0, mu=0, T=0.25, max_chi~0.4

## Calculation Details

### Algorithm

1. Calculate X(w,k) as integral over G(w,k) * G(-w,-k)
2. Calculate V(w,k) using FLEX formula
3. Calculate Sigma(w,k) as integral over G(w,k) * V(w,k)
4. Calculate new G(w,k) using Dyson equation
5. Loop back to step 1 until consistency

### Implementation Notes

- Error is determined by difference in G(w,k) between iterations

## References

1. https://spm-lab.github.io/sparse-ir-tutorial/src/FLEX_jl.html
