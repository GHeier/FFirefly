# superconductor / eliashberg / sparse_ir

## Overview

Solves the linearized Eliashberg equation on the imaginary axis using the Sparse Intermediate Representation (sparse-ir) for efficient Matsubara frequency handling. This approach captures strong-coupling effects and retardation while maintaining computational efficiency through the compact IR basis representation.

## Quick Description

Solves the linearized Eliashberg equation on imaginary axis using convolution and the power iteration / Krylov projection approach.

## Dependencies

- SparseIR.jl for intermediate representation basis
- LinearAlgebra for matrix operations
- FFTW for convolutions
- Firefly Field classes for data I/O

## Install Instructions

```julia
using Pkg
Pkg.add("SparseIR")
```

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | - | k-point mesh dimensions |
| `Temperature` | float | 0.01 | Temperature T in eV |
| `fermi_energy` | float | 0.0 | Chemical potential μ |
| `w_max` | float | 10.0 | Frequency cutoff Λ for IR basis |
| `ir_tol` | float | 1e-10 | IR basis truncation tolerance |

## Results Saved

- `{outdir}/{prefix}_gap.h5` - Gap function Δ(k,iω_n)
- `{outdir}/{prefix}_gap_ir.h5` - Gap in IR basis coefficients
- `{outdir}/{prefix}_eigenvalue.dat` - Leading eigenvalue

## Testing

Expected test behavior:
- Should reproduce BCS limit for weak coupling
- Strong coupling should show enhanced T_c and frequency-dependent gap
- IR basis should converge rapidly with basis size

## Calculation Details

### Algorithm

1. Construct sparse-ir basis for given temperature and frequency cutoff
2. Load vertex V(k,iω) and transform to IR basis
3. Compute anomalous self-energy in IR: Σ_IR = V_IR ⊗ G_IR
4. Transform kernel to act on gap function: K[Δ] = -T Σ_k' V(k-k') × F(k',iω)
5. Use power iteration or Lanczos to find largest eigenvalue
6. Transform solution back to Matsubara frequencies
7. Save gap function in both frequency and IR representations

### Implementation Notes

- IR basis provides exponentially compact representation of Matsubara functions
- Typical calculations need only ~50 IR basis functions vs thousands of frequencies
- Convolutions performed efficiently using FFT in real space
- Captures both weak-coupling BCS and strong-coupling Eliashberg limits

## References

1. H. Shinaoka et al., "Compressing Green's function using intermediate representation between imaginary-time and real-frequency domains", Phys. Rev. B 96, 035147 (2017).
2. M. Wallerberger et al., "sparse-ir: Optimal compression and sparse sampling of many-body propagators", SoftwareX 21, 101266 (2023).
