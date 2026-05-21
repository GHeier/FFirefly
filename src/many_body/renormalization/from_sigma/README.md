# many_body / renormalization / from_sigma

## Overview

Extracts the quasiparticle renormalization factor Z(k) from the frequency derivative of the self-energy at ω→0. This is the standard definition Z = 1 - ∂Σ/∂ω|_{ω=0}, computed numerically from self-energy data on a Matsubara frequency grid. Outputs both the k-resolved Z(k) and Fermi surface averaged values.

## Quick Description

Calculates Z(k) based on the slope of a given Sigma(iω,k) at ω→0

## Dependencies

- Self-energy field from previous calculation (self_energy/triqs or self_energy/sparse_ir)
- Band structure for Fermi surface generation
- HDF5 for I/O

## Install Instructions

```bash
# No special installation required - built automatically by fly-build.sh
```

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | - | k-point mesh dimensions |
| `q_mesh` | int[3] | - | Output mesh (if FS data on mesh) |
| `dimension` | int | 3 | Spatial dimensionality |
| `fermi_energy` | float | 0.0 | Fermi level for FS averaging |
| `filetype` | string | h5 | Output file format |

## Results Saved

- `{outdir}/{prefix}_renormalization.{filetype}` - Z(k) field matching self-energy grid

File format details:
- Complex field with same structure as input self-energy
- If input is on mesh: output is Field on mesh
- If input is scattered points: output matches those points

## Testing

Expected test behavior:
- For non-interacting system (Σ=0), should return Z=1 everywhere
- Z < 1 indicates quasiparticle weight reduction
- Should be consistent with analytic/FS_approx methods

## Calculation Details

### Algorithm

1. Load self-energy Σ(k,iω_n) from HDF5 file
2. Identify lowest Matsubara frequencies for derivative
3. Compute numerical derivative: ∂Σ/∂ω ≈ [Σ(iω_1) - Σ(iω_0)] / (ω_1 - ω_0)
4. Extract imaginary part (gives slope): slope = -Im(∂Σ/∂ω)
5. Compute Z(k) = 1 / (1 - slope) at each k-point
6. Track maximum and average renormalization
7. If Fermi surface available: compute FS-weighted average using dA/|v_F|
8. Save Z(k) field with appropriate mesh/point structure

### Implementation Notes

- Uses finite difference for frequency derivative
- Requires self-energy on fine enough frequency grid near ω=0
- Z < 1 indicates correlation-induced mass enhancement: m*/m = 1/Z
- Fermi surface average gives transport-relevant effective mass

## References

1. G. D. Mahan, "Many-Particle Physics", 3rd ed., Springer (2000), Chapter 5.
2. A. A. Abrikosov, L. P. Gorkov, and I. E. Dzyaloshinski, "Methods of Quantum Field Theory in Statistical Physics", Dover (1975).
