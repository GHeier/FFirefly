# superconductor / bcs / matrix

## Overview

Solves the linearized BCS gap equation using direct matrix diagonalization on the Fermi surface. Computes the superconducting gap function Δ(k) by finding the largest eigenvalue of the pairing kernel matrix. This approach is exact but limited to relatively small Fermi surface discretizations due to O(N²) memory scaling.

## Quick Description

Solves the linearized BCS gap equation on the Fermi Surface using standard matrix diagonalization.

## Dependencies

- Band structure and Fermi surface calculation
- Vertex/interaction field (loaded from HDF5)
- LAPACK/OpenBLAS for eigenvalue solver
- HDF5 for I/O

## Install Instructions

```bash
# No special installation required - built automatically by fly-build.sh
```

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fermi_energy` | float | 0.0 | Chemical potential μ in eV |
| `Temperature` | float | 0.01 | Temperature T in eV |
| `k_mesh` | int[3] | - | k-point mesh for Fermi surface |
| `cutoff_energy` | float | - | Interaction cutoff ω_c in eV |
| `qp_weight` | float | 1.0 | Quasiparticle weight Z |
| `FS_only` | bool | true | Restrict to Fermi surface |
| `num_eigenvalues_to_save` | int | 1 | Number of gap solutions to compute |
| `U0` | float | - | Interaction strength U |

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}/{prefix}_gap.dat` - Gap function Δ(k) at Fermi surface k-points
- `{outdir}/{prefix}_gap.h5` - Gap function in HDF5 format (if enabled)

File format details:
- Text format: columns are kx, ky, kz, Re(Δ), Im(Δ)
- HDF5: complex scalar field on Fermi surface points

## Testing

Expected test behavior:
- For constant attractive interaction V < 0, should recover BCS result: Δ = ω_c * exp(-1/|V|N(0))
- Eigenvalue λ > 1 indicates superconducting instability
- s-wave gap should be constant on Fermi surface for isotropic interaction

## Calculation Details

### Algorithm

1. Generate Fermi surface at chemical potential μ
2. Compute DOS from FS area elements and Fermi velocities: N(k) = dA_k / |v_F(k)|
3. Load interaction vertex V(q) from file or compute from U
4. Construct pairing matrix: P(k,k') = -√(N(k)N(k')) × V(k-k') / (2π)^d
5. Apply frequency cutoff via form factor: f(k) = tanh(ε_k / 2T)
6. Solve eigenvalue problem: P × Δ = λ × Δ using power iteration
7. Extract critical temperature from eigenvalue: T_c ∝ ω_c × exp(-1/λ)
8. Apply quasiparticle renormalization: Δ → Z × Δ

### Implementation Notes

- Memory scales as O(N_FS²) where N_FS is number of Fermi surface points
- Power iteration finds dominant eigenvalue efficiently
- Multiple eigenvalues reveal competing pairing symmetries (s, d, p-wave)
- Quasiparticle weight Z accounts for many-body renormalization

## References

1. J. Bardeen, L. N. Cooper, and J. R. Schrieffer, "Theory of Superconductivity", Phys. Rev. 108, 1175 (1957).
2. D. J. Scalapino, "A common thread: The pairing interaction for unconventional superconductors", Rev. Mod. Phys. 84, 1383 (2012).
