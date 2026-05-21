# many_body / renormalization / FS_approx

## Overview

Calculates the quasiparticle renormalization factor Z and mass enhancement λ_z using Fermi surface averaging. This approach uses the static (ω=0) approximation for the interaction vertex, computing the average coupling strength and effective mass across the Fermi surface. Useful for estimating many-body corrections to superconducting properties.

## Quick Description

Computes quasiparticle weight Z(k) across the Fermi Surface using V(w)=V(0) approximation.

## Dependencies

- NumPy for array operations
- Firefly Field classes for loading susceptibility and renormalization data
- Bands and Surface classes for Fermi surface generation
- HDF5 for I/O

## Install Instructions

```bash
pip install numpy h5py
```

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dimension` | int | 3 | Spatial dimensionality |
| `fermi_energy` | float | 0.0 | Chemical potential μ |
| `mu_from_n` | bool | false | Compute μ from electron number |
| `num_electrons` | float | - | Target electron count if mu_from_n |
| `U0` | float | - | Hubbard interaction strength |
| `qp_weight` | float | 1.0 | Input quasiparticle weight Z |

## Results Saved

- `{outdir}/{prefix}_Z_FS.h5` - Fermi surface averaged renormalization
- Console output: average m*, λ_z = 1/Z - 1, total DOS

## Testing

Expected test behavior:
- For U=0, should return Z=1 (no renormalization)
- λ_z should increase with U and susceptibility enhancement
- Should match analytic results for simple models

## Calculation Details

### Algorithm

1. If mu_from_n enabled: load E_vs_n curve and interpolate to find μ for target n
2. Generate Fermi surface at μ using band structure
3. Compute FS weights: w(k) = dA_k / |v_F(k)| / (2π)^d
4. Load susceptibility χ(q) and compute vertex: V(q) = U² × Z × χ(q)
5. Construct scattering matrix V(k-k') for all FS k-point pairs
6. Compute average coupling: λ = Σ_{k,k'} w(k) × V(k-k') × w(k')
7. Load Z(k) field and compute FS average: Z_avg = Σ_k w(k) × Z(k)
8. Output effective mass m* = 1 + λ_z where λ_z = 1/Z - 1

### Implementation Notes

- Static approximation V(ω=0) valid for low-frequency properties
- Fermi surface averaging appropriate for transport and superconductivity
- Can optionally determine μ self-consistently from electron number
- Coupling λ related to mass enhancement and T_c

## References

1. D. J. Scalapino, "A common thread: The pairing interaction for unconventional superconductors", Rev. Mod. Phys. 84, 1383 (2012).
2. G. D. Mahan, "Many-Particle Physics", 3rd ed., Springer (2000).
