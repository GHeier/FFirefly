# hamiltonian / generate / hk_from_hr

## Overview

Generates the Hamiltonian matrix H(k) on a k-point mesh from tight-binding parameters. Uses Fourier transformation from real-space hopping integrals H(R) to momentum space, producing a complex matrix field that can be used by other calculations for band structure, Fermi surface, and many-body physics.

## Quick Description

Computes Hamiltonian, H(k) based on an H(r) tight-binding construction

## Dependencies

- tbmodels library for tight-binding Fourier transforms
- NumPy for array operations
- HDF5 for output

## Install Instructions

```bash
pip install tbmodels numpy h5py
```

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | - | k-point mesh dimensions [nx, ny, nz] |
| `dimension` | int | 3 | Spatial dimensionality (1, 2, or 3) |
| `hamiltonian` | string | tight_binding | Model type |
| `t0` | float | 1.0 | On-site energy / nearest-neighbor hopping |
| `t1`...`t10` | float | 0.0 | Higher-neighbor hopping parameters |
| `nstates` | int | 1 | Number of orbitals/states per unit cell |

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}/{prefix}_Hk.h5` - Hamiltonian matrix field H(k) in HDF5 format

File format details:
- Complex matrix field (Field_CM) with shape [nx, ny, nz, nstates, nstates]
- Stored with mesh and domain metadata for interpolation
- Can be loaded with Firefly's Field_CM class

## Testing

Expected test behavior:
- For single-band tight-binding, H(k) should match analytical ε(k) = -2t(cos(kx) + cos(ky) + cos(kz))
- Eigenvalues at high-symmetry points should match expected band energies

## Calculation Details

### Algorithm

1. Initialize tight-binding model with on-site energies
2. Add hopping terms for each specified neighbor shell (t0, t1, ..., t10)
3. Define hopping vectors for each direction based on dimensionality
4. Evaluate H(k) on full k-point mesh using FFT from real-space representation
5. Reshape result into [nx, ny, nz, nstates, nstates] tensor
6. Save to HDF5 as Field_CM with mesh and domain metadata

### Implementation Notes

- Uses tbmodels for efficient FFT-based evaluation
- Supports multi-orbital models with nstates > 1
- Hopping parameters follow standard tight-binding conventions
- Output can be used directly by Hamiltonian class for band calculations

## References

1. W. A. Harrison, "Electronic Structure and the Properties of Solids", Dover (1989).
2. tbmodels documentation: https://tbmodels.greschd.ch/
