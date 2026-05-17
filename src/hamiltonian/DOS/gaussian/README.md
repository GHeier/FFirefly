# hamiltonian / DOS / gaussian

## Overview

Computes the electronic density of states (DOS) using Gaussian smearing. Each k-point contributes a Gaussian-broadened delta function at its band energy, yielding a smooth DOS that avoids singularities. Also computes the integrated electron number as a function of chemical potential.

## Quick Description

Computes the Density of States with gaussian spreading for smooth results.

## Dependencies

- NumPy for vectorized operations
- Firefly Bands class for band energies
- HDF5 for output

## Install Instructions

```bash
# No special installation required
pip install numpy h5py
```

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | - | k-point mesh dimensions [nx, ny, nz] |
| `dimension` | int | 3 | Spatial dimensionality (1, 2, or 3) |
| `w_pts` | int | 500 | Number of energy points for DOS |
| `smearing` | float | 0.05 | Gaussian broadening width σ in eV |

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}/{prefix}_DOS.h5` - Density of states D(E) vs energy
- `{outdir}/{prefix}_E_vs_n.h5` - Integrated electron number n(E) vs chemical potential

File format details:
- Both files are 1D real fields stored in HDF5 format
- Energy grid spans from band minimum to maximum with `w_pts` points

## Testing

Expected test behavior:
- DOS should be smooth and positive everywhere
- For 3D tight-binding, should qualitatively match tetrahedra result but with broadened singularities
- Integrated electron number should approach total band filling at high energies

## Calculation Details

### Algorithm

1. Create uniform k-point mesh in the first Brillouin zone
2. Evaluate band energies ε(k) at all k-points
3. Create energy grid from min(ε) to max(ε) with `w_pts` points
4. For each energy E, compute: DOS(E) = (1/N_k) Σ_k exp(-(E-ε_k)²/(2σ²)) / (σ√(2π))
5. Normalize DOS to integrate to 1 over the energy range
6. Integrate DOS to get electron number vs chemical potential

### Implementation Notes

- Gaussian smearing produces inherently smooth DOS
- Larger σ gives smoother results but washes out fine structure
- More computationally intensive than tetrahedra for the same accuracy
- Useful when tetrahedra method has numerical issues

## References

1. M. Methfessel and A. T. Paxton, "High-precision sampling for Brillouin-zone integration in metals", Phys. Rev. B 40, 3616 (1989).
