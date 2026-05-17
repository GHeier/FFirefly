# hamiltonian / DOS / tetrahedra

## Overview

Calculates the electronic density of states (DOS) using the tetrahedron integration method. This approach divides the Brillouin zone into tetrahedra and analytically integrates the band structure to obtain the DOS at each energy point. Also computes the integrated electron number as a function of chemical potential.

## Quick Description

Computes the Density of States using surface construction at discrete w-points.

## Dependencies

### Required
- Band structure model (tight_binding or fermi_gas)
- LAPACK/OpenBLAS for linear algebra
- HDF5 for output

### Optional
- None

## Install Instructions

```bash
# No special installation required - built automatically by fly-build.sh
```

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}/{prefix}_DOS.h5` - Density of states D(E) vs energy
- `{outdir}/{prefix}_E_vs_n.h5` - Integrated electron number n(E) vs chemical potential

File format details:
- Both files are 1D real fields stored in HDF5 format
- Energy grid spans from band minimum to maximum with `w_pts` points

## Testing

Run the test suite:
```bash
fly.x  # Runs all tests including this one
```

Expected test behavior:
- For 3D tight-binding at half-filling, DOS should show van Hove singularities
- Integrated electron number should equal 1.0 at the band center for half-filled single band

## Calculation Details

### Algorithm

1. Determine energy range from band structure minimum/maximum across k-mesh
2. Create energy grid with `w_pts` points spanning the band range
3. For each energy point E, compute the constant-energy surface (isoenergy contour)
4. Integrate over the surface using: DOS(E) = (1/(2π)^d) * ∫ dS / |∇ε(k)|
5. Compute cumulative electron number by integrating DOS over energy

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | - | k-point mesh dimensions [nx, ny, nz] |
| `dimension` | int | 3 | Spatial dimensionality (1, 2, or 3) |
| `w_pts` | int | 500 | Number of energy points for DOS |
| `fermi_energy` | float | 0.0 | Fermi level in eV |

### Implementation Notes

- Uses tetrahedron interpolation for smooth DOS without artificial broadening
- Van Hove singularities are captured accurately
- More efficient than Gaussian smearing for large k-meshes

## References

1. P. E. Blochl, O. Jepsen, and O. K. Andersen, "Improved tetrahedron method for Brillouin-zone integrations", Phys. Rev. B 49, 16223 (1994).
2. G. Lehmann and M. Taut, "On the Numerical Calculation of the Density of States and Related Properties", Phys. Status Solidi B 54, 469 (1972).
