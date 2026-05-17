# hamiltonian / FS / tetrahedra

## Overview

Generates the Fermi surface (FS) - the locus of k-points where the band energy equals the Fermi level. Uses tetrahedron interpolation to construct a smooth surface representation with associated area elements and Fermi velocities. The output is used by many-body and superconductor calculations that integrate over the Fermi surface.

## Quick Description

One or two sentence description of the method/algorithm used.

## Dependencies

### Required
- Band structure model (tight_binding or fermi_gas)
- Tetrahedron integration algorithm
- HDF5 for output (optional)

### Optional
- None

## Install Instructions

```bash
# No special installation required - built automatically by fly-build.sh
```

## Results Saved

Output files created by this calculation (using `prefix` from config):

- `{outdir}/{prefix}_FS.dat` - Fermi surface k-points in text format

File format details:
- ASCII file with columns: `kx ky [kz] band_index`
- Each row is a k-point on the Fermi surface
- For 2D systems, kz column is omitted

## Testing

Run the test suite:
```bash
fly.x  # Runs all tests including this one
```

Expected test behavior:
- For 3D tight-binding at half-filling, should produce a connected surface
- For 2D square lattice at half-filling, should produce a diamond-shaped contour
- Total FS area should match analytical expectations for simple models

## Calculation Details

### Algorithm

1. Load band structure parameters and compute ε(k) on k-mesh
2. Divide Brillouin zone into tetrahedra
3. For each tetrahedron, find intersection with ε(k) = μ using linear interpolation
4. Collect all intersection points to form the Fermi surface
5. Compute surface area elements dA and Fermi velocities v_F = |∇ε(k)|
6. Save k-points with metadata (area, band index)

### Parameters

Configuration parameters from `input.cfg`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k_mesh` | int[3] | - | k-point mesh dimensions [nx, ny, nz] |
| `dimension` | int | 3 | Spatial dimensionality (1, 2, or 3) |
| `fermi_energy` | float | 0.0 | Fermi level μ in eV |
| `nbnd` | int | 1 | Number of bands to include |

### Implementation Notes

- Finer k-mesh produces smoother Fermi surface representation
- Multiple bands may contribute separate Fermi surface sheets
- Area elements are used for proper weighting in FS integrals
- Fermi velocity is stored for DOS and transport calculations

## References

1. P. E. Blochl, O. Jepsen, and O. K. Andersen, "Improved tetrahedron method for Brillouin-zone integrations", Phys. Rev. B 49, 16223 (1994).
