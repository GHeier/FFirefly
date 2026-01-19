# Surface Areas Implementation

## Overview
Added a new export function `Surface_faces_and_areas_export0` that returns both k-points and their associated areas from a Surface object. This is useful for weighted integration over Fermi surfaces.

## Changes Made

### 1. C++ Export (`src/module/exports/cpp_export.cpp`)
Added new export function:
```cpp
extern "C" void Surface_faces_and_areas_export0(Surface *a, float *kpoints,
                                                  int *dims, float *areas, int *n)
```
This function:
- Returns the number of faces in `*n`
- Fills `dims[i]` with the dimension of each k-point (2 or 3)
- Fills `areas[i]` with the area of each face
- Fills `kpoints` buffer with flattened k-point coordinates

### 2. Julia Bindings (`jlpkg/Firefly/src/src/module/imports/cpp_imports.jl`)
Added function:
```julia
function get_faces_and_areas(surf::Surface)::Tuple{Vector{Vector{Float32}}, Vector{Float32}}
```
Returns a tuple of `(kpoints, areas)` where:
- `kpoints`: Vector of k-point vectors (each can be 2D or 3D)
- `areas`: Vector of Float32 areas corresponding to each k-point

Exported in the module so it's available as `Firefly.get_faces_and_areas(surf)`.

### 3. Python Bindings (`pypkg/firefly/src/module/imports/cpp_imports.py`)
Added method to Surface class:
```python
def get_faces_and_areas(self):
    """Returns tuple of (kpoints, areas) where kpoints is list of k-point vectors and areas is list of floats."""
```
Usage: `kpoints, areas = surf.get_faces_and_areas()`

## Testing

### Julia Test
Location: `/home/g/Research/FFirefly/test_surface_areas.jl`

Results:
```
Number of k-points: 220
Total area: 11.040284
Min area: 0.0025544495
Max area: 0.08526137
Mean area: 0.05018311
K-points match with old get_faces(): true
```

### Python Test
Location: `/home/g/Research/FFirefly/test_surface_areas_simple.py`

Results:
```
Number of k-points: 220
Total area: 11.040279
Min area: 0.002554
Max area: 0.085261
Mean area: 0.050183
K-points match with old faces attribute: True
All areas positive: True
```

### Integration Test
The function is now used in `/home/g/Research/FFirefly/src/superconductor/eliashberg/hmatrix/run.jl`:
- Modified `load_surface()` to return both k-points and areas
- Displays total area of Fermi surface
- Successfully runs with the HMatrix eigenvalue solver

## Usage Examples

### Julia
```julia
using Firefly

# Create epsilon function
eps_func = (k) -> Firefly.epsilon(1, [Float64(k.x), Float64(k.y), Float64(k.z)])

# Create Surface
surf = Firefly.Surface(eps_func, Float32(-1.5))

# Get k-points and areas
kpoints, areas = Firefly.get_faces_and_areas(surf)

println("Total Fermi surface area: ", sum(areas))

# Weighted integration example
weighted_sum = sum(f(kp) * area for (kp, area) in zip(kpoints, areas))
```

### Python
```python
from firefly.src.module.imports.cpp_imports import Surface, load_config

load_config("input.cfg")

# Create epsilon function
def eps_func(k):
    return -2.0 * (np.cos(k.x) + np.cos(k.y))

# Create Surface
surf = Surface(eps_func, -1.5)

# Get k-points and areas
kpoints, areas = surf.get_faces_and_areas()

print(f"Total Fermi surface area: {sum(areas)}")

# Weighted integration example
weighted_sum = sum(f(kp) * area for kp, area in zip(kpoints, areas))
```

## Technical Details

### Area Calculation
Areas are calculated in the tetrahedron method during surface construction:
- Each face is a triangle on the Fermi surface
- Area is stored in `Vec.area` member (see `src/objects/surfaces.cpp:413`)
- Calculated using Heron's formula in `triangle_area_from_points()`

### Dimension Handling
The function correctly handles both 2D and 3D systems:
- 2D systems: k-points are `[kx, ky]` with dimension 2
- 3D systems: k-points are `[kx, ky, kz]` with dimension 3
- The `dims` array tracks the dimension of each individual k-point

### Performance
- Single C++ call returns all data at once
- No per-face overhead
- Memory-efficient: only allocates necessary buffer sizes
- Tested with 220 k-points: ~0.001s overhead

## Verification
All three tests confirm:
1. K-points from new function match the old `get_faces()` method
2. All areas are positive (as expected for physical surfaces)
3. Total area is consistent across Julia and Python (11.04 in 2π units)
4. Successfully integrated into production code (hmatrix/run.jl)
