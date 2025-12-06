# BaseData and Field Examples

This directory contains example scripts demonstrating how to use Firefly's BaseData and Field classes in Python and Julia.

## Overview

Firefly's Field system provides interpolatable multi-dimensional tensor fields defined on k-space/frequency grids. The key concepts are:

- **BaseData**: Container class for storing tensor field data with metadata (mesh, domain, frequency points, tensor indices)
- **Field**: Wrapper classes (Field_C, Field_CM, Field_R, Field_RM) that provide interpolation and evaluation
- **inds**: Vector specifying tensor dimensions (NEW in this refactoring!)

## Key Changes in This Refactoring

The `inds` parameter replaces the old `(n_indices, dim_indices)` system:

### Old System
```cpp
n_indices = 2      // Number of indices
dim_indices = 3    // All indices have same dimension
// Could only represent 3×3 matrices
```

### New System
```cpp
inds = {3, 3}      // Each dimension specified individually
inds = {2, 3}      // Non-uniform dimensions now supported!
inds = {1, 2, 2, 1} // Complex vertex tensors
```

## Tensor Ranks and inds

| Rank | inds Example | Description | Field_CM Returns |
|------|--------------|-------------|------------------|
| 0 | `[]` | Scalar | N/A (use Field_C) |
| 1 | `[4]` | 4-component vector | `[4][1]` column matrix |
| 2 | `[3, 3]` | 3×3 matrix | `[3][3]` matrix |
| 2 | `[2, 3]` | 2×3 matrix (non-uniform!) | `[2][3]` matrix |
| 3 | `[2, 2, 2]` | 2×2×2 tensor | `[4][2]` flattened |
| 4 | `[2, 2, 2, 2]` | 2-orbital vertex | `[4][4]` flattened |
| 4 | `[1, 2, 2, 1]` | Non-uniform vertex | `[2][2]` flattened |

## Example Scripts

### Python

1. **example_basedata_fields_simple.py** (Recommended for beginners)
   - Shows basic Field operations using test data
   - Demonstrates load, evaluate, save workflow
   - Good starting point to understand the API

2. **example_basedata_fields.py** (Advanced)
   - Full examples of creating data from scratch
   - Covers all tensor ranks and types
   - Shows numpy array usage

### Julia

1. **example_basedata_fields.jl**
   - Comprehensive Julia examples
   - Covers all tensor types
   - Shows Julia-specific API usage

## Running the Examples

### Python

```bash
# Simple example (recommended)
python3 scripts/example_basedata_fields_simple.py

# Full examples (requires understanding data creation)
python3 scripts/example_basedata_fields.py
```

### Julia

```bash
# Make sure Julia package is precompiled first
julia --project=jlpkg/Firefly -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Run examples
julia scripts/example_basedata_fields.jl
```

## Creating Data from Scratch

For detailed examples of creating BaseData from scratch in different languages, see:

**C++ Examples** (most comprehensive):
```
src/objects/CMField/tests/inds_field_tests.cpp    # All tensor ranks
src/objects/CMField/tests/matrix_field_tests.cpp  # Matrix fields
src/objects/CMField/tests/tensor_field_tests.cpp  # 3D/4D tensors
```

**Python**:
- Use `save_data_scalar()`, `save_data_vector()`, `save_data_matrix()`, `save_data_tensor4()`
- See `src/module/imports/cpp_imports.py` for function signatures

**Julia**:
- Use `save_data_scalar()`, `save_data_vector()`, `save_data_matrix()`, `save_data_tensor4()`
- See `src/module/imports/cpp_imports.jl` for function signatures

## Data Ordering Convention

**Important**: Firefly uses **w-k ordering** (frequency varies slowest):

```python
# Correct w-k ordering for 2D k-space with frequency
for iw in range(nw):           # Frequency (slowest)
    for iky in range(nk_y):    # k_y
        for ikx in range(nk_x): # k_x (fastest)
            data[...] = ...
```

## Field Types

| Class | Type | Use For |
|-------|------|---------|
| `Field_R` | Real scalar | Real-valued functions |
| `Field_C` | Complex scalar | Complex-valued functions |
| `Field_RM` | Real matrix/tensor | Real matrices/tensors |
| `Field_CM` | Complex matrix/tensor | Complex matrices/tensors, vertices |

## Common Patterns

### Load and Evaluate

```python
# Python
import firefly.src.module.imports.cpp_imports as ff

# Load field
field = ff.Field_CM("my_data.h5")

# Evaluate at k-point
k = ff.Vec(0.1, 0.2, 0.3)
matrix = field(k)  # Returns 2D list

# With frequency
w = 0.5
matrix = field(k, w)
```

```julia
# Julia
using Firefly

# Load field
field = Field_CM("my_data.h5")

# Evaluate
k = Vec(0.1, 0.2, 0.3)
matrix = field(k)  # Returns Matrix

# With frequency
w = 0.5
matrix = field(k, w)
```

### Save and Load Round-trip

```python
# Python
field1 = ff.Field_CM("input.h5")
field1.save("output.h5")
field2 = ff.Field_CM("output.h5")

# Verify
k = ff.Vec(0.0)
assert abs(field1(k)[0][0] - field2(k)[0][0]) < 1e-6
```

### Access BaseData Properties

```python
# Python
basedata = ff.BaseData("data.h5")
print(f"Tensor rank: {len(basedata.inds)}")
print(f"Tensor dims: {basedata.inds}")
print(f"K-space dim: {basedata.dimension}")
print(f"K-points: {basedata.nk}")
print(f"Frequencies: {basedata.nw}")
```

## HDF5 File Format

BaseData saves to HDF5 with the following structure:

```
/data          # Flattened tensor data
/mesh          # K-space grid dimensions
/domain        # Brillouin zone vectors
/w_points      # Frequency points (optional)
/inds          # Tensor dimensions (NEW!)
/is_complex    # Boolean flag
/is_vector     # Boolean flag
/is_matrix     # Boolean flag
/with_k        # Boolean flag
/with_w        # Boolean flag
/dimension     # K-space dimensionality
```

You can inspect these files with:
```bash
h5ls -r data.h5
h5dump data.h5
```

## Tips and Best Practices

1. **Always specify inds explicitly** for matrix/tensor fields
2. **Use w-k ordering** when creating data
3. **Check tensor ranks** with `len(basedata.inds)` before evaluation
4. **Use appropriate Field type**: Field_C for scalars, Field_CM for matrices
5. **Test round-trips** to verify data integrity

## Troubleshooting

**"BaseData: type mismatch" error**:
- Check that `inds` is set correctly for your data
- Verify w-k ordering in your data
- Ensure data shape matches `mesh * total_index_size`

**Wrong result shape**:
- Field_CM flattens higher-rank tensors to matrices
- Check expected output shape based on rank (see table above)

**Import errors**:
- Make sure you've built the project: `./scripts/fly-build.sh`
- Check library path for Python/Julia bindings

## Further Reading

- **Architecture**: See `CLAUDE.md` in project root
- **C++ API**: See `src/objects/CMField/fields.hpp`
- **Test Suite**: See `src/objects/CMField/tests/inds_field_tests.cpp`
- **Python API**: See `src/module/imports/cpp_imports.py`
- **Julia API**: See `src/module/imports/cpp_imports.jl`
