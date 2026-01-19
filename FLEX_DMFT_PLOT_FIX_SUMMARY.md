# FLEX_DMFT_plot.py SIGSEGV Fix Summary

## Problem Identified

The SIGSEGV error in `FLEX_DMFT_plot.py` was caused by using the wrong field type (`Field_CM`) for scalar fields.

### Root Cause

In `FLEX_DMFT_plot.py` lines 35 and 43, the code was:
```python
field = fly.Field_CM("data/vertex/FLEX_DMFT_sigma_loc.h5")
...
field = fly.Field_CM("data/vertex/FLEX_DMFT_sigma_imp.h5")
```

However, these files are **scalar fields** (rank=0, `inds=[]`), not matrix fields:
- `FLEX_DMFT_sigma_loc.h5`: frequency-only scalar field (with_k=0, with_w=1, inds=[])
- `FLEX_DMFT_sigma_imp.h5`: frequency-only scalar field (with_k=0, with_w=1, inds=[])

Using `Field_CM` (Complex Matrix Field) to load scalar data causes a segmentation fault when evaluating the field, particularly after other `Field_CM` objects have been loaded.

## Solution

Changed the field type from `Field_CM` to `Field_C` for scalar fields:

**Before:**
```python
field = fly.Field_CM("data/vertex/FLEX_DMFT_sigma_loc.h5")
...
y = np.array(field(w_pts))
...
field = fly.Field_CM("data/vertex/FLEX_DMFT_sigma_imp.h5")
y = np.array(field(w_pts))[:, 0, 0]  # This slicing was unnecessary
```

**After:**
```python
field = fly.Field_C("data/vertex/FLEX_DMFT_sigma_loc.h5")  # Changed to Field_C
...
y = np.array(field(w_pts))
...
field = fly.Field_C("data/vertex/FLEX_DMFT_sigma_imp.h5")  # Changed to Field_C
y = np.array(field(w_pts))  # Removed unnecessary slicing
```

## Field Type Selection Guide

To determine which field type to use:

1. Check the `inds` dataset in the HDF5 file:
   ```bash
   h5dump -d /inds filename.h5
   ```

2. Choose field type based on rank:
   - `inds = []` (rank=0): Use `Field_C` or `Field_R` (scalar)
   - `inds = [n]` (rank=1): Use vector field (if available)
   - `inds = [n, n]` (rank=2): Use `Field_CM` or `Field_RM` (matrix)
   - `inds = [n, n, n]` (rank=3+): Use tensor field (if available)

## Verification

After the fix:
- ✅ `FLEX_DMFT_plot.py` runs without segfault
- ✅ Generates `plots/vertex_shape.png` successfully
- ✅ All four subplots render correctly

## Key Takeaway

**Always match the field loader type to the data rank:**
- Scalar data → `Field_C` or `Field_R`
- Matrix data → `Field_CM` or `Field_RM`

Using the wrong field type will cause segmentation faults, especially when multiple fields are loaded in sequence.
