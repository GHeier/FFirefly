# Vertex Data Ordering Fix - Final Resolution

## Problem

When running `fly.x < 4_vertex.cfg`, the saved vertex file produced incorrect values:
```
Chi(q) = 0.215182
Vertex(q) = 0.133309          ← From saved file (COMPLETELY WRONG)
Expected Vertex(q) = 0.322731  ← Correct value from calculation
```

The saved/reloaded value was off by **142%**.

## Root Cause

**Data ordering mismatch** between how vertex data was computed vs. how it was stored:

### The System's Expectations

The `CMF_search_3d` function in `field_funcs.cpp` uses this indexing formula:
```cpp
f[w_idx * nx * ny + kx_idx * ny + ky_idx]
```

This expects **w-k ordering** where:
- Frequency (w) varies slowest (outermost in memory layout)
- k-points vary fastest (innermost)
- Layout: `[w0,k0], [w0,k1], ..., [w0,kN], [w1,k0], [w1,k1], ..., [w1,kN], ...`

### The Bug

The original vertex.cpp loop structure generated data in **k-w order**:
```cpp
for (int i = 0; i < nx; i++) {        // kx loop OUTER
    for (int j = 0; j < ny; j++) {    // ky loop
        for (int k = 0; k < nz; k++) {  // kz loop
            for (int l = 0; l < wpts.size(); l++) {  // w loop INNER
                vals.push_back(vertex_val);
            }
        }
    }
}
```

This produced: `[k0,w0], [k0,w1], ..., [k0,wN], [k1,w0], [k1,w1], ..., [k1,wN], ...`

### Why It Failed

1. Vertex data was computed and stored in **k-w order**
2. Data was saved with this wrong ordering
3. When reloaded, DataEvaluator assumed **w-k order**
4. Query for point (qx, qy, w) calculated index assuming w-k layout
5. Index pointed to completely wrong position in the k-w ordered array
6. Returned incorrect value (0.133309 instead of 0.322731)

## Solution

**Swap the loop order** to generate data in w-k order (lines 108-133 in vertex.cpp):

```cpp
// w loop outermost to generate data in w-k order (matching CMF_search expectations)
for (int l = 0; l < wpts.size(); l++) {  // Frequency loop OUTER
    float w = wpts[l];
    for (int i = 0; i < nx; i++) {        // kx loop
        for (int j = 0; j < ny; j++) {    // ky loop
            for (int k = 0; k < nz; k++) {  // kz loop INNER
                Vec q = brillouin_zone * Vec(i / nx - 0.5, j / ny - 0.5, k / nz - 0.5);
                q.dimension = chidim;

                cfloat X = chi(q, w);
                cfloat val = (U * U * X) / cfloat(1.0f - U * X) +
                             (U * U * U * X * X) / cfloat(1.0f - U * U * X * X);
                vals.push_back(val);

                if (real(U * X) >= 1.0) {
                    printf("Geometric series not convergent: U*X = %f\n", U * X.real());
                    exit(1);
                }
            }
        }
    }
}
```

Now generates: `[w0,k0], [w0,k1], ..., [w0,kN], [w1,k0], [w1,k1], ..., [w1,kN], ...`

## Verification

**Before fix:**
```
Vertex(q) = 0.133309
Expected Vertex(q) = 0.322731
Difference: 142% (WRONG!)
```

**After fix:**
```
Vertex(q) = 0.320494
Expected Vertex(q) = 0.322731
Difference: 0.69% (correct - interpolation error only)
```

✅ Values now match within interpolation tolerance!

## Important Note

The `as_mesh=false` branch (lines 88-107) already had the correct w-loop-outer structure and did not need fixing. Only the `as_mesh=true` branch (mesh data, lines 108-133) had the bug.

## Files Regenerated

All vertex files have been regenerated with correct ordering:
- ✅ `data/responses/n=0.85_vertex.h5`
- ✅ `data/responses/n=0.65_vertex.h5`

## Key Takeaway

**Always match data ordering to what the evaluation functions expect.** In FFirefly's CMField system:
- Use **w-k ordering** (frequency outermost) for mesh data
- This matches the `CMF_search_*d` indexing formula expectations
- Loop structure should be: `for w { for kx { for ky { for kz { ... } } } }`
