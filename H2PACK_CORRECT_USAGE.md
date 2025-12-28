# H2Pack Correct Usage Guide

## Key Discovery: Column-Major Coordinate Format

**CRITICAL**: H2Pack expects coordinates in **column-major (dimension-contiguous)** format, NOT row-major (interleaved) format.

### Coordinate Layout

```
CORRECT (Column-Major):
coord = [x0, x1, x2, ..., xN,  y0, y1, y2, ..., yN,  z0, z1, z2, ..., zN]
        |<---  all X  --->|  |<---  all Y  --->|  |<---  all Z  --->|

WRONG (Row-Major/Interleaved):
coord = [x0, y0, z0,  x1, y1, z1,  x2, y2, z2, ...]
```

### Code Example

**Generating Points (Column-Major):**
```cpp
int n_point = 500;
int pt_dim = 3;
vector<double> coord(n_point * pt_dim);

for (int i = 0; i < n_point; i++) {
    double x = ..., y = ..., z = ...;

    // Column-major storage
    coord[0 * n_point + i] = x;  // All x coordinates
    coord[1 * n_point + i] = y;  // All y coordinates
    coord[2 * n_point + i] = z;  // All z coordinates
}
```

### Kernel Function (Column-Major Access)

Based on H2Pack's `EXTRACT_3D_COORD()` macro from `H2Pack_3D_kernels.h`:

```cpp
extern "C" void my_kernel(
    const double* coord0, const int ld0, const int n0,
    const double* coord1, const int ld1, const int n1,
    const void* param, double* out_mat, const int ldm
) {
    // Extract coordinates (ld0 = n_points, leading dimension)
    const double *x0 = coord0 + ld0 * 0;  // Offset 0
    const double *y0 = coord0 + ld0 * 1;  // Offset ld0
    const double *z0 = coord0 + ld0 * 2;  // Offset 2*ld0
    const double *x1 = coord1 + ld1 * 0;
    const double *y1 = coord1 + ld1 * 1;
    const double *z1 = coord1 + ld1 * 2;

    for (int i = 0; i < n0; i++) {
        for (int j = 0; j < n1; j++) {
            double dx = x0[i] - x1[j];
            double dy = y0[i] - y1[j];
            double dz = z0[i] - z1[j];
            double r2 = dx*dx + dy*dy + dz*dz;
            out_mat[i * ldm + j] = 1.0 / (1.0 + r2);
        }
    }
}
```

## Standard H2Pack Workflow

```cpp
// 1. Initialize
H2Pack_p h2pack;
int krnl_dim = 1;
H2P_init(&h2pack, pt_dim, krnl_dim, QR_REL_NRM, &rel_tol);

// 2. Calculate enclosing box
H2P_calc_enclosing_box(pt_dim, n_point, coord.data(), nullptr, &h2pack->root_enbox);

// 3. Partition points into tree
H2P_partition_points(h2pack, n_point, coord.data(), 0, 0.0);

// 4. Generate proxy points
H2P_dense_mat_p *pp;
H2P_generate_proxy_point_ID_file(h2pack, nullptr,
    (kernel_eval_fptr)my_kernel, nullptr, &pp);

// 5. Build H2 representation
H2P_build(h2pack, pp, 0, nullptr,
    (kernel_eval_fptr)my_kernel, nullptr, 0);

// 6. Use the H2 matrix (matvec, etc.)
vector<double> x(n_point), y(n_point);
H2P_matvec(h2pack, x.data(), y.data());

// 7. Clean up
H2P_destroy(&h2pack);
```

## Test Results

### test_h2matrix_simple.cpp
- ✅ Compiles without errors
- ✅ Runs without crashes
- ✅ Builds tree properly (max_level=1, n_node=9)
- ✅ Matvec produces valid numeric results
- ⚠️ LAPACK DGEMV warnings appear but results are still correct

### Compilation
```bash
g++ -std=c++17 -fopenmp \\
    -I/home/g/Research/H2Pack/include \\
    test_h2matrix_simple.cpp \\
    -L/home/g/Research/H2Pack/lib \\
    -lH2Pack -lopenblas -llapacke -lm \\
    -o test_h2matrix_simple
```

### Execution
```bash
LD_LIBRARY_PATH=/home/g/Research/H2Pack/lib:$LD_LIBRARY_PATH \\
    ./test_h2matrix_simple
```

### Output
```
H2Pack Direct Test (Column-Major Format)
=========================================
Creating 500 points on a sphere...
Initializing H2Pack...
  Enclosing box calculated
  Partitioned: max_level=1, n_node=9
  Generating proxy points...
  Building H2 matrix...

SUCCESS: H2 matrix built!
  Tree depth: 1 levels
  Tree nodes: 9 nodes

✓ Tree structure verified!

Testing matrix-vector multiplication...
  Computing H2 matvec...
  H2 matvec completed (sum |y|: 141.127)
  Computing dense matvec for comparison...

Comparing H2 vs Dense results...
  Max absolute error: 2.07889e-14
  RMS error: 5.75802e-15
  Max relative error: 3.21924e-12

  Sample comparison (first 5 points):
  i    H2 result        Dense result      Error
  0    0.6313441929     0.6313441929     3.22e-15
  1    0.6814681960     0.6814681960     1.11e-16
  2    0.7117857500     0.7117857500     4.22e-15
  3    0.7188466160     0.7188466160     2.00e-15
  4    0.7021977782     0.7021977782     4.44e-15

✓ H2 matvec matches dense within tolerance!

========================================
ALL TESTS PASSED!
========================================
```

### Numerical Accuracy

The test confirms H2Pack produces **numerically accurate results**:

- **Max error**: ~2e-14 (machine precision)
- **Agreement**: H2 matvec ≈ Dense matvec to ~14-15 digits
- **Validation**: Dense matrix-vector product computed directly for comparison
- **Conclusion**: Column-major format and kernel implementation are correct

## Common Mistakes to Avoid

### ❌ WRONG: Row-major coordinates
```cpp
coord[i * pt_dim + 0] = x;  // Interleaved
coord[i * pt_dim + 1] = y;
coord[i * pt_dim + 2] = z;
```

### ❌ WRONG: Kernel accessing row-major
```cpp
int dim = ld0;  // ld0 is n_points, not dimension!
double diff = coord0[i * ld0 + d] - coord1[j * ld1 + d];
```

### ✅ CORRECT: Column-major coordinates
```cpp
coord[0 * n_point + i] = x;  // All x's together
coord[1 * n_point + i] = y;  // All y's together
coord[2 * n_point + i] = z;  // All z's together
```

### ✅ CORRECT: Kernel accessing column-major
```cpp
const double *x0 = coord0 + ld0 * 0;
const double *y0 = coord0 + ld0 * 1;
const double *z0 = coord0 + ld0 * 2;
double dx = x0[i] - x1[j];
```

## References

- H2Pack repository: `/home/g/Research/H2Pack`
- Official examples: `/home/g/Research/H2Pack/examples/example_H2.c`
- Kernel macros: `/home/g/Research/H2Pack/include/H2Pack_3D_kernels.h`
- Working test: `/home/g/Research/FFirefly/test_h2matrix_simple.cpp`
