# SIMD and Prefetching Explained

## Part 1: SIMD (Single Instruction, Multiple Data)

### What is SIMD?

SIMD allows you to perform the **same operation on multiple values simultaneously** using special CPU instructions.

**Normal (scalar) code:**
```cpp
float a[4] = {1.0, 2.0, 3.0, 4.0};
float b[4] = {5.0, 6.0, 7.0, 8.0};
float c[4];

// Process one at a time (4 operations)
c[0] = a[0] + b[0];  // Takes 1 clock cycle
c[1] = a[1] + b[1];  // Takes 1 clock cycle
c[2] = a[2] + b[2];  // Takes 1 clock cycle
c[3] = a[3] + b[3];  // Takes 1 clock cycle
// Total: 4 clock cycles
```

**SIMD code:**
```cpp
#include <immintrin.h>  // AVX/SSE intrinsics

float a[4] = {1.0, 2.0, 3.0, 4.0};
float b[4] = {5.0, 6.0, 7.0, 8.0};
float c[4];

// Process all four at once (1 operation)
__m128 va = _mm_load_ps(a);      // Load 4 floats into SIMD register
__m128 vb = _mm_load_ps(b);      // Load 4 floats into SIMD register
__m128 vc = _mm_add_ps(va, vb);  // Add all 4 pairs simultaneously!
_mm_store_ps(c, vc);             // Store result

// Total: 1 clock cycle for 4 additions
// Speedup: 4×
```

### How CPU Executes SIMD

**Normal scalar registers:**
```
CPU Register (32-bit):  [  42.0  ]
                         1 float
```

**SIMD registers (SSE/AVX):**
```
SSE Register (128-bit):  [  1.0  |  2.0  |  3.0  |  4.0  ]
                          4 floats packed together

AVX Register (256-bit):  [  1.0  |  2.0  |  3.0  |  4.0  |  5.0  |  6.0  |  7.0  |  8.0  ]
                          8 floats packed together

AVX-512 (512-bit):       [  16 floats packed together  ]
```

When you execute `_mm_add_ps`, the CPU has **dedicated hardware** that adds all 4 pairs in parallel.

### SIMD Instruction Sets

| ISA | Year | Register Size | Floats | Doubles |
|-----|------|--------------|--------|---------|
| SSE | 1999 | 128-bit | 4 | 2 |
| SSE2 | 2001 | 128-bit | 4 | 2 |
| AVX | 2011 | 256-bit | 8 | 4 |
| AVX2 | 2013 | 256-bit | 8 | 4 |
| AVX-512 | 2016 | 512-bit | 16 | 8 |

Modern CPUs have all of these. Your code should target **AVX2** (widely available).

### Example: Interpolation with SIMD

Let's rewrite the trilinear interpolation using SIMD:

**Scalar version (current):**
```cpp
complex<float> result =
    (1-x) * (1-y) * (1-w) * f[0] +  // 3 multiplies, 1 load
    x     * (1-y) * (1-w) * f[1] +  // 3 multiplies, 1 load
    (1-x) * y     * (1-w) * f[2] +  // 3 multiplies, 1 load
    x     * y     * (1-w) * f[3] +  // 3 multiplies, 1 load
    (1-x) * (1-y) * w     * f[4] +  // 3 multiplies, 1 load
    x     * (1-y) * w     * f[5] +  // 3 multiplies, 1 load
    (1-x) * y     * w     * f[6] +  // 3 multiplies, 1 load
    x     * y     * w     * f[7];   // 3 multiplies, 1 load

// Total: 8 loads, 24 multiplies, 7 additions
```

**SIMD version (AVX - real part only for clarity):**
```cpp
#include <immintrin.h>

// Step 1: Load all 8 values at once (in two 256-bit loads)
__m256 values_lo = _mm256_set_ps(
    f[3].real(), f[2].real(), f[1].real(), f[0].real(),
    f[7].real(), f[6].real(), f[5].real(), f[4].real()
);

// Step 2: Compute all 8 weights
__m256 weights = _mm256_set_ps(
    x * y * (1-w),        // f[3]
    (1-x) * y * (1-w),    // f[2]
    x * (1-y) * (1-w),    // f[1]
    (1-x) * (1-y) * (1-w),// f[0]
    x * y * w,            // f[7]
    (1-x) * y * w,        // f[6]
    x * (1-y) * w,        // f[5]
    (1-x) * (1-y) * w     // f[4]
);

// Step 3: Multiply all 8 pairs at once
__m256 products = _mm256_mul_ps(values_lo, weights);

// Step 4: Horizontal sum (reduce to scalar)
// Sum pairs: [a,b,c,d,e,f,g,h] -> [a+b, c+d, e+f, g+h, ...]
__m256 sum1 = _mm256_hadd_ps(products, products);
// Sum pairs again: -> [a+b+c+d, e+f+g+h, ...]
__m256 sum2 = _mm256_hadd_ps(sum1, sum1);
// Extract and add high/low 128-bit lanes
__m128 lo = _mm256_extractf128_ps(sum2, 0);
__m128 hi = _mm256_extractf128_ps(sum2, 1);
__m128 final = _mm_add_ps(lo, hi);

float result_real = _mm_cvtss_f32(final);

// Repeat for imaginary part...
```

**Speedup:**
- Scalar: 24 multiplies sequentially = 24 cycles
- SIMD: 8 multiplies in 1 instruction = 1 cycle
- **Theoretical: 24× faster**
- **Practical: 2-4× faster** (memory bandwidth limits, overhead)

### Practical SIMD Example for Our Use Case

After implementing W-first reorganization, the SIMD code becomes much cleaner:

```cpp
// Helper function: 2D interpolation with SIMD
inline complex<float> interpolate_2d_simd(
    const complex<float>* slice,
    int i, int j, int ny,
    float x_rel, float y_rel
) {
    // Load 4 values (2x2 spatial stencil)
    // Real parts
    __m128 real_vals = _mm_set_ps(
        slice[(i+1)*ny + (j+1)].real(),  // [3]
        slice[i*ny + (j+1)].real(),      // [2]
        slice[(i+1)*ny + j].real(),      // [1]
        slice[i*ny + j].real()           // [0]
    );

    // Compute weights
    __m128 weights = _mm_set_ps(
        x_rel * y_rel,           // w11
        (1-x_rel) * y_rel,       // w01
        x_rel * (1-y_rel),       // w10
        (1-x_rel) * (1-y_rel)    // w00
    );

    // Multiply and sum
    __m128 products = _mm_mul_ps(real_vals, weights);
    __m128 sum1 = _mm_hadd_ps(products, products);  // [a+b, c+d, a+b, c+d]
    __m128 sum2 = _mm_hadd_ps(sum1, sum1);          // [a+b+c+d, ...]
    float result_real = _mm_cvtss_f32(sum2);

    // Repeat for imaginary...
    __m128 imag_vals = _mm_set_ps(
        slice[(i+1)*ny + (j+1)].imag(),
        slice[i*ny + (j+1)].imag(),
        slice[(i+1)*ny + j].imag(),
        slice[i*ny + j].imag()
    );
    __m128 imag_products = _mm_mul_ps(imag_vals, weights);
    __m128 imag_sum1 = _mm_hadd_ps(imag_products, imag_products);
    __m128 imag_sum2 = _mm_hadd_ps(imag_sum1, imag_sum1);
    float result_imag = _mm_cvtss_f32(imag_sum2);

    return complex<float>(result_real, result_imag);
}

// Main interpolation function
complex<Vec> CMF_search_3d_simd(...) {
    // ... (compute indices i, j, k, weights) ...

    // Use SIMD for 2D interpolation at each frequency
    complex<float> val_k  = interpolate_2d_simd(&f[k*nx*ny], i, j, ny, x_rel, y_rel);
    complex<float> val_k1 = interpolate_2d_simd(&f[(k+1)*nx*ny], i, j, ny, x_rel, y_rel);

    // Linear interpolate in frequency
    return (1 - w_rel) * val_k + w_rel * val_k1;
}
```

### Auto-Vectorization (Easier!)

Modern compilers can automatically generate SIMD code if you write simple loops:

```cpp
// Compiler-friendly code
void interpolate_batch(const float* data, const float* weights,
                       float* results, int n) {
    // Simple loop - compiler will vectorize this!
    #pragma omp simd  // Hint to compiler: please vectorize
    for (int i = 0; i < n; i++) {
        results[i] = data[i] * weights[i];
    }
}

// Compiler generates AVX code automatically:
// - Processes 8 elements per iteration
// - Uses _mm256_mul_ps internally
```

**Benefits:**
- No manual intrinsics
- Portable across architectures
- Compiler chooses best instruction set

**To enable:**
```bash
g++ -O3 -march=native -ftree-vectorize your_code.cpp
```

Check if vectorized:
```bash
g++ -O3 -march=native -ftree-vectorize -fopt-info-vec your_code.cpp
# Output: "loop vectorized" if successful
```

---

## Part 2: Prefetching

### What is Prefetching?

Prefetching tells the CPU to **load data into cache before you actually need it**.

### Memory Access Latency

```
L1 Cache:    ~4 cycles   (~1 ns)    ← Very fast
L2 Cache:    ~12 cycles  (~3 ns)    ← Fast
L3 Cache:    ~40 cycles  (~10 ns)   ← Slower
RAM:         ~200 cycles (~50 ns)   ← SLOW!
```

**Problem:** If CPU has to wait for RAM, it sits idle for 200 cycles!

**Solution:** Tell CPU to start loading data early, so it arrives when needed.

### Hardware Prefetcher

Modern CPUs have automatic prefetchers that detect patterns:

```cpp
// Sequential access - hardware prefetcher detects pattern
for (int i = 0; i < n; i++) {
    sum += array[i];  // CPU prefetches array[i+4], array[i+8], etc.
}
// Result: Almost no cache misses!
```

**But hardware prefetcher fails on:**
- Random access patterns
- Complex strided access
- Pointer chasing
- Irregular patterns

### Software Prefetching

You explicitly tell CPU what to load:

```cpp
#include <xmmintrin.h>  // For _mm_prefetch

// Without prefetch
for (int i = 0; i < n; i++) {
    result[i] = compute(data[indices[i]]);  // Cache miss every time!
}

// With prefetch
for (int i = 0; i < n; i++) {
    // Prefetch data that we'll need in ~8 iterations
    if (i + 8 < n) {
        _mm_prefetch(&data[indices[i + 8]], _MM_HINT_T0);
    }

    // By the time we get here, data is (hopefully) in cache
    result[i] = compute(data[indices[i]]);
}
```

### Prefetch Hints

```cpp
_MM_HINT_T0  // Prefetch to L1 cache (use soon, <10 iterations)
_MM_HINT_T1  // Prefetch to L2 cache (use medium-term, 10-100 iterations)
_MM_HINT_T2  // Prefetch to L3 cache (use later, 100+ iterations)
_MM_HINT_NTA // Non-temporal (don't pollute cache, use once)
```

### Example: Prefetch for Interpolation

```cpp
vector<complex<float>> Field_C::operator()(const vector<Vec>& points, float w) {
    vector<complex<float>> results(points.size());

    // Prefetch distance: how many iterations ahead to prefetch
    const int PREFETCH_DISTANCE = 16;

    for (int i = 0; i < points.size(); i++) {
        // Prefetch future point
        if (i + PREFETCH_DISTANCE < points.size()) {
            Vec future_point = points[i + PREFETCH_DISTANCE];

            // Transform and compute base index
            Vec p = transform_point(future_point);
            int base_idx = compute_index(p, w);

            // Prefetch the memory we'll access
            _mm_prefetch(&data[base_idx], _MM_HINT_T0);
            _mm_prefetch(&data[base_idx + 200], _MM_HINT_T0);  // Next row
        }

        // Process current point (data hopefully in cache)
        results[i] = interpolate_single(points[i], w);
    }

    return results;
}
```

### Why PREFETCH_DISTANCE = 16?

```
Memory latency: ~200 cycles
Per interpolation: ~50 cycles

Cycles until we need data: 16 iterations × 50 cycles = 800 cycles
Prefetch latency: ~200 cycles

800 > 200 ✓  Data arrives in time!
```

Too small: Data not ready when needed
Too large: Data evicted from cache before use

### Advanced: Grouped Prefetching

For our case with two frequency slices:

```cpp
// Prefetch both frequency slices for future point
if (i + 16 < points.size()) {
    Vec future = points[i + 16];
    int k = find_frequency_index(future.w);

    // Prefetch from both frequency slices
    int base_k  = k * nx * ny;
    int base_k1 = (k + 1) * nx * ny;

    int idx = compute_spatial_index(future.x, future.y);

    // Prefetch 4 cache lines (covers all 8 neighbor points)
    _mm_prefetch(&f[base_k + idx], _MM_HINT_T0);
    _mm_prefetch(&f[base_k + idx + 8], _MM_HINT_T0);
    _mm_prefetch(&f[base_k1 + idx], _MM_HINT_T0);
    _mm_prefetch(&f[base_k1 + idx + 8], _MM_HINT_T0);
}
```

### Prefetch Performance Impact

**Without prefetch:**
```
Timeline:
  Cycle 0:   Start interpolation for point i
  Cycle 50:  Need data[k*nx*ny + i*ny + j]
  Cycle 50:  Cache miss! Wait for memory...
  Cycle 250: Data arrives (200 cycle penalty)
  Cycle 300: Complete interpolation
```

**With prefetch (distance = 16):**
```
Timeline (point i):
  Cycle 0:   Start interpolation for point i
             (Data was prefetched 16 iterations ago)
  Cycle 50:  Need data[...] - Already in L1 cache!
  Cycle 100: Complete interpolation (no wait)

Timeline (point i+16):
  Cycle 800: Issue prefetch for point i+16
             CPU starts loading in background
```

**Speedup:** Eliminates 200-cycle wait = 2-4× faster on cache-miss-heavy code

### Checking if Prefetch Helps

Use `perf` to measure cache misses:

```bash
# Without prefetch
perf stat -e cache-misses,cache-references ./profile_interp
#   50,000 cache misses (50% miss rate)

# With prefetch
perf stat -e cache-misses,cache-references ./profile_interp_prefetch
#   10,000 cache misses (10% miss rate)  ← 5× reduction!
```

---

## Part 3: Combining W-First + SIMD + Prefetch

### Complete Optimized Implementation

```cpp
#include <immintrin.h>  // AVX intrinsics

// Helper: SIMD 2D interpolation
inline __m128 interpolate_2d_slice_simd(
    const complex<float>* slice,
    int i, int j, int ny,
    float x_rel, float y_rel,
    bool get_real  // true for real, false for imag
) {
    // Extract component from all 4 neighbors
    float v00 = get_real ? slice[i*ny + j].real() : slice[i*ny + j].imag();
    float v10 = get_real ? slice[(i+1)*ny + j].real() : slice[(i+1)*ny + j].imag();
    float v01 = get_real ? slice[i*ny + (j+1)].real() : slice[i*ny + (j+1)].imag();
    float v11 = get_real ? slice[(i+1)*ny + (j+1)].real() : slice[(i+1)*ny + (j+1)].imag();

    // Pack into SIMD register
    __m128 vals = _mm_set_ps(v11, v01, v10, v00);

    // Compute weights
    __m128 weights = _mm_set_ps(
        x_rel * y_rel,
        (1-x_rel) * y_rel,
        x_rel * (1-y_rel),
        (1-x_rel) * (1-y_rel)
    );

    // Multiply and sum
    __m128 prod = _mm_mul_ps(vals, weights);
    __m128 sum1 = _mm_hadd_ps(prod, prod);
    __m128 sum2 = _mm_hadd_ps(sum1, sum1);

    return sum2;  // Contains result in lowest element
}

// Optimized CMF_search_3d with all three optimizations
complex<Vec> CMF_search_3d_optimized(
    float x_val, float y_val, float w_val,
    int nx, int ny,
    vector<float>& w_points,
    vector<complex<Vec>>& f
) {
    // ... (bounds checking, index computation) ...

    float dx = 1.0f / (nx - 1);
    float dy = 1.0f / (ny - 1);

    int i = x_val / dx;
    int j = y_val / dy;
    int k = binary_search(w_val, w_points);

    // Clamp
    if (i >= nx - 1) i = nx - 2;
    if (j >= ny - 1) j = ny - 2;
    if (k >= w_points.size() - 1) k = w_points.size() - 2;

    // Compute weights
    float x_rel = x_val / dx - i;
    float y_rel = y_val / dy - j;
    float w_rel = (w_val - w_points[k]) / (w_points[k+1] - w_points[k]);

    // W-FIRST: Get pointers to two frequency slices
    const complex<float>* slice_k  = &f[k * nx * ny];
    const complex<float>* slice_k1 = &f[(k+1) * nx * ny];

    // SIMD: Interpolate real parts
    __m128 real_k  = interpolate_2d_slice_simd(slice_k, i, j, ny, x_rel, y_rel, true);
    __m128 real_k1 = interpolate_2d_slice_simd(slice_k1, i, j, ny, x_rel, y_rel, true);

    // SIMD: Interpolate imaginary parts
    __m128 imag_k  = interpolate_2d_slice_simd(slice_k, i, j, ny, x_rel, y_rel, false);
    __m128 imag_k1 = interpolate_2d_slice_simd(slice_k1, i, j, ny, x_rel, y_rel, false);

    // Extract scalar results
    float val_k_real = _mm_cvtss_f32(real_k);
    float val_k1_real = _mm_cvtss_f32(real_k1);
    float val_k_imag = _mm_cvtss_f32(imag_k);
    float val_k1_imag = _mm_cvtss_f32(imag_k1);

    // Linear interpolation in frequency
    float result_real = (1 - w_rel) * val_k_real + w_rel * val_k1_real;
    float result_imag = (1 - w_rel) * val_k_imag + w_rel * val_k1_imag;

    return complex<Vec>(complex<float>(result_real, result_imag));
}

// Batch operator with PREFETCHING
vector<complex<float>> operator()(const vector<Vec>& points, float w) {
    vector<complex<float>> results(points.size());

    const int PREFETCH_DIST = 16;

    for (int i = 0; i < points.size(); i++) {
        // PREFETCH: Load future data
        if (i + PREFETCH_DIST < points.size()) {
            Vec future = transform_point(points[i + PREFETCH_DIST]);
            int k_future = find_frequency_index(w);

            int base = k_future * nx * ny;
            int idx = compute_spatial_idx(future);

            // Prefetch both frequency slices
            _mm_prefetch(&f[base + idx], _MM_HINT_T0);
            _mm_prefetch(&f[base + nx*ny + idx], _MM_HINT_T0);
        }

        // Process current point (W-FIRST + SIMD)
        results[i] = CMF_search_3d_optimized(
            points[i].x, points[i].y, w, nx, ny, w_points, f
        );
    }

    return results;
}
```

### Performance Prediction

| Optimization | Speedup | Cumulative |
|--------------|---------|------------|
| Baseline | 1.0× | 2.3M/sec |
| + W-first | 1.7× | 3.9M/sec |
| + SIMD | 2.0× | 7.8M/sec |
| + Prefetch | 1.3× | 10.1M/sec |

**Final: ~4.4× faster overall**

---

## Part 4: Practical Implementation Guide

### Step 1: Enable SIMD in Compiler

```bash
# Add to CMakeLists.txt or compile flags
-O3                  # Aggressive optimization
-march=native        # Use all CPU features
-mavx2               # Enable AVX2 explicitly
-ftree-vectorize     # Auto-vectorization
-ffast-math          # Relaxed FP math (be careful!)
```

### Step 2: Verify SIMD is Being Used

```bash
# Check assembly output
g++ -S -O3 -march=native field_funcs.cpp
grep "vmul" field_funcs.s   # Look for AVX instructions (vmul, vadd, etc.)

# Or use compiler reports
g++ -O3 -march=native -fopt-info-vec-all field_funcs.cpp 2>&1 | grep "vectorized"
```

### Step 3: Benchmark Each Optimization

```cpp
// Add timing to profile_interpolation.cpp
auto bench = [&](auto func, string name) {
    Timer t;
    for (int i = 0; i < 1000; i++) {
        func(test_points[i]);
    }
    cout << name << ": " << t.elapsed_ms() << " ms\n";
};

bench([&](Vec p) { return interpolate_scalar(p); }, "Scalar");
bench([&](Vec p) { return interpolate_wfirst(p); }, "W-first");
bench([&](Vec p) { return interpolate_simd(p); }, "SIMD");
bench([&](Vec p) { return interpolate_prefetch(p); }, "Prefetch");
```

### Step 4: Validate Correctness

```cpp
// Ensure SIMD gives same results as scalar
const float EPSILON = 1e-5;

for (int i = 0; i < 1000; i++) {
    auto scalar = interpolate_scalar(test_points[i]);
    auto simd = interpolate_simd(test_points[i]);

    float diff = abs(scalar - simd);
    if (diff > EPSILON) {
        cerr << "SIMD error at point " << i << ": diff = " << diff << "\n";
    }
}
```

---

## Summary

### SIMD
- **What:** Process multiple values with one instruction
- **How:** Use intrinsics (`_mm_add_ps`, etc.) or compiler auto-vectorization
- **Speedup:** 2-4× for computation-bound code
- **Best for:** Repeated arithmetic operations on arrays

### Prefetching
- **What:** Load data into cache before you need it
- **How:** `_mm_prefetch(&data[future_index], _MM_HINT_T0)`
- **Speedup:** 1.5-2× for memory-bound code
- **Best for:** Random/irregular memory access patterns

### Combined Strategy
1. **W-first:** Reorganize algorithm (1.7× speedup)
2. **SIMD:** Vectorize 2D interpolations (2× more)
3. **Prefetch:** Hide memory latency (1.3× more)

**Total: 4-5× faster than baseline!**

The beauty is they're **orthogonal** - each optimization addresses a different bottleneck, so they multiply together!
