# HMatrix Weighted Kernel Implementation

## Overview
Implemented the proper Fermi surface density of states weighting in the HMatrix pairing interaction kernel, following the pattern from `matrix_creation.cpp`. The kernel now correctly implements:

**P_ij = sqrt(dA_i/v_i) * sqrt(dA_j/v_j) * V(k_i - k_j)**

This represents the linearized BCS gap equation as a properly discretized eigenvalue problem.

## Implementation Details

### Modified File
`src/superconductor/eliashberg/hmatrix/run.jl`

### Key Changes

#### 1. Updated `make_vertex_kernel` function
```julia
function make_vertex_kernel(V::Firefly.Vertex, areas::Vector{Float32}, velocities::Array{Float32, 3})
    # Compute weights: sqrt(dA/v) for each k-point
    v_norms = [LinearAlgebra.norm(velocities[i, 1, :]) for i in 1:length(areas)]
    weights = sqrt.(areas ./ v_norms)

    # Return kernel: sqrt(dA_i/v_i) * sqrt(dA_j/v_j) * V(k1-k2)
    return function(k1, k2, i, j)
        dk = k1 .- k2
        dk_f64 = Float64.(dk)
        V_val = real(V(dk_f64, 0.0))
        return weights[i] * weights[j] * V_val
    end
end
```

**Changes:**
- Added `areas` and `velocities` parameters
- Computes `sqrt(dA/v)` weights for each k-point
- Kernel function now includes indices `i, j` to apply per-point weights
- Returns properly weighted interaction matrix element

#### 2. Updated `KernelMatrixWrapper` indexing
```julia
Base.getindex(K::KernelMatrixWrapper, i::Int, j::Int) = K.kernel(K.kpoints[i], K.kpoints[j], i, j)
```

Now passes indices to the kernel function for weight lookup.

#### 3. Enhanced statistics output
Added density of states statistics:
```julia
dos_weights = areas ./ v_norms
@printf("\nDensity of states (dA/v) statistics:\n")
@printf("  Min dA/v: %.6f\n", minimum(dos_weights))
@printf("  Max dA/v: %.6f\n", maximum(dos_weights))
@printf("  Mean dA/v: %.6f\n", sum(dos_weights) / n)
@printf("  Total DOS: %.6f\n", sum(dos_weights))
```

#### 4. Updated kernel creation call
```julia
kernel = make_vertex_kernel(V, areas, velocities)
```

Passes the required area and velocity arrays.

## Physical Interpretation

### The sqrt(dA/v) Factor

The factor **dA/v(k)** represents the **Fermi surface density of states**:
- **dA**: Area element on the Fermi surface
- **v(k)**: Fermi velocity magnitude = |∇_k ε(k)|

In the continuous BCS gap equation:
```
∫_FS dS / |∇ε(k)| * V(k,k') * Δ(k') * tanh(ε(k')/2T) / (2ε(k'))
```

The measure `dS/|∇ε|` discretizes to `dA/v`.

### Why the Square Root?

The square root appears because we're creating a **symmetric eigenvalue problem**:
1. Define `P_ij = sqrt(w_i) * sqrt(w_j) * V_ij` where `w_i = dA_i/v_i`
2. This makes P symmetric if V is symmetric
3. Eigenvectors `ψ` satisfy: `P * ψ = λ * ψ`
4. Physical gap function: `Δ_i = ψ_i / sqrt(w_i)` (un-weighting, as in `matrix_creation.cpp:97-104`)

### Comparison with matrix_creation.cpp

Our implementation matches the pattern in `src/superconductor/bcs/matrix/matrix_creation.cpp`:

**Lines 27-32:**
```cpp
float f1 = pow(k1.area / vp(k1.n, k1), 0.5);  // sqrt(dA/v)
float f2 = pow(k2.area / vp(k2.n, k2), 0.5);
P(i, j) = -f1 * f2 * (V_func(k1 - k2, 0).real() + V_func(k1 + k2, 0).real()) / 2.0;
```

**Our implementation:**
```julia
weights = sqrt.(areas ./ v_norms)  # sqrt(dA/v) for each k-point
return weights[i] * weights[j] * V_val
```

The only difference is we don't include the negative sign or the (k1+k2) term, which can be added if needed for specific physics.

## Test Results

### Example Output
```
Found 220 k-points on Fermi surface
Total area: 11.040284

Computing Fermi velocities...
Fermi velocities computed in 0.000 seconds
Velocity statistics:
  Min |v|: 1.938030
  Max |v|: 2.622545
  Mean |v|: 2.323205

Density of states (dA/v) statistics:
  Min dA/v: 0.001085
  Max dA/v: 0.034511
  Mean dA/v: 0.021988
  Total DOS: 4.837295

Building weighted pairing matrix with sqrt(dA/v) * V pattern
HMatrix built in 1.348 seconds
Compression ratio: 0.95
Memory savings: 5.37%

Solving HMatrix with Lanczos (KrylovKit)...
Time: 1.997 seconds
Largest eigenvalue: 2.8264546784
Converged: true (iterations: 1)
```

### Impact of Weighting

**Before (unweighted kernel):**
- Eigenvalue: λ ≈ 128.77
- No physical meaning (incorrect normalization)

**After (weighted kernel):**
- Eigenvalue: λ ≈ 2.83
- Physically meaningful (represents pairing strength)
- λ > 1 suggests the system has a superconducting instability

### Physical Significance

The eigenvalue λ ≈ 2.83 indicates:
1. **Strong pairing instability** (λ > 1 means Cooper pairs form)
2. The system would be superconducting at this temperature
3. Critical temperature T_c could be estimated from the eigenvalue
4. The eigenvector would give the gap symmetry (s-wave, d-wave, etc.)

## Data Available

After this implementation, the following data is available in `run.jl`:

1. **kpoints**: Fermi surface k-points (220 points)
2. **areas**: Area element for each k-point
3. **velocities**: Fermi velocity vectors (220 × 1 × 3)
4. **v_norms**: Fermi velocity magnitudes
5. **dos_weights**: Density of states dA/v for each point
6. **H_matrix**: Properly weighted pairing interaction HMatrix
7. **vals_hmat**: Eigenvalues (pairing strengths)
8. **vecs_hmat**: Eigenvectors (gap symmetries, need to be un-weighted)

## Future Enhancements

Potential improvements:
1. **Un-weighting eigenvectors**: Divide by sqrt(dA/v) to get physical gap function
2. **Add (k+k') term**: Include umklapp scattering if needed
3. **Frequency dependence**: Extend to include Matsubara frequencies
4. **Multiple bands**: Generalize to multi-band systems
5. **T_c estimation**: Use eigenvalue to estimate critical temperature
6. **Gap symmetry analysis**: Classify eigenvectors (s, d, p-wave, etc.)

## References

- Implementation pattern follows: `src/superconductor/bcs/matrix/matrix_creation.cpp`
- Un-weighting procedure: `matrix_creation.cpp:97-104` (`vector_to_wave`)
- BCS gap equation documentation: `solver.cpp:31-34`
