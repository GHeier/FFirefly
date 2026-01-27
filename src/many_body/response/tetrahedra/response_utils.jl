"""
Utility functions for tetrahedron-method response function calculation.

This module contains helper functions for:
- K-mesh generation and manipulation
- Energy evaluation on shifted grids E(k+q)
- DOS calculation using tetrahedron method
- Individual response function calculation
"""

using BZIntegral
using BZIntegral.BZInt2D
using Interpolations
using Printf

"""
    get_fractional_mesh(n; centered=true)

Generate fractional coordinate mesh.

# Arguments
- `n`: Tuple (nx, ny, nz) of mesh dimensions
- `centered`: If true, center mesh at [-0.5, 0.5), else [0, 1)

# Returns
- Matrix of size (nkpts, 3) with fractional coordinates
"""
function get_fractional_mesh(n; centered=true)
    nx, ny, nz = n

    if centered
        xs = (0:nx-1) ./ nx .- 0.5
        ys = (0:ny-1) ./ ny .- 0.5
        zs = (0:nz-1) ./ nz .- 0.5
    else
        xs = (0:nx-1) ./ nx
        ys = (0:ny-1) ./ ny
        zs = (0:nz-1) ./ nz
    end

    # fractional k-points: 3 × Nk
    frac = hcat(vec(repeat(xs, inner=(ny*nz))),
                vec(repeat(ys, inner=(nz,), outer=(nx))),
                vec(repeat(zs, outer=(nx*ny))))'
    # Return as Nk × 3 matrix (transpose from 3 × Nk)
    return Matrix{Float64}(frac')
end

"""
    get_kmesh(BZ::AbstractMatrix{<:Real}, n; centered=true)

Generate Cartesian k-point mesh from Brillouin zone vectors.

# Arguments
- `BZ`: Brillouin zone matrix (3×3)
- `n`: Tuple (nx, ny, nz) of mesh dimensions
- `centered`: If true, center mesh at [-0.5, 0.5), else [0, 1)

# Returns
- Matrix of size (nkpts, 3) with Cartesian k-coordinates
"""
function get_kmesh(BZ::AbstractMatrix{<:Real}, n; centered=true)
    nx, ny, nz = n

    if centered
        xs = (0:nx-1) ./ nx .- 0.5
        ys = (0:ny-1) ./ ny .- 0.5
        zs = (0:nz-1) ./ nz .- 0.5
    else
        xs = (0:nx-1) ./ nx
        ys = (0:ny-1) ./ ny
        zs = (0:nz-1) ./ nz
    end

    # fractional k-points: 3 × Nk
    frac = hcat(vec(repeat(xs, inner=(ny*nz))),
                vec(repeat(ys, inner=(nz,), outer=(nx))),
                vec(repeat(zs, outer=(nx*ny))))'
    # Cartesian k-points: 3 × Nk
    kpts = BZ * frac
    # Return as Nk × 3 matrix (transpose from 3 × Nk)
    return Matrix{Float64}(kpts')
end

"""
    evaluate_Ekq_on_grid(itp_Ek, kmesh, q_frac, dim=3)

Evaluate E(k+q) for all k on the grid with periodic boundary conditions.

# Arguments
- `itp_Ek`: Interpolation object from Interpolations.jl
- `kmesh`: (nkx, nky, nkz) grid dimensions
- `q_frac`: [qx, qy, qz] momentum transfer in fractional coordinates
- `dim`: Spatial dimension (2 or 3)

# Returns
- Array of E(k+q) values on the k-grid with shape (nkx, nky, nkz) or (nkx, nky)
"""
function evaluate_Ekq_on_grid(itp_Ek, kmesh, q_frac, dim=3)
    nkx, nky, nkz = kmesh

    # Original k-grid in fractional coordinates (centered: -0.5 to 0.5)
    kx_frac = (0:nkx-1) ./ nkx .- 0.5
    ky_frac = (0:nky-1) ./ nky .- 0.5
    kz_frac = (0:nkz-1) ./ nkz .- 0.5

    # Add q shift and apply periodic boundary conditions: wrap to [-0.5, 0.5)
    kxq_frac = mod.(kx_frac .+ q_frac[1] .+ 0.5, 1.0) .- 0.5
    kyq_frac = mod.(ky_frac .+ q_frac[2] .+ 0.5, 1.0) .- 0.5
    kzq_frac = mod.(kz_frac .+ q_frac[3] .+ 0.5, 1.0) .- 0.5

    # Convert to interpolation indices (1-based)
    Ikxq = (kxq_frac .+ 0.5) .* nkx .+ 1
    Ikyq = (kyq_frac .+ 0.5) .* nky .+ 1
    Ikzq = (kzq_frac .+ 0.5) .* nkz .+ 1

    # Evaluate interpolation with broadcasting
    if dim == 2
        Ix = repeat(reshape(Ikxq, :, 1), 1, nky)
        Iy = repeat(reshape(Ikyq, 1, :), nkx, 1)
        return itp_Ek.(Ix, Iy)
    else
        Ix = repeat(reshape(Ikxq, :, 1, 1), 1, nky, nkz)
        Iy = repeat(reshape(Ikyq, 1, :, 1), nkx, 1, nkz)
        Iz = repeat(reshape(Ikzq, 1, 1, :), nkx, nky, 1)
        return itp_Ek.(Ix, Iy, Iz)
    end
end

"""
    calculate_dos(E, Ek_mesh, iter=2)

Calculate density of states at energy E using tetrahedron method.

# Arguments
- `E`: Energy at which to calculate DOS
- `Ek_mesh`: Energy values on k-mesh
- `iter`: Iteration parameter for tetrahedron integration

# Returns
- DOS value at energy E
"""
function calculate_dos(E, Ek_mesh, iter=2)
    dim = ndims(Ek_mesh)

    if dim == 2
        Wmesh = Quad2DRuleδ(Ek_mesh, E, iter)
    elseif dim == 3
        Wmesh = Quad3DRuleδ(Ek_mesh, E, iter)
    else
        error("Dimension must be 2 or 3")
    end
    DOS = sum(Wmesh)

    return DOS
end

"""
    calculate_response_bzintegral_2(w, Ek_mesh, Ekq_mesh, mu, iter=2)

Calculate response function using BZIntegral tetrahedron method.

χ₀(q,ω) = -2 ∫dk [f(ε(k+q)) - f(ε(k))] / [ω + ε(k) - ε(k+q) + iη]

For the integrand F(k) = 1/[ω + ε(k) - ε(k+q) + iη], we compute:
  ∫ Θ(-ε(k+q)) * F(k) dk  -  ∫ Θ(-ε(k)) * F(k) dk

# Arguments
- `w`: Frequency
- `Ek_mesh`: E(k) on k-mesh
- `Ekq_mesh`: E(k+q) on k-mesh
- `mu`: Chemical potential
- `iter`: Iteration parameter for tetrahedron integration

# Returns
- Complex susceptibility value at (q, ω)
"""
function calculate_response_bzintegral_2(w, Ek_mesh, Ekq_mesh, mu, iter=2)
    # Denominator: ω + ε(k) - ε(k+q) + small eta to avoid division by zero
    eta = 1e-4  # Small broadening parameter
    denom = w .+ Ek_mesh .- Ekq_mesh .+ eta

    dim = ndims(Ek_mesh)

    if dim == 2
        Wmesh = Quad2DRuleΘ𝔇(Ek_mesh, mu, denom, iter) - Quad2DRuleΘ𝔇(Ekq_mesh, mu, denom, iter)
    elseif dim == 3
        Wmesh = Quad3DRuleΘ𝔇(Ek_mesh, mu, denom, iter) - Quad3DRuleΘ𝔇(Ekq_mesh, mu, denom, iter)
    else
        error("Dimension must be 2 or 3")
    end

    out = -1.0 * sum(Wmesh)  # Don't include -2 spin factor
    return out
end

"""
    print_progress(current, total, start_time)

Print progress bar with ETA.

# Arguments
- `current`: Current iteration count
- `total`: Total iterations
- `start_time`: Time when computation started
"""
function print_progress(current, total, start_time)
    elapsed = time() - start_time
    rate = current / elapsed
    remaining = (total - current) / rate
    if remaining > 60
        remaining /= 60
        @printf("  Progress: %d/%d | Rate: %.1f/s | ETA: %.1f min        \r", current, total, rate, remaining)
    else
        @printf("  Progress: %d/%d | Rate: %.1f/s | ETA: %.1f s           \r", current, total, rate, remaining)
    end
    flush(stdout)
end
