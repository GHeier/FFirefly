using Firefly
cfg = Firefly.Config

using HMatrices
using KrylovKit
using LinearAlgebra
using Printf
using StaticArrays

# Load relevant variables from the configuration
dim = cfg.dimension
mu = cfg.fermi_energy
T = cfg.Temperature
w_pts = cfg.w_pts
wc = cfg.cutoff_energy


function load_surface()
    eps_func = (k) -> Firefly.epsilon(1, [Float64(k.x), Float64(k.y), Float64(k.z)])

    # Create Surface at the Fermi energy (mu)
    surf = Firefly.Surface(eps_func, Float32(mu))
    kpoints, areas = get_faces_and_areas(surf)
    println("Found $(length(kpoints)) k-points on Fermi surface")
    println("Total area: $(sum(areas))")

    return kpoints, areas
end

function f(w)
    return 1 / (exp(w / T) + 1)
end

function make_vertex_kernel(V::Firefly.Vertex, areas::Vector{Float32}, velocities::Array{Float32, 3}, frequencies::Vector{Float32}, Z)
    s = 1 # Even variable
    diff = V([0.1, 0.2, 0.3], 0.1) - V([0.1, 2.0, 3.0], -0.1)
    if abs(diff) > 1e-5
        s = -1
        println("Warning: Vertex function is not even in frequency!")
        println("V(k, w) - V(k, -w) at test point = $diff")
    end
    # Compute DOS weights: dA/v for each k-point
    # velocities has shape (npoints, nbands, 3), we want magnitude for band 1
    v_norms = [LinearAlgebra.norm(velocities[i, 1, :]) for i in 1:length(areas)]
    dos_weights = areas ./ v_norms

    # Return a frequency-dependent kernel function V(w-w', k-k')
    # Signature: kernel(k1, k2, w1, w2, i, j, iw, jw)
    # Weighting: sqrt(dA_i/v_i * f(w1)) * sqrt(dA_j/v_j * f(w2)) * V(dk, dw)
    return function(k1, k2, w1, w2, i, j, iw, jw)
        dk = k1 .- k2
        w1 = Float64(w1)
        w2 = Float64(w2)

        # Convert to Float64 for Vertex call
        dk_f64 = Float64.(dk)

        # Get vertex value at momentum and frequency transfer: V(k1-k2, w1-w2)
        # V returns ComplexF32, take real part for kernel matrix
        V_projected = real(V(dk_f64, w1 - w2)) - s * real(V(dk_f64, w1 + w2)) 
        V_val = Float64(V_projected)

        # Fermi function factors at both frequencies
        f1 = f(w1)
        f2 = f(w2)

        # Apply the frequency-dependent DOS weighting from both sides
        # Pattern: sqrt(dA_i/v_i * f(w1)) * sqrt(dA_j/v_j * f(w2)) * V
        weight1 = sqrt(dos_weights[i] * f1)
        weight2 = sqrt(dos_weights[j] * f2)

        return - (Z / 2 * w2) * weight1 * weight2 * V_val
    end
end

struct KernelMatrixWrapper{F, T} <: AbstractMatrix{Float64}
    kernel::F
    kpoints::Vector{T}
end

Base.getindex(K::KernelMatrixWrapper, i::Int, j::Int) = K.kernel(K.kpoints[i], K.kpoints[j], i, j)
Base.size(K::KernelMatrixWrapper) = (length(K.kpoints), length(K.kpoints))
Base.eltype(::Type{<:KernelMatrixWrapper}) = Float64

# Frequency-dependent kernel matrix wrapper for (k,w) pairs
struct FreqKernelMatrixWrapper{F, T} <: AbstractMatrix{Float64}
    kernel::F
    kpoints::Vector{T}
    frequencies::Vector{Float32}
    n_k::Int
    n_w::Int
end

function Base.getindex(K::FreqKernelMatrixWrapper, idx1::Int, idx2::Int)
    # Decompose composite indices into (k, w) pairs
    # Index ordering: k varies faster, so idx = i_k + i_w * n_k
    i_k = ((idx1 - 1) % K.n_k) + 1
    i_w = div(idx1 - 1, K.n_k) + 1
    j_k = ((idx2 - 1) % K.n_k) + 1
    j_w = div(idx2 - 1, K.n_k) + 1

    k1 = K.kpoints[i_k]
    k2 = K.kpoints[j_k]
    w1 = K.frequencies[i_w]
    w2 = K.frequencies[j_w]

    return K.kernel(k1, k2, w1, w2, i_k, j_k, i_w, j_w)
end

Base.size(K::FreqKernelMatrixWrapper) = (K.n_k * K.n_w, K.n_k * K.n_w)
Base.eltype(::Type{<:FreqKernelMatrixWrapper}) = Float64

# Build HMatrix with frequency dependence
function build_hmatrix(kpoints, frequencies, kernel_func; atol=1e-6, rank=20)
    n_k = length(kpoints)
    n_w = length(frequencies)
    n_total = n_k * n_w

    @printf("Building HMatrix for %d k-points × %d frequencies = %d total points\n", n_k, n_w, n_total)

    # Create frequency-dependent kernel matrix wrapper
    K = FreqKernelMatrixWrapper(kernel_func, kpoints, frequencies, n_k, n_w)

    # Build cluster tree based on whether we have frequency dependence
    k_dim = length(kpoints[1])

    if n_w == 1
        # Single frequency: use 3D clustering (k-space only)
        # This avoids degenerate 4D points when all frequencies are the same
        println("Single frequency detected, using 3D k-space clustering...")
        if k_dim == 2
            X = [SVector{3}(Float64(kp[1]), Float64(kp[2]), 0.0) for kp in kpoints]
        elseif k_dim == 3
            X = [SVector{3}(Float64(kp[1]), Float64(kp[2]), Float64(kp[3])) for kp in kpoints]
        else
            error("K-points must be 2D or 3D, got dimension $k_dim")
        end
    else
        # Multiple frequencies: use 4D clustering (k,w)-space
        # Index ordering: k varies faster (idx = i_k + i_w * n_k)
        X = Vector{SVector{4,Float64}}(undef, n_total)
        idx = 1
        for i_w in 1:n_w
            w = Float64(frequencies[i_w])
            for i_k in 1:n_k
                kp = kpoints[i_k]
                if k_dim == 2
                    # 2D k-space: use (kx, ky, 0, w)
                    X[idx] = SVector{4}(Float64(kp[1]), Float64(kp[2]), 0.0, w)
                elseif k_dim == 3
                    # 3D k-space: use (kx, ky, kz, w)
                    X[idx] = SVector{4}(Float64(kp[1]), Float64(kp[2]), Float64(kp[3]), w)
                else
                    error("K-points must be 2D or 3D, got dimension $k_dim")
                end
                idx += 1
            end
        end
    end

    # Build the hierarchical matrix
    println("Building HMatrix...")
    splitter = HMatrices.CardinalitySplitter(; nmax=50)
    Xclt = HMatrices.ClusterTree(X, splitter)
    Yclt = Xclt  # Same cluster tree for symmetric matrix

    # Create compression method
    comp = HMatrices.PartialACA(; atol=atol, rank=rank)

    # Assemble HMatrix from matrix and cluster trees
    H = assemble_hmatrix(K, Xclt, Yclt; comp=comp)

    return H
end

function run()
    # Main function call goes here
    println("="^60)
    println("Initializing HMatrix + Lanczos Eigenvalue Solver")
    println("="^60)

    # Parameters
    mu = -1.0  # Fermi energy
    println("\nGenerating k-points from Fermi surface at μ = $mu")
    H = Firefly.Hamiltonian()
    println("Testing epsilon at a few k-points:")
    for ktest in [[0.0, 0.0, 0.0], [π, π, π], [0.5, 0.5, 0.5]]
        ε = Firefly.epsilon(1, ktest)
        println("  ε($(ktest)) = $ε")
    end

    w_points = Vector{Float32}(undef, w_pts)
    if w_pts == 1
        # Single frequency point: use w=0 (static pairing)
        w_points[1] = 0.0
    else
        # Multiple frequency points: distribute from 0 to wc
        for i in 1:w_pts
            #w_points[i] = -wc + (2 * wc) * (i - 1) / (w_pts - 1)
            w_points[i] = wc * (i - 1) / (w_pts - 1)
        end
    end

    kpoints, areas = load_surface()
    n = length(kpoints)
    println("Total k-points: $n")

    # Load Fermi velocities at k-points
    println("\nComputing Fermi velocities...")
    get_fermi_velocity = Firefly.Imports.get_fermi_velocity
    # Convert Float32 k-points to Float64 for get_fermi_velocity
    kpoints_f64 = [Float64.(kp) for kp in kpoints]
    t_vel = @elapsed velocities = get_fermi_velocity(H, kpoints_f64)
    @printf("Fermi velocities computed in %.3f seconds\n", t_vel)

    # Statistics on velocities
    v_norms = [LinearAlgebra.norm(velocities[i, 1, :]) for i in 1:n]
    @printf("Velocity statistics:\n")
    @printf("  Min |v|: %.6f\n", minimum(v_norms))
    @printf("  Max |v|: %.6f\n", maximum(v_norms))
    @printf("  Mean |v|: %.6f\n", sum(v_norms) / n)

    # Compute density of states weights
    dos_weights = areas ./ v_norms
    @printf("\nDensity of states (dA/v) statistics:\n")
    @printf("  Min dA/v: %.6f\n", minimum(dos_weights))
    @printf("  Max dA/v: %.6f\n", maximum(dos_weights))
    @printf("  Mean dA/v: %.6f\n", sum(dos_weights) / n)
    @printf("  Total DOS: %.6f\n", sum(dos_weights))

    # Load Vertex
    println("\nLoading Vertex...")
    V = Firefly.Vertex()
    println("Vertex loaded.")

    # Create kernel function using Vertex with frequency-dependent weighting
    # Pattern: sqrt(dA_i/v_i * f(w1)) * sqrt(dA_j/v_j * f(w2)) * V(k_i - k_j, w1 - w2)
    kernel = make_vertex_kernel(V, areas, velocities, w_points, 1.0)

    # Build HMatrix with frequency dependence
    println("\n" * "="^60)
    println("Building frequency-dependent pairing matrix: V(w-w', k-k')")
    t_hmat = @elapsed H_matrix = build_hmatrix(kpoints, w_points, kernel; atol=1e-4, rank=30)
    @printf("HMatrix built in %.3f seconds\n", t_hmat)
    @printf("Compression ratio: %.2f\n", compression_ratio(H_matrix))
    @printf("Memory savings: %.2f%%\n", (1 - compression_ratio(H_matrix)) * 100)

    # Solve with HMatrix using Lanczos
    println("\n" * "="^60)
    println("Solving HMatrix with Lanczos (KrylovKit)...")

    # Matrix dimensions: n_k k-points × n_w frequencies
    n_k = length(kpoints)
    n_w = length(w_points)
    n_total = n_k * n_w
    @printf("Matrix size: %d × %d\n", n_total, n_total)

    # Create a function that applies H_matrix to a vector using mul!
    # HMatrices.jl supports matrix-vector multiplication via mul!
    hmv = (v) -> begin
        result = similar(v)
        mul!(result, H_matrix, v)
        return result
    end

    # Use the matrix-vector product function with eigsolve
    # Request more eigenvalues to find the largest positive one
    num_eigs = min(10, n_total)  # Request up to 10 eigenvalues
    t_lanczos_hmat = @elapsed begin
        vals_hmat, vecs_hmat, info_hmat = eigsolve(hmv, n_total, num_eigs, :LM;
                                                    issymmetric=true,
                                                    krylovdim=30,
                                                    maxiter=200,
                                                    tol=1e-8)
    end

    # Find the largest positive eigenvalue
    real_vals = real.(vals_hmat) ./ (2π)^dim 
    positive_vals = filter(x -> x > 0, real_vals)

    if isempty(positive_vals)
        @printf("Warning: No positive eigenvalues found!\n")
        @printf("All %d eigenvalues:\n", length(real_vals))
        for (i, val) in enumerate(real_vals)
            @printf("  λ_%d = %.10f\n", i, val)
        end
        largest_positive = NaN
    else
        largest_positive = maximum(positive_vals)
    end

    @printf("Time: %.3f seconds\n", t_lanczos_hmat)
    @printf("Number of eigenvalues computed: %d\n", length(vals_hmat))
    @printf("Converged: %s (iterations: %d)\n", info_hmat.converged > 0, info_hmat.numiter)

    dos = 0.0
    for i in 1:n
        dos += areas[i] / LinearAlgebra.norm(velocities[i, 1, :])
    end
    dos /= (4 * π^2)
    @printf("\nEstimated Density of States at Fermi level: %.6f states/eV/unit cell\n", dos)

    # Show top eigenvalues
    @printf("\nTop eigenvalues:\n")
    for i in 1:min(5, length(real_vals))
        @printf("  λ_%d = %.10f\n", i, real_vals[i])
    end

    @printf("\nLargest eigenvalue (by magnitude): %.10f\n", real_vals[1])
    @printf("Largest positive eigenvalue: %.10f\n", largest_positive)

    return largest_positive
end

if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    run()
end



