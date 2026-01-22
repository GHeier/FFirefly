using Firefly
cfg = Firefly.Config

using HMatrices
using KrylovKit
using LinearAlgebra
using Printf
using StaticArrays

# Load relevant variables from the configuration
outdir = cfg.outdir
prefix = cfg.prefix
filetype = cfg.filetype

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
    test_k = [0.1, 0.2]
    diff = V(test_k, 0.1) - V(test_k, -0.1)
    if abs(diff / V(test_k, 0.1)) > 1e-2
        s = -1
        println("Warning: Vertex function is not even in frequency!")
        println("V(k, w) - V(k, -w) at test point = $diff")
        println("Percent difference: $(abs(diff / V(test_k, 0.1)) * 100)%")
    end
    # Compute DOS weights: dA/v for each k-point
    # velocities has shape (npoints, nbands, 3), we want magnitude for band 1
    v_norms = [LinearAlgebra.norm(velocities[i, 1, :]) for i in 1:length(areas)]
    dos_weights = areas ./ v_norms

    # Return a frequency-dependent kernel function V(w-w', k-k')
    # Signature: kernel(k1, k2, w1, w2, i, j, iw, jw)
    # Weighting: sqrt(dA_i/v_i * f(w1)) * sqrt(dA_j/v_j * f(w2)) * V(dk, dw)
    return function(k1, k2, w1, w2, i, j, iw, jw)
        # Convert to Float64 for Vertex call
        k1 = Float64.(k1)
        k2 = Float64.(k2)
        w1 = Float64(w1)
        w2 = Float64(w2)

        # Projected even-momentum vertex
        Vp = (V(k1 - k2, 0) + V(k1 + k2, 0)) / 2
        Vm = (V(k1 - k2, 0) + V(k1 + k2, 0)) / 2
        if abs(w2) > 1e-4
            V_val = ( f(-w2) * Vp - f(w2) * Vm ) / (2 * w2)
            #V_val = -(f(-w2) - f(w2)) / (2 * w2)
        else
            #V_val = 1 / (2 * T) * ( V(k1 - k2, 0) + V(k1 + k2, 0) )
            V_val = 1 / (4 * T) * ( Vm )
            V_val = -1 / (4 * T)
        end

        V_val = real(V_val)
        weight = dos_weights[j]

        return -wc / w_pts * weight * V_val
        #return -weight1 * weight2 * V_val
        return - (Z * (2 * wc / w_pts ) / 2 * w2) * weight * V_val
        #return - (Z * (2 * wc / w_pts ) / 2 * w2) * weight1 * weight2 * V_val
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
    Yclt = Xclt  # Same cluster tree for row/column spaces

    # Create compression method
    comp = HMatrices.PartialACA(; atol=atol, rank=rank)

    # Assemble HMatrix from matrix and cluster trees
    H = assemble_hmatrix(K, Xclt, Yclt; comp=comp)

    return H
end

function run()
    # Main function call goes here
    println("="^60)
    println("Initializing HMatrix + Arnoldi Eigenvalue Solver")
    println("="^60)

    # Parameters
    mu = -1.0  # Fermi energy
    println("\nGenerating k-points from Fermi surface at μ = $mu")
    H = Firefly.Hamiltonian()

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
    velocities = get_fermi_velocity(H, kpoints_f64)
    v_norms = [LinearAlgebra.norm(velocities[i, 1, :]) for i in 1:n]
    @printf("Fermi Velocity Min, Max, Ave: %.6f, %.6f, %.6f\n", minimum(v_norms), maximum(v_norms), sum(v_norms) / n)

    # Compute density of states weights
    dos_weights = areas ./ v_norms ./ (2 * π)^dim
    @printf("DOS Min, Max, Ave, Total: %.6f, %.6f, %.6f, %.6f\n", minimum(dos_weights), maximum(dos_weights), sum(dos_weights) / n, sum(dos_weights))

    # Temporary integral for testing
    ival = 0
    for i in 1:w_pts
        if w_points[i] == 0.0
            ival += 1 / (4 * T) 
        else
            ival += tanh(w_points[i] / (2 * T)) / (2 * w_points[i])
        end
    end
    println("Test integral ival over frequencies: ", ival * (wc / w_pts))
    println("Eigenvalue estimate: ", ival * (wc / w_pts) * sum(dos_weights))

    # Load Vertex
    println("\nLoading Vertex...")
    V = Firefly.Vertex()
    println("Vertex loaded.")

    # Create kernel function using Vertex with frequency-dependent weighting
    # Pattern: sqrt(dA_i/v_i * f(w1)) * sqrt(dA_j/v_j * f(w2)) * V(k_i - k_j, w1 - w2)
    kernel = make_vertex_kernel(V, areas, velocities, w_points, 1.0)

    # Build HMatrix with frequency dependence
    println("\n" * "="^60)
    println("Building HMatrix with frequency dependence...")
    t_hmat = @elapsed H_matrix = build_hmatrix(kpoints, w_points, kernel; atol=1e-4, rank=30)
    @printf("HMatrix built in %.3f seconds\n", t_hmat)
    @printf("Compression ratio: %.2f\n", compression_ratio(H_matrix))

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
        vals_hmat, vecs_hmat, info_hmat = eigsolve(hmv, n_total, num_eigs, :LR;
                                                    issymmetric=false,
                                                    krylovdim=30,
                                                    maxiter=200,
                                                    tol=1e-8)
    end

    # Find the largest positive eigenvalue
    # Note: Arnoldi can return complex eigenvalues for non-Hermitian matrices
    # Check for significant imaginary components
    imag_vals = imag.(vals_hmat)
    max_imag = maximum(abs.(imag_vals))
    println("Max imaginary component of eigenvalues: ", max_imag)

    real_vals = real.(vals_hmat) ./ (2π)^dim
    positive_vals = filter(x -> x > 0, real_vals)

    if isempty(positive_vals)
        @printf("Warning: No positive eigenvalues found!\n")
        @printf("All %d eigenvalues:\n", length(real_vals))
        for (i, val) in enumerate(real_vals)
            @printf("  λ_%d = %.10f\n", i, val)
        end
        largest_positive = NaN
        idx_largest_positive = 0
    else
        largest_positive = maximum(positive_vals)
        # Find the index of the largest positive eigenvalue
        idx_largest_positive = findfirst(x -> x == largest_positive, real_vals)
    end

    @printf("Time: %.3f seconds\n", t_lanczos_hmat)
    @printf("Number of eigenvalues computed: %d\n", length(vals_hmat))
    @printf("Converged: %s (iterations: %d)\n", info_hmat.converged > 0, info_hmat.numiter)

    # Show top eigenvalues
    @printf("\nTop eigenvalues:\n")
    for i in 1:min(5, length(real_vals))
        @printf("  λ_%d = %.10f\n", i, real_vals[i])
    end

    @printf("\nLargest eigenvalue (by magnitude): %.10f\n", real_vals[1])
    @printf("Largest positive eigenvalue: %.10f\n", largest_positive)
    println("Eigenvalue estimate: ", ival * (wc / w_pts) * sum(dos_weights))

    # Save the eigenvector corresponding to the largest positive eigenvalue
    if idx_largest_positive > 0
        eigenvec = vecs_hmat[idx_largest_positive]

        # Reshape eigenvector from 1D (n_total) to 2D (n_k, n_w)
        # Index ordering: idx = i_k + i_w * n_k (k varies faster)
        # Julia column-major: reshape fills columns first
        # So reshape(vec, n_k, n_w) will give us [k, w] indexing
        eigenvec_2d = reshape(eigenvec, n_k, n_w)

        # Save gap function with k-points from Fermi surface using save_data!
        filename = outdir * prefix * "_gap." * filetype
        @printf("\nSaving eigenvector to %s\n", filename)

        # Determine k-space dimension
        k_dim = length(kpoints[1])
        has_freq = (w_pts > 1)

        # Prepare k-points matrix [n_k × k_dim]
        points_matrix = Matrix{Float64}(undef, n_k, k_dim)
        for i_k in 1:n_k
            for d in 1:k_dim
                points_matrix[i_k, d] = Float64(kpoints[i_k][d])
            end
        end

        # Prepare data array
        # eigenvec_2d is already [n_k, n_w] from reshape above
        # Convert to real and Float64
        data_array = Float64.(real.(eigenvec_2d))

        # Call save_data! with points parameter
        save_data! = Firefly.Imports.save_data!
        if has_freq
            # With frequency: pass w_points separately
            save_data!(filename, data_array, points=points_matrix, w_points=w_points)
        else
            # No frequency: just spatial points
            save_data!(filename, data_array, points=points_matrix)
        end

        @printf("Eigenvector saved successfully.\n")
        @printf("  Shape: (%d k-points, %d frequencies)\n", n_k, n_w)
        @printf("  Total data points: %d\n", n_k * n_w)
        @printf("  Norm: %.6f\n", LinearAlgebra.norm(eigenvec))
    else
        @printf("\nNo positive eigenvalue found - eigenvector not saved.\n")
    end

    return largest_positive
end

if abspath(PROGRAM_FILE) == @__FILE__ # Runs on file execution
    run()
end



