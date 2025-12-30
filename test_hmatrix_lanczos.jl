using HMatrices
using KrylovKit
using LinearAlgebra
using Printf
using StaticArrays

# Add FFirefly Julia module to path and import
using Firefly

# Define a kernel function f(k-k')
# Example: Coulomb-like interaction in momentum space
function kernel(k1, k2)
    dk = k1 .- k2
    q = LinearAlgebra.norm(dk)
    # Regularized Coulomb: 1/(q^2 + λ^2)
    λ = 0.1
    return 1.0 / (q^2 + λ^2)
end

# Generate k-points using the Surface object with epsilon function
function generate_kpoints_from_surface(mu::Float64)
    # Create the epsilon function for band 1
    # The callback receives RawVec (the struct), not Vec (the wrapper)
    eps_func = (k_raw) -> Firefly.epsilon(1, [Float64(k_raw.x), Float64(k_raw.y), Float64(k_raw.z)])

    # Create Surface at the Fermi energy (mu)
    println("Creating Fermi surface at μ = $mu using tetrahedron method...")
    surf = Firefly.Surface(eps_func, Float32(mu))
    faces = get_faces(surf)

    println("Found $(length(faces)) k-points on Fermi surface")

    # Convert Vec objects to regular vectors
    return faces
end

# Build dense matrix for comparison
function build_dense_matrix(kpoints, kernel_func)
    n = length(kpoints)
    M = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            M[i,j] = kernel_func(kpoints[i], kpoints[j])
        end
    end
    return M
end

# Wrapper struct to make kernel function indexable
struct KernelMatrixWrapper{F, T} <: AbstractMatrix{Float64}
    kernel::F
    kpoints::Vector{T}
end

Base.getindex(K::KernelMatrixWrapper, i::Int, j::Int) = K.kernel(K.kpoints[i], K.kpoints[j])
Base.size(K::KernelMatrixWrapper) = (length(K.kpoints), length(K.kpoints))
Base.eltype(::Type{<:KernelMatrixWrapper}) = Float64

# Build HMatrix
function build_hmatrix(kpoints, kernel_func; atol=1e-6, rank=20)
    n = length(kpoints)

    # Convert kpoints to vector of SVectors for HMatrices
    X = [SVector{3}(kp) for kp in kpoints]

    # Create an indexable kernel matrix wrapper
    K = KernelMatrixWrapper(kernel_func, kpoints)

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

function main()
    println("="^60)
    println("HMatrix + Lanczos Eigenvalue Solver Test")
    println("="^60)

    # Parameters
    mu = -1.0  # Fermi energy
    println("\nGenerating k-points from Fermi surface at μ = $mu")
    println("Testing epsilon at a few k-points:")
    for ktest in [[0.0, 0.0, 0.0], [π, π, π], [0.5, 0.5, 0.5]]
        ε = Firefly.epsilon(1, ktest)
        println("  ε($(ktest)) = $ε")
    end

    kpoints = generate_kpoints_from_surface(mu)
    n = length(kpoints)
    println("Total k-points: $n")

    if n == 0
        @warn "No k-points found on Fermi surface at μ=$mu."
        println("The Fermi surface might be outside the range of the tetrahedron mesh.")
        println("Falling back to a small grid of k-points for demonstration...")
        # Create a small manual grid for demonstration
        kpoints = Vector{Vector{Float64}}()
        for i in 1:8
            for j in 1:8
                for k in 1:8
                    kx = 2π * (i-1) / 8
                    ky = 2π * (j-1) / 8
                    kz = 2π * (k-1) / 8
                    push!(kpoints, [kx, ky, kz])
                end
            end
        end
        n = length(kpoints)
        println("Using $(n) grid k-points instead")
    end

    # Build dense matrix
    println("\n" * "="^60)
    println("Building dense matrix...")
    t_dense = @elapsed M_dense = build_dense_matrix(kpoints, kernel)
    @printf("Dense matrix built in %.3f seconds\n", t_dense)
    @printf("Matrix size: %d x %d\n", size(M_dense)...)
    @printf("Memory: %.2f MB\n", sizeof(M_dense) / 1024^2)

    # Build HMatrix
    println("\n" * "="^60)
    t_hmat = @elapsed H = build_hmatrix(kpoints, kernel; atol=1e-4, rank=30)
    @printf("HMatrix built in %.3f seconds\n", t_hmat)
    @printf("Compression ratio: %.2f\n", compression_ratio(H))
    @printf("Memory savings: %.2f%%\n", (1 - compression_ratio(H)) * 100)

    # Solve with dense matrix using standard eigen
    println("\n" * "="^60)
    println("Solving with dense matrix (full eigendecomposition)...")
    t_eigen = @elapsed begin
        evals_dense, evecs_dense = eigen(Symmetric(M_dense))
    end
    max_eval_dense = evals_dense[end]
    max_evec_dense = evecs_dense[:, end]
    @printf("Time: %.3f seconds\n", t_eigen)
    @printf("Largest eigenvalue: %.10f\n", max_eval_dense)

    # Solve with dense matrix using Lanczos (via KrylovKit)
    println("\n" * "="^60)
    println("Solving dense matrix with Lanczos (KrylovKit)...")
    t_lanczos_dense = @elapsed begin
        vals_dense_kr, vecs_dense_kr, info = eigsolve(M_dense, 1, :LM;
                                                       krylovdim=30,
                                                       maxiter=200,
                                                       tol=1e-8)
    end
    @printf("Time: %.3f seconds\n", t_lanczos_dense)
    @printf("Largest eigenvalue: %.10f\n", real(vals_dense_kr[1]))
    @printf("Converged: %s (iterations: %d)\n", info.converged > 0, info.numiter)

    # Solve with HMatrix using Lanczos
    println("\n" * "="^60)
    println("Solving HMatrix with Lanczos (KrylovKit)...")

    # Create a function that applies H to a vector using mul!
    # HMatrices.jl supports matrix-vector multiplication via mul!
    hmv = (v) -> begin
        result = similar(v)
        mul!(result, H, v)
        return result
    end

    # Use the matrix-vector product function with eigsolve
    t_lanczos_hmat = @elapsed begin
        vals_hmat, vecs_hmat, info_hmat = eigsolve(hmv, n, 1, :LM;
                                                    issymmetric=true,
                                                    krylovdim=30,
                                                    maxiter=200,
                                                    tol=1e-8)
    end
    @printf("Time: %.3f seconds\n", t_lanczos_hmat)
    @printf("Largest eigenvalue: %.10f\n", real(vals_hmat[1]))
    @printf("Converged: %s (iterations: %d)\n", info_hmat.converged > 0, info_hmat.numiter)

    # Compare results
    println("\n" * "="^60)
    println("Comparison of Results:")
    println("="^60)
    @printf("Dense (eigen):           λ_max = %.10f\n", max_eval_dense)
    @printf("Dense (Lanczos):         λ_max = %.10f\n", real(vals_dense_kr[1]))
    @printf("HMatrix (Lanczos):       λ_max = %.10f\n", real(vals_hmat[1]))
    @printf("\nRelative errors:\n")
    @printf("  Dense Lanczos vs eigen:  %.2e\n",
            abs(real(vals_dense_kr[1]) - max_eval_dense) / abs(max_eval_dense))
    @printf("  HMatrix vs eigen:        %.2e\n",
            abs(real(vals_hmat[1]) - max_eval_dense) / abs(max_eval_dense))

    # Compare eigenvectors (normalize and check overlap)
    evec_dense_kr_normalized = real.(vecs_dense_kr[1]) / LinearAlgebra.norm(real.(vecs_dense_kr[1]))
    evec_hmat_normalized = real.(vecs_hmat[1]) / LinearAlgebra.norm(real.(vecs_hmat[1]))
    max_evec_dense_normalized = max_evec_dense / LinearAlgebra.norm(max_evec_dense)

    overlap_kr = abs(dot(evec_dense_kr_normalized, max_evec_dense_normalized))
    overlap_hmat = abs(dot(evec_hmat_normalized, max_evec_dense_normalized))

    @printf("\nEigenvector overlaps (should be ≈ 1):\n")
    @printf("  Dense Lanczos vs eigen:  %.6f\n", overlap_kr)
    @printf("  HMatrix vs eigen:        %.6f\n", overlap_hmat)

    # Performance summary
    println("\n" * "="^60)
    println("Performance Summary:")
    println("="^60)
    @printf("Matrix construction speedup: %.2fx\n", t_dense / t_hmat)
    @printf("Solver speedup (vs full eigen): %.2fx\n", t_eigen / t_lanczos_hmat)
    @printf("Solver speedup (vs dense Lanczos): %.2fx\n", t_lanczos_dense / t_lanczos_hmat)

    println("\n" * "="^60)
    println("Test completed successfully!")
    println("="^60)
end

# Run the test
main()
