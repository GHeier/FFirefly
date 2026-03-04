module HMatrixHelper

export build_hmatrix

using HMatrices
using LinearAlgebra
using Printf
using StaticArrays

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

function build_spoints(kpoints, frequencies)
    n_k = length(kpoints)
    n_w = length(frequencies)
    n_total = n_k * n_w

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
    return X
end

# Build HMatrix with frequency dependence
function build_hmatrix(kpoints, frequencies, kernel_func; atol=1e-6, rank=20)
    n_k = length(kpoints)
    n_w = length(frequencies)
    n_total = n_k * n_w
    @printf("Building HMatrix for %d k-points × %d frequencies = %d total points\n", n_k, n_w, n_total)
    # Create frequency-dependent kernel matrix wrapper
    K = FreqKernelMatrixWrapper(kernel_func, kpoints, frequencies, n_k, n_w)

    X = build_spoints(kpoints, frequencies)
    # Build the hierarchical matrix
    println("Building HMatrix...")
    splitter = HMatrices.CardinalitySplitter(; nmax=25)
    Xclt = HMatrices.ClusterTree(X, splitter)
    Yclt = Xclt  # Same cluster tree for row/column spaces

    # Create compression method
    comp = HMatrices.PartialACA(; atol=atol, rank=rank)
    #comp = HMatrices.TSVD(; atol=atol, rank=rank)

    # Assemble HMatrix from matrix and cluster trees
    H = assemble_hmatrix(K, Xclt, Yclt; comp=comp)

    return H
end

end # module
