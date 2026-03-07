module DenseMatrix

export create_full_matrix

using LinearAlgebra
using Printf

function create_full_matrix(V, dos_weights, kpoints, w_points, kernel_func)
    n_k = length(kpoints)
    n_w = length(w_points)
    n_total = n_k * n_w

    @printf("Building dense matrix for %d k-points × %d frequencies = %d total points\n", n_k, n_w, n_total)

    # Allocate the full matrix
    M = zeros(Float64, n_total, n_total)

    # Fill the matrix using the kernel function
    # Index ordering: k varies faster (idx = i_k + (i_w - 1) * n_k)
    for i_w in 1:n_w
        w1 = w_points[i_w]
        for i_k in 1:n_k
            k1 = kpoints[i_k]
            idx1 = i_k + (i_w - 1) * n_k

            for j_w in 1:n_w
                w2 = w_points[j_w]
                for j_k in 1:n_k
                    k2 = kpoints[j_k]
                    idx2 = j_k + (j_w - 1) * n_k

                    M[idx1, idx2] = kernel_func(k1, k2, w1, w2, i_k, j_k, i_w, j_w)
                end
            end
        end
    end

    @printf("Dense matrix built: %d × %d\n", n_total, n_total)

    return M
end

end # module
