/**
 * @file h2matrix.cpp
 * @brief Implementation of H2Matrix class
 */

#include "h2matrix.hpp"
#include <stdexcept>
#include <cstring>
#include <iostream>
#include <algorithm>

// H2Pack C library headers
#include <H2Pack.h>
#include <H2Pack_config.h>
#include <H2Pack_typedef.h>
#include <H2Pack_build.h>

// C-linkage kernel wrapper for H2Pack
extern "C" {
    static void h2pack_default_kernel(
        const double* coord0, const int ld0, const int n0,
        const double* coord1, const int ld1, const int n1,
        const void* param, double* out_mat, const int ldm
    ) {
        // Simple Yukawa kernel: V(r) = 1/(1 + r²)
        // Coordinates are in row-major (interleaved) format: coord[i*ld + d]
        // ld0/ld1 are the leading dimensions (= pt_dim)
        int dim = ld0;

        for (int i = 0; i < n0; i++) {
            for (int j = 0; j < n1; j++) {
                double r2 = 0.0;
                for (int d = 0; d < dim; d++) {
                    double diff = coord0[i * ld0 + d] - coord1[j * ld1 + d];
                    r2 += diff * diff;
                }
                out_mat[i * ldm + j] = 1.0 / (1.0 + r2);
            }
        }
    }
}

H2Matrix::H2Matrix(
    KernelFunction kernel,
    const std::vector<std::vector<double>>& points,
    int dim,
    double rel_tol
) : kernel_(kernel), dim_(dim), rel_tol_(rel_tol) {

    n_points_ = points.size();
    if (n_points_ == 0) {
        throw std::invalid_argument("H2Matrix: points array is empty");
    }
    if (points[0].size() != static_cast<size_t>(dim)) {
        throw std::invalid_argument("H2Matrix: point dimension mismatch");
    }

    // Flatten coordinates into H2Pack row-major (interleaved) format: [x0,y0,z0, x1,y1,z1, ...]
    // Based on test_h2pack.cpp line 201-206
    coords_.resize(dim * n_points_);
    for (int i = 0; i < n_points_; i++) {
        for (int d = 0; d < dim; d++) {
            coords_[i * dim + d] = points[i][d];
        }
    }

    // Initialize H2Pack with correct parameters (matching official example_H2.c)
    H2Pack_p h2p = nullptr;
    int krnl_dim = 1;  // Scalar kernel
    H2P_init(&h2p, dim, krnl_dim, QR_REL_NRM, &rel_tol_);
    h2pack_ = (void*)h2p;

    // Calculate enclosing box and assign to H2Pack structure
    H2P_calc_enclosing_box(dim, n_points_, coords_.data(), nullptr, &h2p->root_enbox);

    // Partition points into tree structure
    // Use default max_leaf_points (0 = use H2Pack defaults: 200 for 2D, 400 for 3D)
    int max_leaf_points = 0;
    double max_leaf_size = 0.0;
    H2P_partition_points(h2p, n_points_, coords_.data(), max_leaf_points, max_leaf_size);

    // Generate proxy points for faster construction
    // Following official H2Pack example pattern from test_H2_scalar.c
    H2P_dense_mat_p *pp;
    H2P_generate_proxy_point_ID_file(
        h2p,
        nullptr,  // kernel param
        (kernel_eval_fptr) h2pack_default_kernel,
        nullptr,  // fname (no file output)
        &pp       // H2Pack will allocate and assign proxy points array
    );

    // Build H2 representation using the generated proxy points
    int BD_JIT = 0;  // Don't use just-in-time compression
    int krnl_bimv_flops = 0;  // No special BIMV kernel

    H2P_build(
        h2p,
        pp,  // Pass proxy points array
        BD_JIT,
        nullptr,  // kernel param
        (kernel_eval_fptr) h2pack_default_kernel,
        nullptr,  // krnl_bimv (no special BIMV kernel)
        krnl_bimv_flops
    );

    // Note: pp is now stored in h2pack and will be freed by H2P_destroy
}

H2Matrix::~H2Matrix() {
    if (h2pack_) {
        H2Pack_p h2p = (H2Pack_p)h2pack_;
        H2P_destroy(&h2p);
        h2pack_ = nullptr;
    }
}

std::vector<double> H2Matrix::matvec(const std::vector<double>& x) const {
    if (x.size() != static_cast<size_t>(n_points_)) {
        throw std::invalid_argument("H2Matrix::matvec: input vector size mismatch");
    }

    std::vector<double> y(n_points_);

    // H2Pack matvec operates in-place on non-const pointers
    // We need to copy x to avoid modifying it
    std::vector<double> x_copy = x;

    H2Pack_p h2p = (H2Pack_p)h2pack_;
    H2P_matvec(h2p, x_copy.data(), y.data());

    return y;
}

double H2Matrix::compression_ratio() const {
    // Estimate: original matrix has N² elements
    // Compressed form typically has O(N log N) elements
    // This is a rough estimate - exact count requires traversing H2Pack structure

    size_t original_size = static_cast<size_t>(n_points_) * n_points_;

    // Approximate compressed size based on tree depth and node count
    int max_lvl = max_level();
    int n_node = num_nodes();

    // Rough estimate: each node stores O(k²) where k is rank
    // Typical rank ~ 20-50 for smooth kernels
    size_t estimated_compressed_size = n_node * 50 * 50;

    if (estimated_compressed_size == 0) {
        return 1.0;  // No compression
    }

    return static_cast<double>(original_size) / estimated_compressed_size;
}

int H2Matrix::max_level() const {
    H2Pack_p h2p = (H2Pack_p)h2pack_;
    return h2p->max_level;
}

int H2Matrix::num_nodes() const {
    H2Pack_p h2p = (H2Pack_p)h2pack_;
    return h2p->n_node;
}

// Note: h2pack_kernel_wrapper moved to free function h2pack_default_kernel
// at top of file with extern "C" linkage
