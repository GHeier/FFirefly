/**
 * H2Pack wrapper implementation
 */
#include "h2pack_wrapper.hpp"
#include "../hamiltonian/band_structure.hpp"
#include "../objects/CMField/vertex.hpp"
#include "cfg.hpp"
#include <cmath>
#include <iostream>
#include <chrono>

// H2Pack headers
#include <H2Pack.h>
#include <H2Pack_config.h>

using namespace std;

// Global parameter structure for kernel evaluation
struct BCS_Kernel_Params {
    Vertex* V_func;
    float renorm;
    vector<Vec>* FS_points;
};

H2PackMatrix::H2PackMatrix(int dim, double rel_tol)
    : h2pack(nullptr), dim(dim), rel_tol(rel_tol),
      original_nnz(0), compressed_nnz(0), compression_ratio(0.0),
      max_rank(0), build_time(0.0), matvec_time(0.0) {
}

H2PackMatrix::~H2PackMatrix() {
    if (h2pack != nullptr) {
        H2P_destroy(&h2pack);
    }
}

void H2PackMatrix::build_from_kernel(const vector<Vec>& FS, float renorm) {
    auto start = chrono::high_resolution_clock::now();

    n_points = FS.size();
    points = FS;
    original_nnz = (long)n_points * n_points;

    cout << "\n================================" << endl;
    cout << "H2Pack Matrix Construction" << endl;
    cout << "================================" << endl;
    cout << "Number of points: " << n_points << endl;
    cout << "Dimension: " << dim << endl;
    cout << "Relative tolerance: " << rel_tol << endl;

    // Prepare coordinate array for H2Pack
    vector<double> coords(n_points * dim);
    for (int i = 0; i < n_points; i++) {
        coords[i * dim + 0] = points[i].x;
        coords[i * dim + 1] = points[i].y;
        if (dim == 3) coords[i * dim + 2] = points[i].z;
    }

    // Setup kernel parameters
    Vertex V_func;
    BCS_Kernel_Params params;
    params.V_func = &V_func;
    params.renorm = renorm;
    params.FS_points = &points;

    // ===== H2Pack construction =====

    // Initialize H2Pack
    H2P_init(&h2pack, n_points, dim, QR_REL_NRM, &rel_tol);

    // Set point coordinates
    DTYPE *root_enbox = nullptr;
    H2P_calc_enclosing_box(dim, n_points, coords.data(), nullptr, &root_enbox);

    // Partition the point set
    H2P_partition_points(h2pack, n_points, coords.data(), 0, 0);

    // Build the H2 representation using proxy point method
    H2P_generate_proxy_point_ID_file(
        h2pack, (void*)&params, bcs_pairing_kernel, nullptr, nullptr
    );

    H2P_build(
        h2pack,
        nullptr,
        0,
        (void*)&params,
        bcs_pairing_kernel,
        nullptr,
        0  // use_proxy_file = 0 (generate on the fly)
    );

    // ===== End H2Pack construction =====

    auto end = chrono::high_resolution_clock::now();
    build_time = chrono::duration<double>(end - start).count();

    compute_stats();
    print_stats();
}

void H2PackMatrix::build_from_matrix(const Matrix& P, const vector<Vec>& FS) {
    auto start = chrono::high_resolution_clock::now();

    n_points = FS.size();
    points = FS;
    original_nnz = (long)n_points * n_points;

    cout << "\n================================" << endl;
    cout << "H2Pack Matrix from Dense Matrix" << endl;
    cout << "================================" << endl;

    // Prepare coordinates
    vector<double> coords(n_points * dim);
    for (int i = 0; i < n_points; i++) {
        coords[i * dim + 0] = points[i].x;
        coords[i * dim + 1] = points[i].y;
        if (dim == 3) coords[i * dim + 2] = points[i].z;
    }

    // ===== H2Pack from dense matrix =====

    // Initialize H2Pack
    H2P_init(&h2pack, n_points, dim, QR_REL_NRM, &rel_tol);

    // Set coordinates and partition
    DTYPE *root_enbox = nullptr;
    H2P_calc_enclosing_box(dim, n_points, coords.data(), nullptr, &root_enbox);
    H2P_partition_points(h2pack, n_points, coords.data(), 0, 0);

    // Convert Matrix to dense array for H2Pack
    vector<double> P_dense(n_points * n_points);
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < n_points; j++) {
            P_dense[i * n_points + j] = const_cast<Matrix&>(P)(i, j);
        }
    }

    // Build H2 from dense matrix using sample points
    H2P_build_with_sample_point(
        h2pack,
        nullptr,
        0,
        (void*)P_dense.data(),
        nullptr,
        nullptr,
        0
    );

    // ===== End H2Pack from dense matrix =====

    auto end = chrono::high_resolution_clock::now();
    build_time = chrono::duration<double>(end - start).count();

    compute_stats();
    print_stats();
}

void H2PackMatrix::matvec(const vector<float>& x, vector<float>& y) {
    auto start = chrono::high_resolution_clock::now();

    if (x.size() != n_points || y.size() != n_points) {
        cerr << "Error: vector size mismatch in H2Pack matvec" << endl;
        return;
    }

    // ===== H2Pack matvec =====

    // Convert float to double for H2Pack
    vector<double> x_d(n_points), y_d(n_points);
    for (int i = 0; i < n_points; i++) x_d[i] = x[i];

    // Perform H2 matrix-vector multiplication
    H2P_matvec(h2pack, x_d.data(), y_d.data());

    // Convert back to float
    for (int i = 0; i < n_points; i++) y[i] = y_d[i];

    // ===== End H2Pack matvec =====

    auto end = chrono::high_resolution_clock::now();
    matvec_time = chrono::duration<double>(end - start).count();
}

void H2PackMatrix::compute_stats() {
    // ===== Compute statistics =====

    // Get storage info from H2Pack structure
    // mat_size[0] = U matrices, [1] = B matrices, [2] = D matrices
    compressed_nnz = h2pack->mat_size[0] + h2pack->mat_size[1] + h2pack->mat_size[2];

    // Get max rank from QR_stop_rank (max rank used in compression)
    max_rank = h2pack->QR_stop_rank;

    compression_ratio = (double)original_nnz / (double)compressed_nnz;

    // ===== End statistics =====
}

void H2PackMatrix::print_stats() const {
    cout << "\n================================" << endl;
    cout << "H2Pack Compression Statistics" << endl;
    cout << "================================" << endl;
    cout << "Matrix size:           " << n_points << " x " << n_points << endl;
    cout << "Original storage:      " << original_nnz << " elements" << endl;
    cout << "Compressed storage:    " << compressed_nnz << " elements" << endl;
    cout << "Compression ratio:     " << compression_ratio << "x" << endl;
    cout << "Maximum block rank:    " << max_rank << endl;
    cout << "Build time:            " << build_time << " seconds" << endl;
    cout << "================================" << endl;

    // Memory savings
    double orig_mem_mb = (original_nnz * sizeof(float)) / (1024.0 * 1024.0);
    double comp_mem_mb = (compressed_nnz * sizeof(double)) / (1024.0 * 1024.0);
    cout << "Memory: " << orig_mem_mb << " MB → " << comp_mem_mb << " MB" << endl;
    cout << "Savings: " << (orig_mem_mb - comp_mem_mb) << " MB ("
         << 100.0 * (1.0 - comp_mem_mb/orig_mem_mb) << "%)" << endl;
    cout << "================================\n" << endl;
}

// BCS pairing kernel function for H2Pack
void bcs_pairing_kernel(
    const double* coord0, const int ld0, const int n0,
    const double* coord1, const int ld1, const int n1,
    const void* krnl_param, double* mat, const int ldm
) {
    BCS_Kernel_Params* p = (BCS_Kernel_Params*)krnl_param;
    Vertex& V_func = *(p->V_func);
    vector<Vec>& FS = *(p->FS_points);

    // Get dimension from first point
    int dim_pt = (FS.size() > 0 && FS[0].z != 0.0) ? 3 : 2;

    // Evaluate kernel for all pairs
    // coord0 and coord1 are stored in column-major format (pt_dim x n_points)
    for (int i = 0; i < n0; i++) {
        Vec k1;
        k1.x = coord0[0 * ld0 + i];  // First dimension
        k1.y = coord0[1 * ld0 + i];  // Second dimension
        if (dim_pt == 3) k1.z = coord0[2 * ld0 + i];

        // Find corresponding Fermi surface point for form factors
        int idx1 = i;  // Assume ordering matches
        float f1 = pow(FS[idx1].area / vp(FS[idx1].n, FS[idx1]), 0.5);

        for (int j = 0; j < n1; j++) {
            Vec k2;
            k2.x = coord1[0 * ld1 + j];
            k2.y = coord1[1 * ld1 + j];
            if (dim_pt == 3) k2.z = coord1[2 * ld1 + j];

            int idx2 = j;
            float f2 = pow(FS[idx2].area / vp(FS[idx2].n, FS[idx2]), 0.5);

            // BCS pairing kernel: P(k1, k2) = -f1 * f2 * [V(k1-k2) + V(k1+k2)]/2
            double V_val = (V_func(k1 - k2, 0).real() + V_func(k1 + k2, 0).real()) / 2.0;
            mat[i * ldm + j] = -f1 * f2 * V_val;
        }
    }
}
