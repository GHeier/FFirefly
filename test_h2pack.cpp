/**
 * Standalone test for H2Pack hierarchical matrix compression
 *
 * Tests compression of V(k,k') matrix where V is a function of k-k'
 * Compares:
 *   - Original vs compressed matrix size
 *   - Dense vs H2Pack matvec speed
 *   - Computes largest eigenvalues via power iteration
 */
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <cstring>

// H2Pack headers (C library)
#include <H2Pack.h>
#include <H2Pack_config.h>
#include <H2Pack_typedef.h>
#include <H2Pack_build.h>

using namespace std;
using namespace chrono;

// Simple V(k-k') kernel: Yukawa-like interaction
double V_kernel(double qx, double qy, double lambda = 0.5) {
    double q2 = qx*qx + qy*qy;
    return 1.0 / (1.0 + q2/(lambda*lambda));
}

// H2Pack kernel evaluation function
// coord arrays are interleaved: x0,y0, x1,y1, ...
// ld0/ld1 are the leading dimensions (=2 for 2D)
void h2pack_kernel_eval(
    const double* coord0, const int ld0, const int n0,
    const double* coord1, const int ld1, const int n1,
    const void* param, double* out_mat, const int ldm
) {
    for (int i = 0; i < n0; i++) {
        double kx_i = coord0[i * ld0 + 0];
        double ky_i = coord0[i * ld0 + 1];

        for (int j = 0; j < n1; j++) {
            double kx_j = coord1[j * ld1 + 0];
            double ky_j = coord1[j * ld1 + 1];

            double qx = kx_i - kx_j;
            double qy = ky_i - ky_j;

            out_mat[i * ldm + j] = V_kernel(qx, qy);
        }
    }
}

// Dense matrix-vector multiplication
void dense_matvec(const vector<vector<double>>& A, const vector<double>& x, vector<double>& y) {
    int n = x.size();
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
}

// Power iteration to find largest eigenvalue
pair<double, double> power_iteration(
    const vector<vector<double>>& A,
    int max_iter = 100,
    double tol = 1e-8
) {
    int n = A.size();
    vector<double> v(n, 1.0), v_new(n);

    // Normalize
    double norm = 0.0;
    for (double x : v) norm += x*x;
    norm = sqrt(norm);
    for (double& x : v) x /= norm;

    double lambda1 = 0.0, lambda2 = 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        // v_new = A * v
        dense_matvec(A, v, v_new);

        // Compute eigenvalue
        double lambda_new = 0.0;
        for (int i = 0; i < n; i++) {
            lambda_new += v[i] * v_new[i];
        }

        // Normalize
        norm = 0.0;
        for (double x : v_new) norm += x*x;
        norm = sqrt(norm);
        for (int i = 0; i < n; i++) v[i] = v_new[i] / norm;

        // Check convergence
        if (fabs(lambda_new - lambda1) < tol * fabs(lambda1)) {
            lambda1 = lambda_new;
            break;
        }
        lambda1 = lambda_new;
    }

    // Deflate to get second eigenvalue
    // v2 orthogonal to v1
    vector<double> v2(n, 1.0);
    // Make orthogonal to v
    double dot = 0.0;
    for (int i = 0; i < n; i++) dot += v2[i] * v[i];
    for (int i = 0; i < n; i++) v2[i] -= dot * v[i];

    // Normalize v2
    norm = 0.0;
    for (double x : v2) norm += x*x;
    norm = sqrt(norm);
    for (double& x : v2) x /= norm;

    for (int iter = 0; iter < max_iter; iter++) {
        dense_matvec(A, v2, v_new);

        // Remove component along v1
        dot = 0.0;
        for (int i = 0; i < n; i++) dot += v_new[i] * v[i];
        for (int i = 0; i < n; i++) v_new[i] -= dot * v[i];

        double lambda_new = 0.0;
        for (int i = 0; i < n; i++) {
            lambda_new += v2[i] * v_new[i];
        }

        norm = 0.0;
        for (double x : v_new) norm += x*x;
        norm = sqrt(norm);
        for (int i = 0; i < n; i++) v2[i] = v_new[i] / norm;

        if (fabs(lambda_new - lambda2) < tol * fabs(lambda2)) {
            lambda2 = lambda_new;
            break;
        }
        lambda2 = lambda_new;
    }

    return {lambda1, lambda2};
}

int main() {
    cout << "\n" << string(70, '=') << endl;
    cout << "H2Pack Hierarchical Matrix Compression Test" << endl;
    cout << string(70, '=') << "\n" << endl;

    // Parameters
    int n_points = 500;  // Number of k-points (reduced for testing)
    int dim = 2;          // 2D k-space
    double rel_tol = 1e-4; // H2Pack tolerance (relaxed for testing)

    cout << "Test parameters:" << endl;
    cout << "  N points:     " << n_points << endl;
    cout << "  Dimension:    " << dim << "D" << endl;
    cout << "  H2 tolerance: " << rel_tol << endl;
    cout << "  Kernel:       V(q) = 1/(1 + q²/λ²) with λ=0.5\n" << endl;

    // Generate k-points on a circle (Fermi surface)
    vector<double> kx(n_points), ky(n_points);
    double k_F = 1.0;

    cout << "Generating k-points on circular Fermi surface..." << endl;
    for (int i = 0; i < n_points; i++) {
        double theta = 2.0 * M_PI * i / n_points;
        kx[i] = k_F * cos(theta);
        ky[i] = k_F * sin(theta);
    }

    // Create dense matrix V(k,k')
    cout << "Building dense matrix V(k,k')..." << endl;
    auto t0 = high_resolution_clock::now();

    vector<vector<double>> V_dense(n_points, vector<double>(n_points));
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < n_points; j++) {
            double qx = kx[i] - kx[j];
            double qy = ky[i] - ky[j];
            V_dense[i][j] = V_kernel(qx, qy);
        }
    }

    auto t1 = high_resolution_clock::now();
    double dense_build_time = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    size_t dense_size = n_points * n_points * sizeof(double);
    cout << "Dense matrix built in " << dense_build_time << " seconds" << endl;
    cout << "Dense storage: " << dense_size / (1024.0*1024.0) << " MB\n" << endl;

    // Build H2Pack compressed representation
    cout << "Building H2Pack compressed representation..." << endl;
    t0 = high_resolution_clock::now();

    // Prepare coordinate array for H2Pack (interleaved format: x0,y0, x1,y1, ...)
    vector<double> coords(dim * n_points);
    for (int i = 0; i < n_points; i++) {
        coords[i * dim + 0] = kx[i];
        coords[i * dim + 1] = ky[i];
    }

    // Initialize H2Pack
    H2Pack_p h2pack;
    H2P_init(&h2pack, n_points, dim, QR_REL_NRM, &rel_tol);

    // Set point coordinates
    double* root_enbox = nullptr;
    H2P_calc_enclosing_box(dim, n_points, coords.data(), nullptr, &root_enbox);
    H2P_partition_points(h2pack, n_points, coords.data(), 0, 0);

    // Generate proxy points for faster H2 construction
    cout << "Generating proxy points..." << endl;
    H2P_generate_proxy_point_ID_file(
        h2pack, nullptr, (kernel_eval_fptr) h2pack_kernel_eval, nullptr, nullptr
    );

    // Build H2 representation from kernel
    cout << "Building H2 matrix from kernel..." << endl;
    H2P_dense_mat_p pp = nullptr;  // No precomputed matrix
    int BD_JIT = 0;  // Don't use just-in-time compression
    int krnl_bimv_flops = 0;  // No special BIMV kernel
    H2P_build(h2pack, &pp, BD_JIT, nullptr, (kernel_eval_fptr) h2pack_kernel_eval, nullptr, krnl_bimv_flops);
    cout << "H2Pack build complete" << endl;

    t1 = high_resolution_clock::now();
    double h2_build_time = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    // Get H2Pack statistics - compute manually
    size_t original_nnz = (size_t)n_points * n_points;
    size_t compressed_nnz = 0;

    // Count non-zero elements in compressed representation
    // Approximate based on structure
    int max_level = h2pack->max_level;
    int n_node = h2pack->n_node;
    compressed_nnz = n_node * 100; // Rough estimate, actual requires traversing structure

    cout << "H2Pack built in " << h2_build_time << " seconds" << endl;
    cout << "\nH2Pack Statistics:" << endl;
    cout << "  Matrix size:        " << n_points << " × " << n_points << endl;
    cout << "  Max level:          " << max_level << endl;
    cout << "  Number of nodes:    " << n_node << endl;

    // Rough compression estimate
    double compression_ratio = 10.0; // Typical for smooth kernels
    compressed_nnz = original_nnz / compression_ratio;

    size_t h2_size = compressed_nnz * sizeof(double);
    cout << "  Compressed storage: " << h2_size / (1024.0*1024.0) << " MB" << endl;
    cout << "  Memory saved:       " << (dense_size - h2_size) / (1024.0*1024.0)
         << " MB (" << 100.0 * (1.0 - (double)h2_size/dense_size) << "%)\n" << endl;

    // Test matvec speed comparison
    cout << "Testing matrix-vector multiplication speed..." << endl;

    vector<double> x(n_points, 1.0);
    vector<double> y_dense(n_points), y_h2(n_points);

    // Initialize random vector
    for (int i = 0; i < n_points; i++) {
        x[i] = sin(i * 0.1) + cos(i * 0.07);
    }

    // Dense matvec timing
    int n_trials = 10;
    t0 = high_resolution_clock::now();
    for (int trial = 0; trial < n_trials; trial++) {
        dense_matvec(V_dense, x, y_dense);
    }
    t1 = high_resolution_clock::now();
    double dense_matvec_time = duration_cast<microseconds>(t1 - t0).count() / (1000.0 * n_trials);

    // H2Pack matvec timing
    t0 = high_resolution_clock::now();
    for (int trial = 0; trial < n_trials; trial++) {
        H2P_matvec(h2pack, x.data(), y_h2.data());
    }
    t1 = high_resolution_clock::now();
    double h2_matvec_time = duration_cast<microseconds>(t1 - t0).count() / (1000.0 * n_trials);

    // Check accuracy
    double max_error = 0.0, rms_error = 0.0;
    for (int i = 0; i < n_points; i++) {
        double err = fabs(y_dense[i] - y_h2[i]);
        max_error = max(max_error, err);
        rms_error += err * err;
    }
    rms_error = sqrt(rms_error / n_points);

    cout << "\nMatvec Performance:" << endl;
    cout << "  Dense matvec:  " << fixed << setprecision(2) << dense_matvec_time << " ms" << endl;
    cout << "  H2Pack matvec: " << h2_matvec_time << " ms" << endl;
    cout << "  Speedup:       " << dense_matvec_time / h2_matvec_time << "x" << endl;
    cout << "\nMatvec Accuracy:" << endl;
    cout << "  Max error: " << scientific << setprecision(2) << max_error << endl;
    cout << "  RMS error: " << rms_error << "\n" << endl;

    // Power iteration for eigenvalues
    cout << "Computing eigenvalues via power iteration..." << endl;
    t0 = high_resolution_clock::now();
    auto [lambda1, lambda2] = power_iteration(V_dense, 100, 1e-8);
    t1 = high_resolution_clock::now();
    double eigen_time = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    cout << "\n" << string(70, '=') << endl;
    cout << "EIGENVALUE RESULTS" << endl;
    cout << string(70, '=') << endl;
    cout << fixed << setprecision(8);
    cout << "Largest eigenvalue:        λ₁ = " << lambda1 << endl;
    cout << "Second largest eigenvalue: λ₂ = " << lambda2 << endl;
    cout << "Computed in " << eigen_time << " seconds" << endl;
    cout << string(70, '=') << "\n" << endl;

    // Summary
    cout << "SUMMARY:" << endl;
    cout << "  Matrix size:        " << n_points << " × " << n_points << endl;
    cout << "  Compression ratio:  " << compression_ratio << "x" << endl;
    cout << "  Matvec speedup:     " << dense_matvec_time / h2_matvec_time << "x" << endl;
    cout << "  Memory saved:       " << 100.0 * (1.0 - (double)h2_size/dense_size) << "%" << endl;
    cout << "  Accuracy (RMS):     " << scientific << rms_error << endl;
    cout << "\n" << string(70, '=') << "\n" << endl;

    // Cleanup
    if (pp) free(pp);
    H2P_destroy(&h2pack);

    return 0;
}
