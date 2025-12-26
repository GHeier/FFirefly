/**
 * Standalone test for V(k,k') matrix compression
 *
 * Compares dense vs compressed (low-rank approximation) representations
 * - Building a V(k,k') interaction matrix
 * - Matrix size analysis
 * - Dense vs compressed matvec timing
 * - Power iteration to find largest eigenvalues
 */
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace chrono;

// Simple low-rank compressed matrix (simulates H2Pack-style compression)
struct CompressedMatrix {
    int n;
    int rank;  // Compression rank
    vector<vector<double>> U;  // n × rank
    vector<vector<double>> V;  // rank × n

    CompressedMatrix(int n, int rank) : n(n), rank(rank) {
        U.resize(n, vector<double>(rank));
        V.resize(rank, vector<double>(n));
    }

    // Matvec: y = (U * V) * x
    void matvec(const vector<double>& x, vector<double>& y) {
        // tmp = V * x  (rank × n) * (n × 1) = (rank × 1)
        vector<double> tmp(rank, 0.0);
        #pragma omp parallel for
        for (int i = 0; i < rank; i++) {
            for (int j = 0; j < n; j++) {
                tmp[i] += V[i][j] * x[j];
            }
        }

        // y = U * tmp  (n × rank) * (rank × 1) = (n × 1)
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            y[i] = 0.0;
            for (int j = 0; j < rank; j++) {
                y[i] += U[i][j] * tmp[j];
            }
        }
    }

    size_t storage_bytes() const {
        return (n * rank + rank * n) * sizeof(double);
    }

    double compression_ratio(size_t dense_size) const {
        return (double)dense_size / storage_bytes();
    }
};

// Construct low-rank approximation using randomized SVD
CompressedMatrix compress_matrix(const vector<vector<double>>& A, int target_rank) {
    int n = A.size();
    CompressedMatrix C(n, target_rank);

    // Simple rank-k approximation using power iteration
    // This simulates what H2Pack does with hierarchical low-rank blocks

    // Start with random matrix
    vector<vector<double>> Q(n, vector<double>(target_rank));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < target_rank; j++) {
            Q[i][j] = (double)rand() / RAND_MAX - 0.5;
        }
    }

    // Power iteration to find dominant subspace
    for (int iter = 0; iter < 3; iter++) {
        vector<vector<double>> Y(n, vector<double>(target_rank, 0.0));

        // Y = A * Q
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < target_rank; k++) {
                for (int j = 0; j < n; j++) {
                    Y[i][k] += A[i][j] * Q[j][k];
                }
            }
        }

        // QR decomposition (simplified - just normalize columns)
        for (int k = 0; k < target_rank; k++) {
            // Normalize column k
            double norm = 0.0;
            for (int i = 0; i < n; i++) {
                norm += Y[i][k] * Y[i][k];
            }
            norm = sqrt(norm);

            if (norm > 1e-10) {
                for (int i = 0; i < n; i++) {
                    Q[i][k] = Y[i][k] / norm;
                }
            }
        }
    }

    // U = Q
    C.U = Q;

    // V = Q^T * A
    for (int i = 0; i < target_rank; i++) {
        for (int j = 0; j < n; j++) {
            C.V[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                C.V[i][j] += Q[k][i] * A[k][j];
            }
        }
    }

    return C;
}

// V(k-k') kernel: Yukawa-type interaction with regularization
double V_kernel(double qx, double qy, double qz = 0.0, double lambda = 0.5) {
    double q2 = qx*qx + qy*qy + qz*qz;
    double q = sqrt(q2);
    // Regularized at q=0 to avoid divergence
    return exp(-q/lambda) / (q + 0.1);  // Screened Coulomb
}

// Dense matrix-vector multiplication
void dense_matvec(const vector<vector<double>>& A, const vector<double>& x, vector<double>& y) {
    int n = x.size();
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
}

// Power iteration to find largest eigenvalue (dense matrix version)
pair<double, double> power_iteration_dense(
    const vector<vector<double>>& A,
    int max_iter = 100,
    double tol = 1e-9,
    bool verbose = true
) {
    int n = A.size();
    vector<double> v(n, 1.0), v_new(n);

    // Normalize
    double norm = 0.0;
    for (double x : v) norm += x*x;
    norm = sqrt(norm);
    for (double& x : v) x /= norm;

    double lambda1 = 0.0;

    if (verbose) {
        cout << "  Iteration  |  Eigenvalue  |  Change\n";
        cout << "  -----------|--------------|------------\n";
    }

    for (int iter = 0; iter < max_iter; iter++) {
        // v_new = A * v
        dense_matvec(A, v, v_new);

        // Compute eigenvalue (Rayleigh quotient)
        double lambda_new = 0.0;
        for (int i = 0; i < n; i++) {
            lambda_new += v[i] * v_new[i];
        }

        // Normalize
        norm = 0.0;
        for (double x : v_new) norm += x*x;
        norm = sqrt(norm);
        for (int i = 0; i < n; i++) v[i] = v_new[i] / norm;

        double change = fabs(lambda_new - lambda1);
        if (verbose && (iter % 10 == 0 || change < tol * fabs(lambda1))) {
            cout << "  " << setw(9) << iter << "  |  "
                 << fixed << setprecision(8) << setw(12) << lambda_new << "  |  "
                 << scientific << setprecision(2) << change << "\n";
        }

        // Check convergence
        if (change < tol * fabs(lambda1) && iter > 5) {
            lambda1 = lambda_new;
            if (verbose) cout << "  Converged after " << iter << " iterations\n";
            break;
        }
        lambda1 = lambda_new;
    }

    // Deflate to get second eigenvalue
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

    double lambda2 = 0.0;
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

        if (fabs(lambda_new - lambda2) < tol * fabs(lambda2) && iter > 5) {
            lambda2 = lambda_new;
            break;
        }
        lambda2 = lambda_new;
    }

    return {lambda1, lambda2};
}

// Power iteration for compressed matrix
pair<double, double> power_iteration_compressed(
    CompressedMatrix& C,
    int max_iter = 100,
    double tol = 1e-9,
    bool verbose = true
) {
    int n = C.n;
    vector<double> v(n, 1.0), v_new(n);

    // Normalize
    double norm = 0.0;
    for (double x : v) norm += x*x;
    norm = sqrt(norm);
    for (double& x : v) x /= norm;

    double lambda1 = 0.0;

    if (verbose) {
        cout << "  Iteration  |  Eigenvalue  |  Change\n";
        cout << "  -----------|--------------|------------\n";
    }

    for (int iter = 0; iter < max_iter; iter++) {
        // v_new = C * v
        C.matvec(v, v_new);

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

        double change = fabs(lambda_new - lambda1);
        if (verbose && (iter % 10 == 0 || change < tol * fabs(lambda1))) {
            cout << "  " << setw(9) << iter << "  |  "
                 << fixed << setprecision(8) << setw(12) << lambda_new << "  |  "
                 << scientific << setprecision(2) << change << "\n";
        }

        if (change < tol * fabs(lambda1) && iter > 5) {
            lambda1 = lambda_new;
            if (verbose) cout << "  Converged after " << iter << " iterations\n";
            break;
        }
        lambda1 = lambda_new;
    }

    // Second eigenvalue (deflation)
    vector<double> v2(n, 1.0);
    double dot = 0.0;
    for (int i = 0; i < n; i++) dot += v2[i] * v[i];
    for (int i = 0; i < n; i++) v2[i] -= dot * v[i];

    norm = 0.0;
    for (double x : v2) norm += x*x;
    norm = sqrt(norm);
    for (double& x : v2) x /= norm;

    double lambda2 = 0.0;
    for (int iter = 0; iter < max_iter; iter++) {
        C.matvec(v2, v_new);

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

        if (fabs(lambda_new - lambda2) < tol * fabs(lambda2) && iter > 5) {
            lambda2 = lambda_new;
            break;
        }
        lambda2 = lambda_new;
    }

    return {lambda1, lambda2};
}

int main() {
    cout << "\n" << string(70, '=') << endl;
    cout << "V(k,k') Interaction Matrix Test" << endl;
    cout << string(70, '=') << "\n" << endl;

    // Parameters
    int n_points = 2000;  // Number of k-points on Fermi surface
    int dim = 2;          // 2D Fermi surface
    double lambda = 0.5;  // Screening length

    cout << "Test parameters:" << endl;
    cout << "  N points:        " << n_points << endl;
    cout << "  Dimension:       " << dim << "D" << endl;
    cout << "  Kernel:          V(q) = exp(-q/λ) / q" << endl;
    cout << "  Screening λ:     " << lambda << endl;
    cout << "\n";

    // Generate k-points on a circular Fermi surface
    vector<double> kx(n_points), ky(n_points);
    double k_F = 1.0;  // Fermi momentum

    cout << "Generating k-points on circular Fermi surface (k_F = " << k_F << ")..." << endl;
    for (int i = 0; i < n_points; i++) {
        double theta = 2.0 * M_PI * i / n_points;
        kx[i] = k_F * cos(theta);
        ky[i] = k_F * sin(theta);
    }
    cout << "Generated " << n_points << " points\n" << endl;

    // Create dense matrix V(k,k')
    cout << "Building dense interaction matrix V(k,k')..." << endl;
    auto t0 = high_resolution_clock::now();

    vector<vector<double>> V(n_points, vector<double>(n_points));
    for (int i = 0; i < n_points; i++) {
        if (i % 200 == 0) {
            cout << "\r  Progress: " << (100 * i / n_points) << "%" << flush;
        }
        for (int j = 0; j < n_points; j++) {
            double qx = kx[i] - kx[j];
            double qy = ky[i] - ky[j];
            V[i][j] = V_kernel(qx, qy, 0.0, lambda);
        }
    }
    cout << "\r  Progress: 100%" << endl;

    auto t1 = high_resolution_clock::now();
    double build_time = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    size_t dense_size = (size_t)n_points * n_points * sizeof(double);
    cout << "Dense matrix built in " << build_time << " seconds" << endl;
    cout << "Matrix size: " << n_points << " × " << n_points << endl;
    cout << "Storage: " << dense_size / (1024.0*1024.0) << " MB\n" << endl;

    // Test matvec performance
    cout << "Testing matrix-vector multiplication performance..." << endl;
    vector<double> x(n_points), y(n_points);

    // Initialize test vector
    for (int i = 0; i < n_points; i++) {
        x[i] = sin(i * 0.1) + cos(i * 0.07);
    }

    // Dense matvec timing
    int n_trials = 10;
    cout << "  Running " << n_trials << " trials..." << endl;
    t0 = high_resolution_clock::now();
    for (int trial = 0; trial < n_trials; trial++) {
        dense_matvec(V, x, y);
    }
    t1 = high_resolution_clock::now();
    double matvec_time = duration_cast<microseconds>(t1 - t0).count() / (1000.0 * n_trials);

    cout << "  Dense matvec time: " << fixed << setprecision(2) << matvec_time << " ms" << endl;
    cout << "  Performance: " << scientific << setprecision(2)
         << 2.0 * n_points * n_points / (matvec_time * 1e-3) / 1e9
         << " GFLOP/s\n" << endl;

    // Matrix statistics
    cout << "Matrix statistics:" << endl;
    double min_val = V[0][0], max_val = V[0][0], avg_val = 0.0;
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < n_points; j++) {
            min_val = min(min_val, V[i][j]);
            max_val = max(max_val, V[i][j]);
            avg_val += V[i][j];
        }
    }
    avg_val /= (n_points * n_points);

    cout << "  Min value: " << scientific << setprecision(4) << min_val << endl;
    cout << "  Max value: " << max_val << endl;
    cout << "  Avg value: " << avg_val << "\n" << endl;

    // Build compressed representation
    cout << string(70, '=') << endl;
    cout << "Building compressed matrix representation" << endl;
    cout << string(70, '=') << "\n" << endl;

    int compression_rank = 300;  // Rank for low-rank approximation
    cout << "Using low-rank approximation with rank = " << compression_rank << endl;
    cout << "(This simulates hierarchical compression like H2Pack)\n" << endl;

    t0 = high_resolution_clock::now();
    CompressedMatrix V_compressed = compress_matrix(V, compression_rank);
    t1 = high_resolution_clock::now();
    double compress_time = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    size_t compressed_size = V_compressed.storage_bytes();
    double comp_ratio = V_compressed.compression_ratio(dense_size);

    cout << "Compressed matrix built in " << compress_time << " seconds" << endl;
    cout << "Compressed storage: " << compressed_size / (1024.0*1024.0) << " MB" << endl;
    cout << "Compression ratio: " << fixed << setprecision(1) << comp_ratio << "x" << endl;
    cout << "Memory saved: " << (dense_size - compressed_size) / (1024.0*1024.0)
         << " MB (" << 100.0 * (1.0 - (double)compressed_size/dense_size) << "%)\n" << endl;

    // Test compressed matvec
    cout << "Testing compressed matrix-vector multiplication..." << endl;
    vector<double> y_compressed(n_points);

    t0 = high_resolution_clock::now();
    for (int trial = 0; trial < n_trials; trial++) {
        V_compressed.matvec(x, y_compressed);
    }
    t1 = high_resolution_clock::now();
    double compressed_matvec_time = duration_cast<microseconds>(t1 - t0).count() / (1000.0 * n_trials);

    // Check accuracy
    double max_error = 0.0, rms_error = 0.0;
    double y_norm = 0.0;
    for (int i = 0; i < n_points; i++) {
        double err = fabs(y[i] - y_compressed[i]);
        max_error = max(max_error, err);
        rms_error += err * err;
        y_norm += y[i] * y[i];
    }
    rms_error = sqrt(rms_error / n_points);
    y_norm = sqrt(y_norm);
    double relative_error = rms_error / y_norm;

    cout << "  Compressed matvec time: " << fixed << setprecision(2)
         << compressed_matvec_time << " ms" << endl;
    cout << "  Speedup vs dense: " << matvec_time / compressed_matvec_time << "x" << endl;
    cout << "  Accuracy (RMS error): " << scientific << setprecision(2)
         << relative_error << "\n" << endl;

    // Power iteration for eigenvalues - DENSE MATRIX
    cout << string(70, '=') << endl;
    cout << "Computing eigenvalues: DENSE MATRIX" << endl;
    cout << string(70, '=') << "\n" << endl;

    t0 = high_resolution_clock::now();
    auto [lambda1_dense, lambda2_dense] = power_iteration_dense(V, 100, 1e-9, true);
    t1 = high_resolution_clock::now();
    double eigen_time_dense = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    cout << "\nDense matrix eigenvalues computed in " << eigen_time_dense << " seconds\n" << endl;

    // Power iteration for eigenvalues - COMPRESSED MATRIX
    cout << string(70, '=') << endl;
    cout << "Computing eigenvalues: COMPRESSED MATRIX" << endl;
    cout << string(70, '=') << "\n" << endl;

    t0 = high_resolution_clock::now();
    auto [lambda1_compressed, lambda2_compressed] = power_iteration_compressed(V_compressed, 100, 1e-9, true);
    t1 = high_resolution_clock::now();
    double eigen_time_compressed = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    cout << "\nCompressed matrix eigenvalues computed in " << eigen_time_compressed << " seconds\n" << endl;

    // Final comparison
    cout << "\n" << string(70, '=') << endl;
    cout << "FINAL COMPARISON: DENSE vs COMPRESSED" << endl;
    cout << string(70, '=') << "\n" << endl;

    cout << "DENSE MATRIX:\n";
    cout << "  Matrix size:          " << n_points << " × " << n_points << endl;
    cout << "  Storage:              " << fixed << setprecision(2)
         << dense_size / (1024.0*1024.0) << " MB" << endl;
    cout << "  Matvec time:          " << setprecision(3) << matvec_time << " ms" << endl;
    cout << "  Largest eigenvalue:   " << setprecision(10) << lambda1_dense << endl;
    cout << "  2nd eigenvalue:       " << lambda2_dense << "\n" << endl;

    cout << "COMPRESSED MATRIX:\n";
    cout << "  Matrix size:          " << n_points << " × " << n_points
         << " (rank " << compression_rank << ")" << endl;
    cout << "  Storage:              " << setprecision(2)
         << compressed_size / (1024.0*1024.0) << " MB" << endl;
    cout << "  Matvec time:          " << setprecision(3) << compressed_matvec_time << " ms" << endl;
    cout << "  Largest eigenvalue:   " << setprecision(10) << lambda1_compressed << endl;
    cout << "  2nd eigenvalue:       " << lambda2_compressed << "\n" << endl;

    cout << "IMPROVEMENT:\n";
    cout << "  Compression ratio:    " << setprecision(1) << comp_ratio << "x" << endl;
    cout << "  Matvec speedup:       " << setprecision(2)
         << matvec_time / compressed_matvec_time << "x" << endl;
    cout << "  Memory saved:         " << setprecision(2)
         << (dense_size - compressed_size) / (1024.0*1024.0) << " MB ("
         << setprecision(1) << 100.0 * (1.0 - (double)compressed_size/dense_size) << "%)" << endl;
    cout << "  Eigenvalue error:     " << scientific << setprecision(2)
         << fabs(lambda1_dense - lambda1_compressed) / fabs(lambda1_dense) * 100.0 << "%" << endl;
    cout << "\n" << string(70, '=') << "\n" << endl;

    return 0;
}
