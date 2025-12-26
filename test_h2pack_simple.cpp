/**
 * Simple standalone H2Pack test - based on H2Pack examples
 */
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <cstring>

#include <H2Pack.h>

using namespace std;
using namespace chrono;

// Simple Yukawa kernel: V(r) = exp(-r) / r
void kernel_eval(
    const double* coord0, const int ld0, const int n0,
    const double* coord1, const int ld1, const int n1,
    const void* param, double* out_mat, const int ldm
) {
    int dim = ld0;
    for (int i = 0; i < n0; i++) {
        for (int j = 0; j < n1; j++) {
            double dx = coord0[i * dim + 0] - coord1[j * dim + 0];
            double dy = coord0[i * dim + 1] - coord1[j * dim + 1];
            double r = sqrt(dx*dx + dy*dy) + 1e-10;
            out_mat[i * ldm + j] = exp(-r) / r;
        }
    }
}

int main() {
    cout << "\n====================================================================\n";
    cout << "Simple H2Pack Test - V(k,k') Compression\n";
    cout << "====================================================================\n\n";

    // Parameters
    int n_points = 1000;
    int dim = 2;
    double rel_tol = 1e-4;

    cout << "Parameters:\n";
    cout << "  N points:  " << n_points << "\n";
    cout << "  Dimension: " << dim << "D\n";
    cout << "  Tolerance: " << rel_tol << "\n\n";

    // Generate circular k-points
    vector<double> coords(n_points * dim);
    double k_F = 1.0;
    for (int i = 0; i < n_points; i++) {
        double theta = 2.0 * M_PI * i / n_points;
        coords[i * dim + 0] = k_F * cos(theta);
        coords[i * dim + 1] = k_F * sin(theta);
    }

    // Build dense matrix for reference
    cout << "Building dense matrix...\n";
    auto t0 = high_resolution_clock::now();
    vector<vector<double>> V_dense(n_points, vector<double>(n_points));
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < n_points; j++) {
            double dx = coords[i*dim+0] - coords[j*dim+0];
            double dy = coords[i*dim+1] - coords[j*dim+1];
            double r = sqrt(dx*dx + dy*dy) + 1e-10;
            V_dense[i][j] = exp(-r) / r;
        }
    }
    auto t1 = high_resolution_clock::now();
    double dense_time = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;
    size_t dense_size = (size_t)n_points * n_points * sizeof(double);
    cout << "  Dense matrix built in " << dense_time << " s\n";
    cout << "  Storage: " << dense_size / (1024.0*1024.0) << " MB\n\n";

    // Initialize H2Pack
    cout << "Initializing H2Pack...\n";
    H2Pack_p h2pack;
    H2P_init(&h2pack, n_points, dim, QR_REL_NRM, &rel_tol);

    double* root_enbox = nullptr;
    H2P_calc_enclosing_box(dim, n_points, coords.data(), nullptr, &root_enbox);
    H2P_partition_points(h2pack, n_points, coords.data(), 0, 0.0);
    cout << "  Partitioning complete\n";
    cout << "  Max level: " << h2pack->max_level << "\n";
    cout << "  Num nodes: " << h2pack->n_node << "\n\n";

    // Select sample points
    cout << "Selecting sample points...\n";
    t0 = high_resolution_clock::now();
    H2P_dense_mat_p* sample_pt = nullptr;
    double tau = 0.7;
    H2P_select_sample_point(h2pack, nullptr, kernel_eval, tau, &sample_pt);
    t1 = high_resolution_clock::now();
    cout << "  Sample point selection: " << duration_cast<milliseconds>(t1-t0).count()/1000.0 << " s\n\n";

    // Build H2 representation
    cout << "Building H2 representation...\n";
    t0 = high_resolution_clock::now();
    H2P_build_with_sample_point(h2pack, sample_pt, 0, nullptr, kernel_eval, nullptr, 0);
    t1 = high_resolution_clock::now();
    double h2_build_time = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;
    cout << "  H2 build time: " << h2_build_time << " s\n\n";

    // Print H2Pack statistics
    cout << "====================================================================\n";
    cout << "H2Pack Compression Statistics\n";
    cout << "====================================================================\n";
    H2P_print_statistic(h2pack);
    cout << "\n";

    // Test matvec
    cout << "Testing matrix-vector multiplication...\n";
    vector<double> x(n_points), y_dense(n_points), y_h2(n_points);
    for (int i = 0; i < n_points; i++) {
        x[i] = sin(i * 0.1) + cos(i * 0.07);
    }

    // Dense matvec
    t0 = high_resolution_clock::now();
    for (int trial = 0; trial < 5; trial++) {
        for (int i = 0; i < n_points; i++) {
            y_dense[i] = 0.0;
            for (int j = 0; j < n_points; j++) {
                y_dense[i] += V_dense[i][j] * x[j];
            }
        }
    }
    t1 = high_resolution_clock::now();
    double dense_matvec_time = duration_cast<microseconds>(t1 - t0).count() / (1000.0 * 5);

    // H2 matvec
    t0 = high_resolution_clock::now();
    for (int trial = 0; trial < 5; trial++) {
        H2P_matvec(h2pack, x.data(), y_h2.data());
    }
    t1 = high_resolution_clock::now();
    double h2_matvec_time = duration_cast<microseconds>(t1 - t0).count() / (1000.0 * 5);

    // Check accuracy
    double max_error = 0.0, rms_error = 0.0;
    for (int i = 0; i < n_points; i++) {
        double err = fabs(y_dense[i] - y_h2[i]);
        max_error = max(max_error, err);
        rms_error += err * err;
    }
    rms_error = sqrt(rms_error / n_points);

    cout << "\nMatvec Performance:\n";
    cout << "  Dense: " << fixed << setprecision(2) << dense_matvec_time << " ms\n";
    cout << "  H2:    " << h2_matvec_time << " ms\n";
    cout << "  Speedup: " << dense_matvec_time / h2_matvec_time << "x\n";
    cout << "\nAccuracy:\n";
    cout << "  Max error: " << scientific << setprecision(2) << max_error << "\n";
    cout << "  RMS error: " << rms_error << "\n\n";

    // Power iteration for eigenvalues
    cout << "Computing eigenvalues via power iteration...\n";
    t0 = high_resolution_clock::now();

    vector<double> v(n_points, 1.0), v_new(n_points);
    double norm = 0.0;
    for (double val : v) norm += val*val;
    norm = sqrt(norm);
    for (double& val : v) val /= norm;

    double lambda1 = 0.0;
    for (int iter = 0; iter < 50; iter++) {
        for (int i = 0; i < n_points; i++) {
            v_new[i] = 0.0;
            for (int j = 0; j < n_points; j++) {
                v_new[i] += V_dense[i][j] * v[j];
            }
        }

        double lambda_new = 0.0;
        for (int i = 0; i < n_points; i++) {
            lambda_new += v[i] * v_new[i];
        }

        norm = 0.0;
        for (double val : v_new) norm += val*val;
        norm = sqrt(norm);
        for (int i = 0; i < n_points; i++) v[i] = v_new[i] / norm;

        if (fabs(lambda_new - lambda1) < 1e-8 * fabs(lambda1)) {
            lambda1 = lambda_new;
            break;
        }
        lambda1 = lambda_new;
    }

    // Second eigenvalue
    vector<double> v2(n_points, 1.0);
    double dot = 0.0;
    for (int i = 0; i < n_points; i++) dot += v2[i] * v[i];
    for (int i = 0; i < n_points; i++) v2[i] -= dot * v[i];
    norm = 0.0;
    for (double val : v2) norm += val*val;
    norm = sqrt(norm);
    for (double& val : v2) val /= norm;

    double lambda2 = 0.0;
    for (int iter = 0; iter < 50; iter++) {
        for (int i = 0; i < n_points; i++) {
            v_new[i] = 0.0;
            for (int j = 0; j < n_points; j++) {
                v_new[i] += V_dense[i][j] * v2[j];
            }
        }

        dot = 0.0;
        for (int i = 0; i < n_points; i++) dot += v_new[i] * v[i];
        for (int i = 0; i < n_points; i++) v_new[i] -= dot * v[i];

        double lambda_new = 0.0;
        for (int i = 0; i < n_points; i++) {
            lambda_new += v2[i] * v_new[i];
        }

        norm = 0.0;
        for (double val : v_new) norm += val*val;
        norm = sqrt(norm);
        for (int i = 0; i < n_points; i++) v2[i] = v_new[i] / norm;

        if (fabs(lambda_new - lambda2) < 1e-8 * fabs(lambda2)) {
            lambda2 = lambda_new;
            break;
        }
        lambda2 = lambda_new;
    }

    t1 = high_resolution_clock::now();
    double eigen_time = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    cout << "\n====================================================================\n";
    cout << "EIGENVALUE RESULTS\n";
    cout << "====================================================================\n";
    cout << fixed << setprecision(8);
    cout << "Largest eigenvalue:        λ₁ = " << lambda1 << "\n";
    cout << "Second largest eigenvalue: λ₂ = " << lambda2 << "\n";
    cout << "Computed in " << eigen_time << " seconds\n";
    cout << "====================================================================\n\n";

    // Summary
    cout << "SUMMARY:\n";
    cout << "  Matrix size:       " << n_points << " × " << n_points << "\n";
    cout << "  Dense storage:     " << dense_size / (1024.0*1024.0) << " MB\n";
    cout << "  Matvec speedup:    " << dense_matvec_time / h2_matvec_time << "x\n";
    cout << "  Accuracy (RMS):    " << scientific << rms_error << "\n";
    cout << "  Top eigenvalue:    " << fixed << lambda1 << "\n";
    cout << "====================================================================\n\n";

    // Cleanup
    H2P_destroy(&h2pack);
    if (root_enbox) free(root_enbox);

    return 0;
}
