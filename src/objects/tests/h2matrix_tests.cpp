/**
 * @file h2matrix_tests.cpp
 * @brief Implementation of H2Matrix tests
 */

#include "../../config/load/c_config.h"
#include "../h2matrix.hpp"
#include "h2matrix_tests.hpp"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

const double H2MATRIX_TEST_TOL = 1e-3;  // Relaxed tolerance for H2Pack compression

/**
 * @brief Test H2Matrix creation from kernel function
 *
 * Creates an H2Matrix from a simple Yukawa-like kernel and verifies
 * that the matrix is built correctly with reasonable compression.
 */
bool test_h2matrix_creation() {
    printf("  Starting test_h2matrix_creation...\n");
    // Define simple kernel: V(r) = 1/(1 + r²)
    auto kernel = [](const double* x1, const double* x2, int dim) -> double {
        double r2 = 0.0;
        for (int i = 0; i < dim; i++) {
            double d = x1[i] - x2[i];
            r2 += d * d;
        }
        return 1.0 / (1.0 + r2);
    };

    // Create points on a sphere (3D) - H2Pack has issues with 2D
    int n_points = 500;
    int dim = 3;
    vector<vector<double>> points(n_points, vector<double>(dim));

    double radius = 1.0;
    for (int i = 0; i < n_points; i++) {
        double theta = 2.0 * M_PI * i / n_points;
        double phi = M_PI * (i % 20) / 20.0;
        points[i][0] = radius * sin(phi) * cos(theta);
        points[i][1] = radius * sin(phi) * sin(theta);
        points[i][2] = radius * cos(phi);
    }

    // Build H2Matrix
    H2Matrix mat(kernel, points, dim, 1e-3);

    // Check basic properties
    if (mat.size() != n_points) {
        printf("H2Matrix creation test failed: wrong size (%d vs %d)\n",
               mat.size(), n_points);
        return false;
    }

    if (mat.dimension() != dim) {
        printf("H2Matrix creation test failed: wrong dimension (%d vs %d)\n",
               mat.dimension(), dim);
        return false;
    }

    // Check that tree was built
    if (mat.max_level() <= 0) {
        printf("H2Matrix creation test failed: tree not built (max_level = %d)\n",
               mat.max_level());
        return false;
    }

    if (mat.num_nodes() <= 0) {
        printf("H2Matrix creation test failed: no tree nodes (num_nodes = %d)\n",
               mat.num_nodes());
        return false;
    }

    // Check that compression occurred
    double comp_ratio = mat.compression_ratio();
    if (comp_ratio < 1.5) {
        printf("H2Matrix creation test failed: insufficient compression (ratio = %.2f)\n",
               comp_ratio);
        return false;
    }

    return true;
}

/**
 * @brief Test H2Matrix matvec operation
 *
 * Compares H2Matrix matvec with dense matrix-vector multiplication
 * to verify correctness within tolerance.
 */
bool test_h2matrix_matvec() {
    printf("  Starting test_h2matrix_matvec...\n");

    // Define kernel: V(r) = exp(-r)
    auto kernel = [](const double* x1, const double* x2, int dim) -> double {
        double r2 = 0.0;
        for (int i = 0; i < dim; i++) {
            double d = x1[i] - x2[i];
            r2 += d * d;
        }
        double r = sqrt(r2);
        return exp(-r);
    };

    // Create points on a 3D grid - H2Pack has issues with 2D
    int n_per_dim = 5;  // 5×5×5 = 125 points
    int dim = 3;
    int n_points = n_per_dim * n_per_dim * n_per_dim;
    vector<vector<double>> points(n_points, vector<double>(dim));

    int idx = 0;
    for (int i = 0; i < n_per_dim; i++) {
        for (int j = 0; j < n_per_dim; j++) {
            for (int k = 0; k < n_per_dim; k++) {
                points[idx][0] = i * 0.2;
                points[idx][1] = j * 0.2;
                points[idx][2] = k * 0.2;
                idx++;
            }
        }
    }

    printf("    Building H2Matrix with %d points...\n", n_points);
    // Build H2Matrix
    H2Matrix mat(kernel, points, dim, 1e-3);

    // Create test vector: x[i] = sin(i * 0.1)
    vector<double> x(n_points);
    for (int i = 0; i < n_points; i++) {
        x[i] = sin(i * 0.1);
    }

    printf("    Computing H2 matvec...\n");
    // Compute H2 matvec
    vector<double> y_h2 = mat.matvec(x);

    printf("    Computing dense matvec for comparison...\n");
    // Compute dense matvec for comparison
    vector<double> y_dense(n_points, 0.0);
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < n_points; j++) {
            double val = kernel(points[i].data(), points[j].data(), dim);
            y_dense[i] += val * x[j];
        }
    }

    printf("    Comparing results...\n");
    // Compare results
    double max_error = 0.0;
    double rms_error = 0.0;
    for (int i = 0; i < n_points; i++) {
        double err = fabs(y_h2[i] - y_dense[i]);
        max_error = max(max_error, err);
        rms_error += err * err;
    }
    rms_error = sqrt(rms_error / n_points);

    printf("    Max error: %.2e, RMS error: %.2e (tolerance: %.2e)\n",
           max_error, rms_error, H2MATRIX_TEST_TOL);

    // Check accuracy
    if (max_error > H2MATRIX_TEST_TOL) {
        printf("H2Matrix matvec test failed: max error = %.2e (tolerance = %.2e)\n",
               max_error, H2MATRIX_TEST_TOL);
        return false;
    }

    if (rms_error > H2MATRIX_TEST_TOL) {
        printf("H2Matrix matvec test failed: RMS error = %.2e (tolerance = %.2e)\n",
               rms_error, H2MATRIX_TEST_TOL);
        return false;
    }

    printf("  test_h2matrix_matvec completed successfully\n");
    return true;
}

/**
 * @brief Run all H2Matrix tests
 */
bool h2matrix_tests() {
    printf("Starting H2Matrix test suite...\n");
    int num_tests = 2;
    bool all_tests[2];

    all_tests[0] = test_h2matrix_creation();
    printf("  Finished creation test\n");

    all_tests[1] = test_h2matrix_matvec();
    printf("  Finished matvec test\n");

    return print_test_results(all_tests, num_tests, "H2Matrix tests");
}
