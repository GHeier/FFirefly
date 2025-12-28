#include "hmatrix_tests.hpp"
#include "../hmatrix.hpp"
#include "all.hpp"
#include "../../config/load/c_config.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>

using namespace std;

/**
 * Test 1: H-matrix creation and basic properties
 *
 * Creates an H-matrix from points on a circle with an exponential kernel
 * and verifies that the matrix is created successfully with reasonable
 * compression ratios.
 */
bool test_hmatrix_creation() {
    printf("  Test 1: H-matrix creation... ");

    try {
        // Create points on a circle
        uint n_points = 1000;
        vector<vector<double>> points(n_points);

        for (uint i = 0; i < n_points; i++) {
            double angle = 2.0 * M_PI * i / n_points;
            points[i] = {cos(angle), sin(angle)};
        }

        // Define exponential kernel: k(x,y) = exp(-|x-y|^2)
        auto kernel = [](const double* x, const double* y) -> double {
            double dx = x[0] - y[0];
            double dy = x[1] - y[1];
            double dist2 = dx * dx + dy * dy;
            return exp(-dist2);
        };

        // Build H-matrix
        HMatrix hm(points, kernel, 4, 32, 1e-6, 2.0);

        // Verify basic properties
        if (hm.size() != n_points) {
            printf("FAILED: Wrong matrix size\n");
            return false;
        }

        // Check that compression is achieved
        double compression = hm.compression_ratio();
        if (compression < 2.0) {
            printf("FAILED: Insufficient compression (%.2fx)\n", compression);
            return false;
        }

        // Print statistics (optional, comment out for quiet tests)
        // hm.print_stats();

        printf("PASSED (compression: %.2fx)\n", compression);
        return true;

    } catch (const exception& e) {
        printf("FAILED: Exception thrown: %s\n", e.what());
        return false;
    }
}

/**
 * Test 2: Matrix-vector multiplication accuracy
 *
 * Compares H-matrix matvec with dense matrix matvec to verify correctness.
 * Uses a smaller matrix size for feasibility of dense computation.
 */
bool test_hmatrix_matvec() {
    printf("  Test 2: Matrix-vector multiplication... ");

    try {
        // Create smaller point set for dense comparison
        uint n_points = 200;
        vector<vector<double>> points(n_points);

        for (uint i = 0; i < n_points; i++) {
            double angle = 2.0 * M_PI * i / n_points;
            points[i] = {cos(angle), sin(angle)};
        }

        // Define exponential kernel
        auto kernel = [](const double* x, const double* y) -> double {
            double dx = x[0] - y[0];
            double dy = x[1] - y[1];
            double dist2 = dx * dx + dy * dy;
            return exp(-dist2);
        };

        // Build H-matrix
        HMatrix hm(points, kernel, 4, 16, 1e-6, 2.0);

        // Create test vector
        vector<double> x(n_points);
        for (uint i = 0; i < n_points; i++) {
            x[i] = sin(2.0 * M_PI * i / n_points);
        }

        // H-matrix matvec
        vector<double> y_hmat = hm.matvec(x);

        // Dense matrix matvec (ground truth)
        vector<vector<double>> dense = hm.to_dense();
        vector<double> y_dense(n_points, 0.0);

        for (uint i = 0; i < n_points; i++) {
            for (uint j = 0; j < n_points; j++) {
                y_dense[i] += dense[i][j] * x[j];
            }
        }

        // Compute relative error
        double error_norm = 0.0;
        double y_norm = 0.0;

        for (uint i = 0; i < n_points; i++) {
            double diff = y_hmat[i] - y_dense[i];
            error_norm += diff * diff;
            y_norm += y_dense[i] * y_dense[i];
        }

        error_norm = sqrt(error_norm);
        y_norm = sqrt(y_norm);
        double relative_error = error_norm / y_norm;

        // Check that relative error is small
        if (relative_error > 1e-4) {
            printf("FAILED: Relative error too large: %.2e\n", relative_error);
            return false;
        }

        printf("PASSED (relative error: %.2e)\n", relative_error);
        return true;

    } catch (const exception& e) {
        printf("FAILED: Exception thrown: %s\n", e.what());
        return false;
    }
}

/**
 * Test 3: Different kernel functions
 *
 * Tests H-matrix with various kernel types to ensure generality
 */
bool test_hmatrix_kernels() {
    printf("  Test 3: Different kernel functions... ");

    try {
        // Create points
        uint n_points = 500;
        vector<vector<double>> points(n_points);

        for (uint i = 0; i < n_points; i++) {
            double angle = 2.0 * M_PI * i / n_points;
            points[i] = {cos(angle), sin(angle)};
        }

        // Test 1: Newton kernel (1/r)
        auto kernel_newton = [](const double* x, const double* y) -> double {
            double dx = x[0] - y[0];
            double dy = x[1] - y[1];
            double r = sqrt(dx * dx + dy * dy);
            return (r < 1e-10) ? 0.0 : 1.0 / r;
        };

        HMatrix hm_newton(points, kernel_newton, 4, 32, 1e-6, 2.0);

        // Test 2: Gaussian kernel
        auto kernel_gauss = [](const double* x, const double* y) -> double {
            double dx = x[0] - y[0];
            double dy = x[1] - y[1];
            return exp(-0.5 * (dx * dx + dy * dy));
        };

        HMatrix hm_gauss(points, kernel_gauss, 4, 32, 1e-6, 2.0);

        // Verify both matrices were created successfully
        if (hm_newton.size() != n_points || hm_gauss.size() != n_points) {
            printf("FAILED: Wrong matrix size\n");
            return false;
        }

        // Quick matvec test
        vector<double> x(n_points, 1.0);
        vector<double> y1 = hm_newton.matvec(x);
        vector<double> y2 = hm_gauss.matvec(x);

        printf("PASSED\n");
        return true;

    } catch (const exception& e) {
        printf("FAILED: Exception thrown: %s\n", e.what());
        return false;
    }
}

/**
 * Main test runner for H-matrix tests
 */
bool hmatrix_tests() {
    printf("\nRunning H-matrix tests\n");

    bool all_tests[3] = {
        test_hmatrix_creation(),
        test_hmatrix_matvec(),
        test_hmatrix_kernels()
    };

    return print_test_results(all_tests, 3, "H-matrix tests");
}
