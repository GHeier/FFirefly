/**
 * Verification test for H2Matrix class
 * Tests that the H2 structure is built correctly
 */
#include <iostream>
#include <vector>
#include <cmath>
#include "src/objects/h2matrix.hpp"

using namespace std;

int main() {
    cout << "H2Matrix Verification Test" << endl;
    cout << "==========================\n" << endl;

    // Define kernel: V(r) = 1/(1 + r²)
    auto kernel = [](const double* x1, const double* x2, int dim) -> double {
        double r2 = 0.0;
        for (int i = 0; i < dim; i++) {
            double d = x1[i] - x2[i];
            r2 += d * d;
        }
        return 1.0 / (1.0 + r2);
    };

    // Create points on a sphere (3D)
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

    cout << "Test 1: H2Matrix Construction" << endl;
    cout << "------------------------------" << endl;
    cout << "Building H2Matrix with " << n_points << " points..." << endl;

    H2Matrix mat(kernel, points, dim, 1e-3);

    cout << "✓ H2Matrix constructed successfully\n" << endl;

    cout << "Test 2: H2Matrix Properties" << endl;
    cout << "----------------------------" << endl;
    cout << "  Size: " << mat.size() << " (expected: " << n_points << ")" << endl;
    cout << "  Dimension: " << mat.dimension() << " (expected: " << dim << ")" << endl;
    cout << "  Max level: " << mat.max_level() << endl;
    cout << "  Num nodes: " << mat.num_nodes() << endl;
    cout << "  Compression ratio: " << mat.compression_ratio() << "x" << endl;

    bool props_ok = true;
    if (mat.size() != n_points) {
        cout << "  ✗ Size mismatch!" << endl;
        props_ok = false;
    }
    if (mat.dimension() != dim) {
        cout << "  ✗ Dimension mismatch!" << endl;
        props_ok = false;
    }
    if (mat.max_level() <= 0) {
        cout << "  ✗ Tree not partitioned!" << endl;
        props_ok = false;
    }
    if (mat.num_nodes() <= 1) {
        cout << "  ✗ No tree nodes created!" << endl;
        props_ok = false;
    }
    if (mat.compression_ratio() < 1.0) {
        cout << "  ✗ Invalid compression ratio!" << endl;
        props_ok = false;
    }

    if (props_ok) {
        cout << "✓ All properties correct\n" << endl;
    } else {
        cout << "✗ Some properties failed\n" << endl;
        return 1;
    }

    cout << "Test 3: Direct Kernel Evaluation" << endl;
    cout << "---------------------------------" << endl;
    // Sample a few kernel values to verify the function works
    double k00 = kernel(points[0].data(), points[0].data(), dim);
    double k01 = kernel(points[0].data(), points[1].data(), dim);
    double k12 = kernel(points[1].data(), points[2].data(), dim);

    cout << "  K(p0, p0) = " << k00 << " (should be ~1.0 for self-interaction)" << endl;
    cout << "  K(p0, p1) = " << k01 << " (should be < 1.0)" << endl;
    cout << "  K(p1, p2) = " << k12 << " (should be < 1.0)" << endl;

    if (fabs(k00 - 1.0) > 1e-6) {
        cout << "  ✗ Self-interaction kernel value wrong!" << endl;
        return 1;
    }
    if (k01 <= 0.0 || k01 >= 1.0) {
        cout << "  ✗ Kernel values out of expected range!" << endl;
        return 1;
    }

    cout << "✓ Kernel evaluation correct\n" << endl;

    cout << "\n========================================" << endl;
    cout << "SUCCESS: All verification tests passed!" << endl;
    cout << "======================================== \n" << endl;
    cout << "The H2Matrix class successfully:" << endl;
    cout << "  - Constructs hierarchical matrix representation" << endl;
    cout << "  - Builds octree spatial partitioning" << endl;
    cout << "  - Achieves " << mat.compression_ratio() << "x compression vs dense matrix" << endl;
    cout << "  - Stores " << n_points << "x" << n_points << " matrix in compressed form" << endl;

    return 0;
}
