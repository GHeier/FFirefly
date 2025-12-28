/**
 * Standalone test for H2Matrix matvec operation
 * Compares H2Matrix matvec with dense matrix-vector multiplication
 */
#include <iostream>
#include <vector>
#include <cmath>
#include "src/objects/h2matrix.hpp"

using namespace std;

int main() {
    cout << "H2Matrix matvec test" << endl;
    cout << "===================" << endl;

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

    cout << "Building H2Matrix with " << n_points << " points..." << endl;
    H2Matrix mat(kernel, points, dim, 1e-3);

    cout << "  Size: " << mat.size() << endl;
    cout << "  Dimension: " << mat.dimension() << endl;
    cout << "  Max level: " << mat.max_level() << endl;
    cout << "  Num nodes: " << mat.num_nodes() << endl;
    cout << "  Compression ratio: " << mat.compression_ratio() << endl;

    // Create test vector: x[i] = sin(i * 0.1)
    vector<double> x(n_points);
    for (int i = 0; i < n_points; i++) {
        x[i] = sin(i * 0.1);
    }

    cout << "\nComputing H2 matvec..." << endl;
    vector<double> y_h2 = mat.matvec(x);

    cout << "Computing dense matvec for comparison..." << endl;
    vector<double> y_dense(n_points, 0.0);
    for (int i = 0; i < n_points; i++) {
        for (int j = 0; j < n_points; j++) {
            double val = kernel(points[i].data(), points[j].data(), dim);
            y_dense[i] += val * x[j];
        }
    }

    cout << "\nComparing results..." << endl;
    double max_error = 0.0;
    double rms_error = 0.0;
    for (int i = 0; i < n_points; i++) {
        double err = fabs(y_h2[i] - y_dense[i]);
        max_error = max(max_error, err);
        rms_error += err * err;
    }
    rms_error = sqrt(rms_error / n_points);

    cout << "  Max error: " << scientific << max_error << endl;
    cout << "  RMS error: " << rms_error << endl;
    cout << "  Tolerance: 1e-3" << endl;

    // Sample comparison at a few points
    cout << "\nSample values (first 5 points):" << endl;
    cout << "  i    H2 result    Dense result    Error" << endl;
    for (int i = 0; i < min(5, n_points); i++) {
        cout << "  " << i << "    " << y_h2[i] << "    " << y_dense[i]
             << "    " << fabs(y_h2[i] - y_dense[i]) << endl;
    }

    // Check accuracy
    double tolerance = 1e-3;
    if (max_error > tolerance) {
        cout << "\nFAILED: max error exceeds tolerance!" << endl;
        return 1;
    }

    if (rms_error > tolerance) {
        cout << "\nFAILED: RMS error exceeds tolerance!" << endl;
        return 1;
    }

    cout << "\nSUCCESS: H2Matrix matvec agrees with dense within tolerance!" << endl;
    return 0;
}
