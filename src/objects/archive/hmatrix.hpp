#pragma once

#include <vector>
#include <functional>
#include <complex>
#include <stdexcept>
#include <cmath>

extern "C" {
#include <hmatrix.h>
#include <kernelmatrix.h>
#include <basic.h>
#include <cluster.h>
#include <clustergeometry.h>
#include <block.h>
#include <aca.h>
#include <matrixnorms.h>
#include <avector.h>
#include <amatrix.h>
}

using namespace std;

/**
 * HMatrix class
 *
 * C++ wrapper for H2Lib's hierarchical matrix (hmatrix) implementation.
 * H-matrices provide a memory-efficient and fast representation of dense matrices
 * arising from kernel functions or other problems with inherent low-rank structure.
 *
 * Key features:
 * - Hierarchical block structure with admissible (low-rank) and inadmissible (dense) blocks
 * - Compressed storage: O(n log n) or O(n) instead of O(n^2) for dense matrices
 * - Fast matrix-vector multiplication: O(n log n) instead of O(n^2)
 * - Supports various kernel functions (Newton, exponential, custom)
 *
 * Usage example:
 *   // Create points on a circle
 *   vector<vector<double>> points = ...;
 *
 *   // Define kernel function k(x,y) = exp(-|x-y|^2)
 *   auto kernel = [](const double* x, const double* y) {
 *       double dx = x[0] - y[0], dy = x[1] - y[1];
 *       return exp(-(dx*dx + dy*dy));
 *   };
 *
 *   // Build H-matrix
 *   HMatrix hm(points, kernel, 4, 32, 1e-6, 2.0);
 *
 *   // Matrix-vector multiplication
 *   vector<double> x(points.size(), 1.0);
 *   vector<double> y = hm.matvec(x);
 */
class HMatrix {
public:
    /**
     * Constructor: Build H-matrix from kernel function and points
     *
     * @param points List of spatial points (each point is a vector of coordinates)
     * @param kernel Kernel function k(x,y) that takes two coordinate arrays
     * @param interpolation_order Interpolation order m (number of interpolation points)
     * @param leafsize Maximum cluster size for leaf nodes
     * @param eps Tolerance for low-rank approximation (e.g., 1e-6)
     * @param eta Admissibility parameter (typically 2.0)
     */
    HMatrix(const vector<vector<double>>& points,
            function<double(const double*, const double*)> kernel,
            uint interpolation_order = 4,
            uint leafsize = 32,
            double eps = 1e-6,
            double eta = 2.0);

    /**
     * Destructor: Clean up H2Lib structures
     */
    ~HMatrix();

    // Delete copy constructor and assignment to prevent double-free
    HMatrix(const HMatrix&) = delete;
    HMatrix& operator=(const HMatrix&) = delete;

    /**
     * Matrix-vector multiplication: y = H * x
     *
     * @param x Input vector (size must match number of points)
     * @return Output vector y
     */
    vector<double> matvec(const vector<double>& x) const;

    /**
     * Matrix-vector multiplication with alpha scaling: y = alpha * H * x + y
     *
     * @param x Input vector
     * @param y Output vector (modified in place)
     * @param alpha Scaling factor (default 1.0)
     * @param transpose Use H^T instead of H (default false)
     */
    void matvec(const vector<double>& x, vector<double>& y,
                double alpha = 1.0, bool transpose = false) const;

    /**
     * Convert to dense matrix (for testing/comparison)
     * WARNING: This defeats the purpose of H-matrices and uses O(n^2) memory!
     *
     * @return Dense matrix as vector of vectors (row-major)
     */
    vector<vector<double>> to_dense() const;

    /**
     * Get number of rows/columns (matrix is square)
     */
    uint size() const { return n_points; }

    /**
     * Get memory usage in bytes
     */
    size_t memory_size() const;

    /**
     * Get memory usage for nearfield (dense blocks) in bytes
     */
    size_t nearfield_size() const;

    /**
     * Get memory usage for farfield (low-rank blocks) in bytes
     */
    size_t farfield_size() const;

    /**
     * Get compression ratio compared to dense matrix
     */
    double compression_ratio() const;

    /**
     * Approximate spectral norm ||H||_2
     */
    double norm() const;

    /**
     * Print statistics about the H-matrix structure
     */
    void print_stats() const;

private:
    // H2Lib objects
    phmatrix hm;                    // The H-matrix
    pclustergeometry cg;           // Cluster geometry
    pcluster root;                 // Cluster tree
    pblock broot;                  // Block tree
    pkernelmatrix km;              // Kernel matrix description

    // Parameters
    uint n_points;                 // Number of points
    uint dim;                      // Spatial dimension
    uint m;                        // Interpolation order
    double tolerance;              // Approximation tolerance

    // Helper function to fill H-matrix from kernel matrix
    static void fill_hmatrix_kernelmatrix(pckernelmatrix km, phmatrix h);

    // Store kernel function (need to keep it alive for callbacks)
    function<double(const double*, const double*)> kernel_func;

    // Static wrapper for C callback (using :: to disambiguate from std::real)
    static field kernel_callback(const ::real* xx, const ::real* yy, void* data);
};
