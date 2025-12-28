/**
 * @file h2matrix.hpp
 * @brief H2Matrix class for hierarchical matrix compression using H2Pack
 *
 * Provides a C++ wrapper around H2Pack for efficient representation and
 * matrix-vector multiplication of kernel matrices V(x_i, x_j) where
 * V is a smooth function.
 */

#pragma once

#include <vector>
#include <functional>
#include <memory>

/**
 * @brief Kernel evaluation function type
 *
 * Takes two points (x1, x2) in d-dimensional space and returns V(x1, x2)
 * Signature: (const double* x1, const double* x2, int dim) -> double
 */
using KernelFunction = std::function<double(const double*, const double*, int)>;

/**
 * @brief Hierarchical matrix class using H2Pack compression
 *
 * Compresses kernel matrices V(x_i, x_j) where V is a smooth function
 * of the distance/difference between points. Provides fast matrix-vector
 * multiplication with O(N log N) complexity instead of O(N²).
 *
 * Usage:
 * @code
 *   // Define kernel: V(r) = 1/(1 + r²)
 *   auto kernel = [](const double* x1, const double* x2, int dim) {
 *       double r2 = 0.0;
 *       for (int i = 0; i < dim; i++) {
 *           double d = x1[i] - x2[i];
 *           r2 += d * d;
 *       }
 *       return 1.0 / (1.0 + r2);
 *   };
 *
 *   // Create points
 *   std::vector<std::vector<double>> points = {...};  // N points × dim
 *
 *   // Build H2 matrix
 *   H2Matrix mat(kernel, points, 2, 1e-4);
 *
 *   // Matrix-vector multiply
 *   std::vector<double> x(N, 1.0);
 *   std::vector<double> y = mat.matvec(x);
 * @endcode
 */
class H2Matrix {
public:
    /**
     * @brief Construct H2Matrix from kernel function and points
     *
     * @param kernel Kernel function V(x1, x2)
     * @param points N × dim array of point coordinates
     * @param dim Spatial dimension
     * @param rel_tol Relative tolerance for compression (default: 1e-4)
     */
    H2Matrix(KernelFunction kernel,
             const std::vector<std::vector<double>>& points,
             int dim,
             double rel_tol = 1e-4);

    /**
     * @brief Destructor - cleans up H2Pack structure
     */
    ~H2Matrix();

    // Disable copy (H2Pack structures are non-copyable)
    H2Matrix(const H2Matrix&) = delete;
    H2Matrix& operator=(const H2Matrix&) = delete;

    /**
     * @brief Matrix-vector multiplication: y = A * x
     *
     * @param x Input vector (size N)
     * @return Output vector y (size N)
     */
    std::vector<double> matvec(const std::vector<double>& x) const;

    /**
     * @brief Get number of points (matrix dimension)
     */
    int size() const { return n_points_; }

    /**
     * @brief Get spatial dimension
     */
    int dimension() const { return dim_; }

    /**
     * @brief Get compression statistics
     *
     * @return Approximate compression ratio (original size / compressed size)
     */
    double compression_ratio() const;

    /**
     * @brief Get maximum tree level
     */
    int max_level() const;

    /**
     * @brief Get number of tree nodes
     */
    int num_nodes() const;

private:
    KernelFunction kernel_;      ///< User-provided kernel function (currently unused, default kernel in use)
    void* h2pack_;              ///< H2Pack structure (opaque pointer)
    int n_points_;              ///< Number of points
    int dim_;                   ///< Spatial dimension
    double rel_tol_;            ///< Relative tolerance
    std::vector<double> coords_; ///< Flattened coordinates (dim * n_points)
};
