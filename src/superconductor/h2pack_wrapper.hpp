/**
 * H2Pack wrapper for BCS pairing matrix compression
 *
 * This provides an interface to use H2Pack hierarchical matrix compression
 * for the BCS pairing matrix P, which can achieve O(N) storage and O(N) matvec.
 *
 * Author: Griffin Heier
 */
#pragma once

#include <vector>
#include "../objects/matrix.hpp"
#include "../objects/vec.hpp"

// Workaround: H2Pack headers need cstring but don't include it
#include <cstring>

// H2Pack headers
#include <H2Pack_typedef.h>

/**
 * H2Pack matrix wrapper class
 *
 * Provides hierarchical low-rank compression of the BCS pairing matrix
 */
class H2PackMatrix {
public:
    H2Pack_s* h2pack;          // H2Pack structure
    std::vector<Vec> points;   // k-points on Fermi surface
    int n_points;              // Number of points
    int dim;                   // Dimension (2 or 3)
    double rel_tol;            // Relative tolerance for compression

    // Statistics
    long original_nnz;         // Original matrix size (N²)
    long compressed_nnz;       // Compressed storage
    double compression_ratio;  // Compression factor
    int max_rank;              // Maximum block rank
    double build_time;         // Construction time (seconds)
    double matvec_time;        // Matvec time (seconds)

    H2PackMatrix(int dim = 2, double rel_tol = 1e-6);
    ~H2PackMatrix();

    /**
     * Build H2 matrix from BCS pairing kernel
     *
     * @param FS Fermi surface k-points
     * @param renorm Renormalization factor
     */
    void build_from_kernel(const std::vector<Vec>& FS, float renorm);

    /**
     * Build H2 matrix from existing dense matrix P
     *
     * @param P Dense pairing matrix
     * @param FS Fermi surface k-points
     */
    void build_from_matrix(const Matrix& P, const std::vector<Vec>& FS);

    /**
     * Matrix-vector multiplication: y = H2 * x
     */
    void matvec(const std::vector<float>& x, std::vector<float>& y);

    /**
     * Print compression statistics
     */
    void print_stats() const;

    /**
     * Get compression ratio
     */
    double get_compression_ratio() const { return compression_ratio; }

    /**
     * Get maximum rank in hierarchical representation
     */
    int get_max_rank() const { return max_rank; }

private:
    void compute_stats();
};


/**
 * BCS pairing kernel function for H2Pack
 *
 * Computes V(k1 - k2) for two points on the Fermi surface
 * This is the kernel function passed to H2Pack
 *
 * Signature matches kernel_eval_fptr from H2Pack
 */
void bcs_pairing_kernel(
    const double* coord0, const int ld0, const int n0,
    const double* coord1, const int ld1, const int n1,
    const void* krnl_param, double* mat, const int ldm
);

