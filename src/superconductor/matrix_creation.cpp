#include <cassert>
#include <iostream>
#include <math.h>
#include <vector>

#include "../algorithms/integration.hpp"
#include "../config/load/cpp_config.hpp"
#include "../hamiltonian/band_structure.hpp"
#include "../hamiltonian/interaction.hpp"
#include "../objects/CMField/vertex.hpp"
#include "../objects/matrix.hpp"
#include "../objects/vec.hpp"
#include "cfg.hpp"
#include "matrix_creation.hpp"
#include "solver.hpp"
#include "utilities.hpp"
#include "h2pack_wrapper.hpp"

using namespace std;

// Create V matrix
// Picks the potential based on the global variable "interaction"
void create_P(Matrix &P, vector<Vec> &k) {
    Vertex V_func;
    cout << "Creating P Matrix\n";
    for (int i = 0; i < P.size; i++) {
        Vec k1 = k[i];
        float f1 = pow(k1.area / vp(k1.n, k1), 0.5);
        #pragma omp parallel for
        for (int j = 0; j < P.size; j++) {
            Vec k2 = k[j];
            float f2 = pow(k2.area / vp(k2.n, k2), 0.5);
            P(i, j) = -f1 * f2 * (V_func(k1 - k2, 0).real() + V_func(k1 + k2, 0).real()) / 2.0;
            //P(i, j) = f1 * f2 * (cos(k1.x) - cos(k1.y)) * (cos(k2.x) - cos(k2.y));
            assert(isnan(P(i, j)) == false);
        }
        progress_bar(1.0 * i / (P.size - 1));
    }
    P *= (1 / pow(2 * M_PI, dim));
    cout << "\nP Matrix Created\n";
}

// Creates the P matrix based around the multiple energy surfaces calculated
// above
void create_P_freq(Matrix &P, vector<vector<Vec>> &k, double T) {
    Vertex V_func;
    cout << "Creating P Matrix with frequency\n";
    for (int i = 0; i < k.size(); i++) {

        int ind1 = 0;
        for (int temp = 0; temp < i; temp++)
            ind1 += k[temp].size();

        for (int j = 0; j < k[i].size(); j++) {
            Vec k1 = k[i][j];
            float d1 = k1.area / vp(k1.n, k1);
            float w1 = wc * points[l - 1][i];
            float fde1 = f_singlet(w1, T) * weights[l - 1][i];
            for (int x = 0; x < k.size(); x++) {

                int ind2 = 0;
                for (int temp = 0; temp < x; temp++)
                    ind2 += k[temp].size();

#pragma omp parallel for
                for (int y = 0; y < k[x].size(); y++) {
                    Vec k2 = k[x][y];
                    float d2 = k2.area / vp(k2.n, k2);
                    float w2 = wc * points[l - 1][x];
                    // f * d_epsilon
                    float fde2 = f_singlet(w2, T) * weights[l - 1][x];

                    double V_int = (V_func(k1 - k2, w1 - w2).real() +
                                    V_func(k1 + k2, w1 + w2).real()) /
                                   2.0;
                    double prefactor = pow(d1 * d2 * fde1 * fde2, 0.5);
                    P(ind1 + j, ind2 + y) = (float)(-prefactor * V_int);
                }
            }
            string message =
                "Portion " + to_string(i) + " of " + to_string(k.size());
            progress_bar(1.0 * (ind1 + j) / (P.size - 1), message);
        }
    }
    cout << "P Matrix Created\n";
    P *= wc * (1 / pow(2 * M_PI, dim));
}

int matrix_size_from_freq_FS(vector<vector<Vec>> &freq_FS) {
    int size = 0;
    for (int i = 0; i < freq_FS.size(); i++) {
        size += freq_FS[i].size();
    }
    return size;
}

// Un-shifting the area-shifted eigenvectors in order to find wavefunction
void vector_to_wave(vector<Vec> &FS, Eigenvector *vectors) {
    for (unsigned int i = 0; i < num_eigenvalues_to_save; i++) {
        for (unsigned int j = 0; j < vectors[i].size; j++) {
            Vec k = FS[j];
            vectors[i][j] /= pow(k.area / vp(k.n, k), 0.5);
        }
    }
}

void freq_vector_to_wave(vector<vector<Vec>> &freq_FS, Eigenvector *vectors) {
    int size = matrix_size_from_freq_FS(freq_FS);
    for (unsigned int x = 0; x < num_eigenvalues_to_save; x++) {
        int ind = 0;
        for (unsigned int i = 0; i < freq_FS.size(); i++) {
            for (unsigned int j = 0; j < freq_FS[i].size(); j++) {
                Vec k = freq_FS[i][j];
                vectors[x][ind] /= pow(k.area / vp(k.n, k), 0.5);
                ind++;
            }
        }
    }
}

/**
 * Create and test H2Pack hierarchical matrix compression of P
 *
 * This function demonstrates H2Pack compression on the BCS pairing matrix.
 * It builds both the dense matrix P and the H2Pack compressed version,
 * then compares them and outputs compression statistics.
 *
 * @param P Dense pairing matrix (output)
 * @param k Fermi surface k-points
 * @param renorm Renormalization factor
 */
void create_P_h2pack(Matrix &P, vector<Vec> &k, float renorm) {
    cout << "\n" << string(70, '=') << endl;
    cout << "H2Pack Hierarchical Matrix Compression Test" << endl;
    cout << string(70, '=') << endl;

    // First create the dense P matrix using standard method
    create_P(P, k);

    // Create H2Pack compressed version
    H2PackMatrix h2_matrix(dim, 1e-6);  // 2D/3D, relative tolerance 1e-6

    cout << "\nBuilding H2Pack representation from kernel..." << endl;
    h2_matrix.build_from_kernel(k, renorm);

    // Alternative: build from dense matrix
    // cout << "\nBuilding H2Pack representation from dense matrix..." << endl;
    // h2_matrix.build_from_matrix(P, k);

    // Print compression statistics
    h2_matrix.print_stats();

    // Test matrix-vector multiplication
    cout << "\nTesting H2Pack matrix-vector multiplication..." << endl;
    int n = k.size();
    vector<float> x(n, 1.0);  // Test vector (all ones)
    vector<float> y_dense(n, 0.0);
    vector<float> y_h2(n, 0.0);

    // Dense matvec
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            y_dense[i] += P(i, j) * x[j];
        }
    }

    // H2Pack matvec
    h2_matrix.matvec(x, y_h2);

    // Compare results
    float max_diff = 0.0;
    float rel_err = 0.0;
    for (int i = 0; i < n; i++) {
        float diff = abs(y_dense[i] - y_h2[i]);
        max_diff = max(max_diff, diff);
        rel_err += diff * diff;
    }
    rel_err = sqrt(rel_err / n);

    cout << "\nMatrix-vector multiplication comparison:" << endl;
    cout << "  Max absolute error:  " << max_diff << endl;
    cout << "  RMS relative error:  " << rel_err << endl;

    // Summary
    cout << "\n" << string(70, '=') << endl;
    cout << "Summary:" << endl;
    cout << string(70, '=') << endl;
    cout << "  Original matrix:       " << n << " x " << n << endl;
    cout << "  Original storage:      " << n*n << " elements" << endl;
    cout << "  Compressed storage:    " << h2_matrix.compressed_nnz << " elements" << endl;
    cout << "  Compression ratio:     " << h2_matrix.get_compression_ratio() << "x" << endl;
    cout << "  Maximum block rank:    " << h2_matrix.get_max_rank() << endl;
    cout << "  Build time:            " << h2_matrix.build_time << " seconds" << endl;
    cout << "  Matvec accuracy:       " << rel_err << " (relative error)" << endl;
    cout << string(70, '=') << "\n" << endl;

    // Memory comparison
    float dense_mb = (n * n * sizeof(float)) / (1024.0 * 1024.0);
    float h2_mb = (h2_matrix.compressed_nnz * sizeof(double)) / (1024.0 * 1024.0);
    cout << "Memory savings: " << dense_mb - h2_mb << " MB ("
         << 100.0 * (1.0 - h2_mb/dense_mb) << "%)" << endl;

    cout << "\nNote: To enable H2Pack, uncomment the H2Pack code in" << endl;
    cout << "      h2pack_wrapper.cpp and link against libH2Pack.a" << endl;
}
