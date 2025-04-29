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

using namespace std;

// Create V matrix
// Picks the potential based on the global variable "interaction"
void create_P(Matrix &P, vector<Vec> &k) {
    Vertex V_func;
    cout << "Creating P Matrix\n";
    for (int i = 0; i < P.size; i++) {
        Vec k1 = k[i];
#pragma omp parallel for
        for (int j = 0; j < P.size; j++) {
            Vec k2 = k[j];
            P(i, j) = (float)(-pow(k1.area / vp(k1.n, k1), 0.5) *
                              (V_func(k1 - k2, 0, "up", "down").real() +
                               V_func(k1 + k2, 0, "up", "down").real()) *
                              pow(k2.area / vp(k2.n, k2), 0.5)) /
                      2.0;
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
