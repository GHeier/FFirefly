#include <iostream>
#include <vector>
#include <math.h>
#include <cassert>

#include "solver.hpp"
#include "matrix_creation.hpp"
#include "utilities.hpp"
#include "../hamiltonian/potential.hpp"
#include "../hamiltonian/band_structure.hpp"
#include "../algorithms/integration.hpp"
#include "../objects/matrix.hpp"
#include "../objects/vec.hpp"
#include "../config/load/cpp_config.hpp"

using namespace std;

// Create V matrix
// Picks the potential based on the global variable "interaction"
void create_P(Matrix &P, vector<Vec> &k) {
    cout << "Creating P Matrix\n";
    for (int i = 0; i < P.size; i++) {
        Vec k1 = k[i];
        #pragma omp parallel for
        for (int j = 0; j < P.size; j++) {
            Vec k2 = k[j];
            P(i,j) = (float)(-pow(k1.area/vp(k1.n, k1),0.5) * V(k1, k2) * pow(k2.area/vp(k2.n, k2),0.5));
            assert(isnan(P(i,j)) == false);
        }
        progress_bar(1.0 * i / P.size);
    }
    P *= (2 / pow(2*M_PI, dim));
    cout << "\nP Matrix Created\n";
}

// Creates the P matrix based around the multiple energy surfaces calculated above
void create_P_freq(Matrix &P, vector<vector<Vec>> &k, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube2) {

    cout << "Creating P Matrix with frequency\n";
    for (int i = 0; i < k.size(); i++) {

        int ind1 = 0;
        for (int temp = 0; temp < i; temp++)
            ind1 += k[temp].size();

        for (int j = 0; j < k[i].size(); j++) {
            Vec k1 = k[i][j];
            for (int x = 0; x < k.size(); x++) {

                int ind2 = 0;
                for (int temp = 0; temp < x; temp++)
                    ind2 += k[temp].size();

                #pragma omp parallel for
                for (int y = 0; y < k[x].size(); y++) {
                    Vec k2 = k[x][y];
                    float d1 = pow(k1.area/vp(k1.n, k1),0.5); 
                    float d2 = pow(k2.area/vp(k2.n, k2),0.5); 
                    // f * d_epsilon
                    float fde1 = f_singlet(wc * points[l-1][i], T) * weights[l-1][i];
                    float fde2 = f_singlet(wc * points[l-1][x], T) * weights[l-1][x];
                    float w = wc * (points[l-1][x] - points[l-1][i]);

                    P(ind1 + j,ind2 + y) = (float)(- d1 * d2 * pow(fde1*fde2,0.5) * V(k1, k2)); 
                }
            }
            string message = "Portion " + to_string(i) + " of " + to_string(k.size());
            progress_bar(1.0 * (ind1 + j) / P.size, message);
        }
    }
    cout << "P Matrix Created\n";

    //return P * 2 * wc / (l * k_size);
    P *= wc * (2 / pow(2*M_PI, dim)); 
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
            vectors[i][j] /= pow(k.area/vp(k.n, k),0.5);
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
                vectors[x][ind] /= pow(k.area/vp(k.n, k),0.5);
                ind++;
            }
        }
    }
}

