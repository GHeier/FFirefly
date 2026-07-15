#include "preloaded_hamiltonians.hpp"
#include "band_structure.hpp"
#include "src/objects/vec.hpp"
#include "src/config/load/cpp_config.hpp"
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

using namespace std;
using namespace std::complex_literals;

using mat = vector<vector<complex<float>>>;

// Create a 1×1 matrix with a single real value
mat make_1x1_matrix(float value) {
    mat H(1);
    H[0].resize(1);
    H[0][0] = complex<float>(value, 0.0);
    return H;
}

// Create an N×N matrix initialized to zero
static mat make_matrix(int n) {
    mat H(n);
    for (int i = 0; i < n; i++) {
        H[i].resize(n, complex<float>(0.0, 0.0));
    }
    return H;
}

void fill_conjugate(mat &H) {
    for (int i = 0; i < H.size(); i++) {
        for (int j = 0; j < H[i].size(); j++) {
            if (H[i][j] == complex<float>(0.0,0.0)) {
                H[i][j] = conj(H[j][i]);
            }
        }
    }
}

mat hamiltonian_emery(float kx, float ky, float kz,
                  float eps_d, float eps_p,
                  float t_pd, float t_pp,
                  float mu) {

    // Create 3×3 matrix
    auto H = make_matrix(3);

    // Form factors
    // f_x(k) = 2i·sin(kx/2)
    // f_y(k) = 2i·sin(ky/2)
    complex<float> f_x(0.0, 2.0 * sin(kx/2.0));
    complex<float> f_y(0.0, 2.0 * sin(ky/2.0));

    // g(k) = 2[cos(kx/2) - cos(ky/2)]
    complex<float> g(2.0 * (cos(kx/2.0) - cos(ky/2.0)), 0.0);

    // Diagonal elements (on-site energies)
    H[0][0] = complex<float>(eps_d - mu, 0.0);  // Cu d orbital
    H[1][1] = complex<float>(eps_p - mu, 0.0);  // O px orbital
    H[2][2] = complex<float>(eps_p - mu, 0.0);  // O py orbital

    // Off-diagonal elements (hopping terms)
    // Cu-O hopping
    H[0][1] = t_pd * f_x;
    H[0][2] = t_pd * f_y;
    H[1][0] = conj(H[0][1]);  // Hermiticity
    H[2][0] = conj(H[0][2]);

    // O-O hopping
    H[1][2] = t_pp * g;
    H[2][1] = conj(H[1][2]);  // Hermiticity

    return H;
}

// ============================================================================
// Utility Functions
// ============================================================================

mat hamiltonian_dp(float kx, float ky, float kz, float eps_dx2y2, float eps_dz, float delta_dp) {
    float eps_p = eps_dx2y2 - delta_dp;

    auto H = make_matrix(4);
    H[0][0] = eps_dx2y2;
    H[1][0] = 0.0;
    H[1][1] = eps_dz - 2 * t4 * (cos(kx) + cos(ky));
    H[2][0] = 2i * (double)t0 * sin(0.5 * kx);
    H[2][1] = -2i * (double)t3 * sin(0.5 * kx);
    H[2][2] = eps_p + 2 * t2 * cos(kx) + 2 * t5 * (cos(kx + ky) + cos(kx - ky));
    H[3][0] = -2i * (double)t0 * sin(0.5 * ky);
    H[3][1] = -2i * (double)t3 * sin(0.5 * ky);
    H[3][2] = 2 * t1 * (cos(0.5 * kx + 0.5 * ky) - cos(0.5 * kx - 0.5 * ky));
    H[3][3] = eps_p + 2 * t2 * cos(ky) + 2 * t5 * (cos(kx + ky) + cos(kx - ky));

    fill_conjugate(H);
    return H;
}

mat H(Vec k) {
    if (hamiltonian == "tight_binding") {
        float e = epsilon(1, k);
        return make_1x1_matrix(e);
    }
    else if (hamiltonian == "emery") {
        float eps_d = 0.0;
        float eps_p = 3.6;
        float t_pd = 1.3;
        float t_pp = 0.65;
        float mu = 0.0;
        return hamiltonian_emery(k.x, k.y, k.z, eps_d, eps_p, t_pd, t_pp, mu);
    }
    else if (hamiltonian == "dp") {
        // LSCO values from HIROSHI WATANABE et al. PHYSICAL REVIEW RESEARCH 3, 033157 (2021)
        float eps_dx2y2 = -0.87;
        float eps_dz = -0.11;
        float delta_dp = 2.26;
        return hamiltonian_dp(k.x, k.y, k.z, eps_dx2y2, eps_dz, delta_dp);
    }
    else {
        printf("Hamiltonian type: %s with %d bands is not supported.\n", hamiltonian.c_str(), nbnd);
        throw runtime_error("Unsupported number of bands for H(k)");
    }
}
