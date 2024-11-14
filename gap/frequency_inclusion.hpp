#pragma once
#ifndef FREQUENCY_INCLUSION_H
#define FREQUENCY_INCLUSION_H

#include <fstream>
#include <complex>
#include "calculations.h"
#include "interpolate.h"
#include "vec.h"
#include "susceptibility.h"

struct MatCube {
    vector<vector<vector<vector<complex<float>>>>> cube;
    float x_min;
    float x_max;
    float y_min;
    float y_max;
    float z_min;
    float z_max;
    float w_min;
    float w_max;
    int num_integral_pts;

    // Constructor
    MatCube(int mx, int my, int mz, int w_pts, float x_min, float x_max, 
            float y_min, float y_max, float z_min, float z_max, 
            float w_min, float w_max, int num_integral_pts)
        : x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), z_min(z_min), 
        z_max(z_max), w_min(w_min), w_max(w_max), num_integral_pts(num_integral_pts)
    {
        // Initialize 'cube' to the specified dimensions
        cube = vector<vector<vector<vector<complex<float>>>>>(
            mx, vector<vector<vector<complex<float>>>>(
                my, vector<vector<complex<float>>>(
                    mz, vector<complex<float>>(w_pts)
                )
            )
        );
    }

    // Accessor
    complex<float> operator()(Vec q, complex<float> w) {
        q = to_IBZ_2(q);
        return interpolate_4D_complex(q(0), q(1), q(2), w.imag(), 
                x_min, x_max, y_min, y_max, z_min, z_max, w_min, w_max, cube);
    }

    // Data Loading
    void load_data(string filename) {
        ifstream file(filename);
        float x_val, y_val, z_val, w_val, re, im;
        // Now fill the cube with the values. The cube is filled over the entirety of the BZ
        for (int i = 0; i < cube.size(); i++) {
            for (int j = 0; j < cube[0].size(); j++) {
                for (int k = 0; k < cube[0][0].size(); k++) {
                    for (int l = 0; l < cube[0][0][0].size(); l++) {
                        file >> x_val >> y_val >> z_val >> w_val >> re >> im;
                        cube[i][j][k][l] = complex<float>(re, im);
                    }
                }
            }
        }
    }
};

/**
 * @brief Calculates the size of the matrix from the surfaces defined, including those beyond
 * the Fermi Surface
 * 
 * @param freq_FS The frequencies and Fermi surface
 * @return int The size of the matrix
 */
int matrix_size_from_freq_FS(vector<vector<Vec>> &freq_FS);

/**
 * @brief Calculates the energy surfaces off of the Fermi Surface using the tetrahedron method
 * 
 * @param mu The chemical potential
 * @return vector<vector<Vec>> The surfaces
 */
vector<vector<Vec>> freq_tetrahedron_method(float mu);

/**
 * @brief Creates the P matrix for these surfaces
 *
 * @param P The matrix to be filled
 * @param k The surfaces 
 * @param T The temperature
 * @param chi_cube2 The susceptibility cube
 */
void create_P_freq(Matrix &P, vector<vector<Vec>> &k, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube2);

/**
 * @brief Creates the susceptibility cubes for each difference in energy between surfaces
 * 
 * @param T The temperature
 * @param mu The chemical potential
 * @return unordered_map<float, vector<vector<vector<float>>>>
 */
unordered_map <float, vector<vector<vector<float>>>> chi_cube_freq(float T, float mu);

/**
 * @brief Creates the susceptibility cube in 4D using the Matsubara frequencies
 * 
 * @param T The temperature
 * @param mu The chemical potential
 * @param m_pts The number of k points
 * @param w_pts The number of Matsubara frequencies
 * @param w_min The minimum Matsubara frequency
 * @param w_max The maximum Matsubara frequency
 *
 * @return MatCube The Matsubara cube
 */
MatCube create_matsubara_cube(float T, float MU, int m_pts, int w_pts, float w_min, float w_max, int num_integral_pts);

/**
 * @brief Interpolate chi from the map for a given frequency
 *
 * @param chi_cube_map The map of chi values
 * @param q The momentum transfer
 * @param w The frequency
 *
 * @return float The interpolated chi value
 */
float calculate_chi_from_cube_map(const unordered_map<float, vector<vector<vector<float>>>> &chi_cube_map, Vec q, float w);

#endif
