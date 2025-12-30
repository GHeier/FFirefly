/*
 * Field interpolation function declarations.
 * These functions are used by DataEvaluator for multi-dimensional interpolation.
 *
 * Author: Griffin Heier
 */
#pragma once

#include <complex>
#include <vector>

#include "src/objects/vec.hpp"

using namespace std;

// 1D interpolation with frequency
complex<Vec> CMF_search_1d(float w_val, vector<float> &w_points,
                           vector<complex<Vec>> &f);

// 2D interpolation (1 spatial + 1 frequency)
complex<Vec> CMF_search_2d(float x_val, float w_val, int nx,
                           vector<float> &w_points, vector<complex<Vec>> &f);

// 3D interpolation (2 spatial + 1 frequency)
complex<Vec> CMF_search_3d(float x_val, float y_val, float w_val, int nx,
                           int ny, vector<float> &w_points,
                           vector<complex<Vec>> &f);

// 4D interpolation (3 spatial + 1 frequency)
complex<Vec> CMF_search_4d(float x_val, float y_val, float z_val, float w_val,
                           int nx, int ny, int nz, vector<float> &w_points,
                           vector<complex<Vec>> &f);
