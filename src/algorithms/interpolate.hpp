// Description: Header file for interpolation functions
#pragma once

#include <vector>
#include <complex>
#include "src/objects/vec.hpp"

using namespace std;

float sanitize_within_bounds(float value, float min_bound, float max_bound, float tolerance = 1e-5);
int binary_search(float x_val, vector<float> &x);
extern float interpolate_1D(float x_val, float x_min, float x_max, vector<float> &f);
extern float interpolate_2D(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, vector<vector<float>> &f);
extern float interpolate_3D(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, vector<vector<vector<float>>> &f);
extern float interpolate_4D(float x_val, float y_val, float z_val, float w_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, float w_min, float w_max, vector<vector<vector<vector<float>>>> &f);

complex<float> interpolate_1D_complex(float x_val, float x_min, float x_max, vector<complex<float>> &f);
complex<float> interpolate_2D_complex(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, vector<vector<complex<float>>> &f);
complex<float> interpolate_3D_complex(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, vector<vector<vector<complex<float>>> > &f);
complex<float> interpolate_4D_complex(float x_val, float y_val, float z_val, float w_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, float w_min, float w_max, vector<vector<vector<vector<complex<float>>>>> &f);

float interpolate_2D(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, int nx, int ny, vector<float> &f);
float interpolate_3D(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, int nx, int ny, int nz, vector<float> &f);
float interpolate_4D(float x_val, float y_val, float z_val, float w_val, 
        float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, 
        float w_min, float w_max, int nx, int ny, int nz, int nw, vector<float> &f);

complex<Vec> interpolate_1D(float x_val, float x_min, float x_max, vector<complex<Vec>> &f);
complex<Vec> interpolate_2D(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, int nx, int ny, vector<complex<Vec>> &f);
complex<Vec> interpolate_3D(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, int nx, int ny, int nz, vector<complex<Vec>> &f);
complex<Vec> interpolate_4D(float x_val, float y_val, float z_val, float w_val,
        float x_min, float x_max, float y_min, float y_max, float z_min, float z_max,
        float w_min, float w_max, int nx, int ny, int nz, int nw, vector<complex<Vec>> &f);

// Interpolation for indexed fields (vectors and matrices)
using cfloat = std::complex<float>;

// 1D spatial interpolation returning vector
vector<cfloat> interpolate_1D_vec(float x_val, float x_min, float x_max, int nx, const vector<vector<cfloat>>& f);
vector<cfloat> interpolate_2D_vec(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, int nx, int ny, const vector<vector<cfloat>>& f);
vector<cfloat> interpolate_3D_vec(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, int nx, int ny, int nz, const vector<vector<cfloat>>& f);

// Spatial interpolation returning matrix (2D tensor)
vector<vector<cfloat>> interpolate_1D_mat(float x_val, float x_min, float x_max, int nx, const vector<vector<vector<cfloat>>>& f);
vector<vector<cfloat>> interpolate_2D_mat(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, int nx, int ny, const vector<vector<vector<cfloat>>>& f);
vector<vector<cfloat>> interpolate_3D_mat(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, int nx, int ny, int nz, const vector<vector<vector<cfloat>>>& f);

// Spatial interpolation returning 3D tensor
vector<vector<vector<cfloat>>> interpolate_1D_ten3(float x_val, float x_min, float x_max, int nx, const vector<vector<vector<vector<cfloat>>>>& f);
vector<vector<vector<cfloat>>> interpolate_2D_ten3(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, int nx, int ny, const vector<vector<vector<vector<cfloat>>>>& f);
vector<vector<vector<cfloat>>> interpolate_3D_ten3(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, int nx, int ny, int nz, const vector<vector<vector<vector<cfloat>>>>& f);

// Spatial interpolation returning 4D tensor
vector<vector<vector<vector<cfloat>>>> interpolate_1D_ten4(float x_val, float x_min, float x_max, int nx, const vector<vector<vector<vector<vector<cfloat>>>>>& f);
vector<vector<vector<vector<cfloat>>>> interpolate_2D_ten4(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, int nx, int ny, const vector<vector<vector<vector<vector<cfloat>>>>>& f);
vector<vector<vector<vector<cfloat>>>> interpolate_3D_ten4(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, int nx, int ny, int nz, const vector<vector<vector<vector<vector<cfloat>>>>>& f);
