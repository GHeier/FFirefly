// Description: Header file for interpolation functions
#pragma once

#include <vector>
#include <complex>
#include "../objects/vec.hpp"

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
