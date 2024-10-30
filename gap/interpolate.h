// Description: Header file for interpolation functions
#pragma once
#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <vector>

using namespace std;

float interpolate_1D(float x_val, float x_min, float x_max, vector<float> &f);
float interpolate_2D(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, vector<vector<float>> &f);
float interpolate_3D(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, vector<vector<vector<float>>> &f);
float interpolate_4D(float x_val, float y_val, float z_val, float w_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, float w_min, float w_max, vector<vector<vector<vector<float>>>> &f);

#endif
