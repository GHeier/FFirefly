#pragma once

#include <unordered_map>
#include <vector>

#include "../objects/vec.hpp"
#include "../objects/matrix.hpp"

float f_singlet(float x, float T);
float f_singlet_integral(float T);
float get_Tc_FS_only(double eig);
float f(vector<Vec> k, float T,
        const unordered_map<float, vector<vector<vector<float>>>> &cube_map);
float get_Tc(
    vector<Vec> k,
    const unordered_map<float, vector<vector<vector<float>>>> &cube_map);
float get_renormalization(vector<Vec> &FS);
float get_renormalization_off_FS(vector<vector<Vec>> &FS);
vector<float> matrix_projections(vector<Vec> &FS, Matrix &P, float renorm);
