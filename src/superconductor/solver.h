#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <unordered_map>

#include "../objects/vec.h"


float f_singlet(float x, float T);
float f_singlet_integral(float T);
float f(vector<Vec> k, float T, const unordered_map<float, vector<vector<vector<float>>> > &cube_map);
float get_Tc(vector<Vec> k, const unordered_map<float, vector<vector<vector<float>>> > &cube_map);
float get_DOS(vector<Vec> &FS);
float coupling_calc(vector<Vec> &FS, float T);

#endif
