#pragma once
#ifndef MATRIX_CREATION_H
#define MATRIX_CREATION_H

#include <vector>
#include <unordered_map>

#include "../objects/vec.h"
#include "../objects/matrix.hpp"

void create_P(Matrix &P, vector<Vec> &k, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube2);

#endif
