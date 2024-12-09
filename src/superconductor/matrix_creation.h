#pragma once
#ifndef MATRIX_CREATION_H
#define MATRIX_CREATION_H

#include <vector>
#include <unordered_map>

#include "../objects/vec.h"
#include "../objects/matrix.hpp"

void create_P(Matrix &P, vector<Vec> &k);
void create_P_freq(Matrix &P, vector<vector<Vec>> &k, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube2);
int matrix_size_from_freq_FS(vector<vector<Vec>> &freq_FS);
void vector_to_wave(vector<Vec> &FS, Eigenvector *vectors);
void freq_vector_to_wave(vector<vector<Vec>> &freq_FS, Eigenvector *vectors);

#endif
