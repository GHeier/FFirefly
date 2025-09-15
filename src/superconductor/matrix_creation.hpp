#pragma once

#include <unordered_map>
#include <vector>

#include "../objects/matrix.hpp"
#include "../objects/vec.hpp"

void create_P(Matrix &P, vector<Vec> &k);
void create_P_freq(Matrix &P, vector<vector<Vec>> &k, double T);
int matrix_size_from_freq_FS(vector<vector<Vec>> &freq_FS);
void vector_to_wave(vector<Vec> &FS, Eigenvector *vectors);
void freq_vector_to_wave(vector<vector<Vec>> &freq_FS, Eigenvector *vectors);
