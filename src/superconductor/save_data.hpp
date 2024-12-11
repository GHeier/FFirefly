#pragma once

#include <vector>
#include <complex>

#include "frequency_inclusion.hpp"
#include "../objects/vec.hpp"
#include "../objects/matrix.hpp"
#include "../config/load/cpp_config.hpp"

using std::endl;
using std::cout;
using std::vector;

/**
 * @brief Save the Delta(kx, ky, kz)'s and their eigenvalues to the file
 */
void save(string file_name, float T, vector<Vec> FS, Eigenvector* solutions);
/**
 * @brief Save the Delta(kx, ky, kz)'s and their eigenvalues to the file in the case where the 
 * interaction is not rounded to just the FS.
 */
void save_with_freq(string file_name, float T, vector<vector<Vec>> &freq_FS, Eigenvector* solutions);
/**
 * @brief Save the Fermi surface to the file
 */
void save_FS(vector<Vec> FS);
/**
 * @brief Save the potential vs q to the file
 */
void save_potential_vs_q(vector<Vec> &FS, Matrix &P, string filename);
/**
 * @brief Save the chi vs q to the file
 */
void save_chi_vs_q(const vector<vector<vector<float>>> &cube, vector<Vec> &FS, string filename);

void save_matsubara_cube(const MatCube &cube, float wmin, float wmax, string filename);

