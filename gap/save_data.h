#pragma once
#ifndef SAVE_DATA_H_
#define SAVE_DATA_H_

#include <vector>

#include "calculations.h"
#include "vec.h"
#include "matrix.hpp"
#include "cfg.h"

using std::endl;
using std::cout;
using std::vector;

extern int iteration;
extern float mu;

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
void save_FS(vector<Vec> FS, int iteration, float mu);
/**
 * @brief Save the potential vs q to the file
 */
void save_potential_vs_q(vector<Vec> &FS, Matrix &P, string filename);
/**
 * @brief Save the chi vs q to the file
 */
void save_chi_vs_q(const vector<vector<vector<float>>> &cube, vector<Vec> &FS, string filename);

void fermi_velocity_average(float mu); // averages. being tested because of failure to run on recursion
void interaction_potential_average(float mu);


void save_vp_dir(float mu); // saves vp_file in directory

void save_V_dir(float mu); // saves V_file in dir

void archive_and_clear_data(const std::string& data_dir); // archives and clears the directories used for the recursive data saving. not in save_data_recursive.cpp because one might desire to archive on every iteration


#endif
