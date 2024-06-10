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

void save(string file_name, double T, vector<Vec> k, std::vector<Eigenvector> solutions);
void save_FS(vector<Vec> FS);
void save_potential_vs_q(vector<Vec> &FS, Matrix &P, string filename);
void save_chi_vs_q(const vector<vector<vector<double>>> &cube, vector<Vec> &FS, string filename);

#endif
