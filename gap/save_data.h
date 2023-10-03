#pragma once
#ifndef SAVE_DATA_H_
#define SAVE_DATA_H_

#include <vector>
#include <Eigen/Dense>

#include "calculations.h"
#include "vec.h"
#include "cfg.h"

using namespace Eigen;
using std::endl;
using std::cout;
using std::vector;

void save(string file_name, double T, vector<Vec> k, std::vector<EigAndVec> solutions);
void save_FS(vector<Vec> FS);
void save_potential_vs_q(vector<Vec> &FS, MatrixXd &P, string filename);
void save_chi_vs_q(const vector<vector<vector<double>>> &cube, vector<Vec> &FS, string filename);

#endif
