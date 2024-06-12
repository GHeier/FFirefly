#pragma once
#ifndef CALCULATIONS_H_
#define CALCULATIONS_H_

#include <vector>

#include "cfg.h"
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"

using std::vector;

vector<Eigenvector> power_iteration(Matrix A, double error);
double f_singlet(double x, double T);
double f_singlet_integral(double T);
void create_P(Matrix &P, vector<Vec> &k, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube2);
double f(vector<Vec> k, double T);
double get_Tc(vector<Vec> k);
void vector_to_wave(vector<Vec> &FS, vector<Eigenvector> &vectors);
double get_DOS(vector<Vec> &FS);
double coupling_calc(vector<Vec> &FS, double T);

#endif
