#pragma once
#ifndef FREQUENCY_INCLUSION_H
#define FREQUENCY_INCLUSION_H

#include "fermi_surface.h"
#include "calculations.h"


vector<vector<Vec>> freq_tetrahedron_method(double mu);
MatrixXd create_P_freq(vector<vector<Vec>> &k, double T, const vector<vector<vector<double>>> &chi_cube2);
MatrixXd create_P_freq2(vector<vector<Vec>> &k, double T, const vector<vector<vector<double>>> &chi_cube2);
double f_singlet_integral_test(double T);

#endif
