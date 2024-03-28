#pragma once
#ifndef FREQUENCY_INCLUSION_H
#define FREQUENCY_INCLUSION_H

#include "fermi_surface.h"
#include "calculations.h"


vector<vector<Vec>> freq_tetrahedron_method(double mu);
void create_P_freq(Matrix &P, vector<vector<Vec>> &k, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube2);
unordered_map<double, vector<vector<vector<double>>>> chi_cube_freq(double T, double mu, double DOS);
double f_singlet_integral_test(double T);

#endif
