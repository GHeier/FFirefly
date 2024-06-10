#pragma once
#ifndef FREQUENCY_INCLUSION_H
#define FREQUENCY_INCLUSION_H

#include "fermi_surface.h"
#include "calculations.h"


vector<vector<Vec>> freq_tetrahedron_method(double mu);
void create_P_freq(Matrix &P, vector<vector<Vec>> &k, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube2);
unordered_map<double, vector<vector<vector<double>>>> chi_cube_freq(double T, double mu, double DOS);
double integrand(Vec k, Vec q, double w, double T);
double denominator(Vec k, Vec q);
double denominator_diff(Vec k, Vec q);
void get_spacing_vec(vector<double> &spacing, double w, double a, double b, int pts);
void get_spacing_curve_consts(double w, double a, double b, double &A, double &upr, double &lwr);
double calculate_chi_from_cube_map(const unordered_map<double, vector<vector<vector<double>>>> &chi_cube_map, Vec q, double w);
void get_bounds3(Vec q, double &upper, double &lower, double (*func)(Vec k, Vec q));

#endif
