#pragma once
#ifndef FREQUENCY_INCLUSION_H
#define FREQUENCY_INCLUSION_H

#include "fermi_surface.h"
#include "calculations.h"


vector<vector<Vec>> freq_tetrahedron_method(double mu);
void create_P_freq(Matrix P, vector<vector<Vec>> &k, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube2);
unordered_map<double, vector<vector<vector<double>>>> chi_cube_freq(double T, double mu, double DOS);
double f_singlet_integral_test(double T);
double zero_temp_func(double w, double dE);
double nonzero_ratio(double w, double dE, Vec k, Vec q, double T);
vector<vector<Vec>> get_singularity_freq_surfaces(double (*func)(Vec k, Vec q), Vec q, double center, double width);
vector<vector<Vec>> get_smooth_freq_surfaces(double (*func)(Vec k, Vec q), Vec q, double lower, double upper, int num_points);
double bound_sign(Vec k, Vec q);
double other_bound_sign(Vec k, Vec q);
double gauss_chi_sum(vector<vector<Vec>> &surfaces, Vec q, double w, double T, double width);
double bound_chi_sum(Vec q, double w, double T, double de, double width, double b, double a);
double chi_integrate_freq(Vec q, double width, double w, double T);
void get_bounds(Vec q, double w, double &upper, double &lower);
double imaginary_chi_integrate(Vec q, double w);
void shift_layers(vector<Vec> &layers, Vec shift, double &max, double &min);
double modified_e_diff(Vec k, Vec q);
double chi_ep_integrate(Vec q, double w, double T);
void get_bounds2(Vec q, double &upper, double &lower);
double bound_chi_sum2(Vec q, double w, double T, int pts, double b, double a);
double sphere_func(Vec k, Vec q);
void get_bounds3(Vec q, double &upper, double &lower, double (*func)(Vec k, Vec q));
double singularity_func(Vec k, Vec q);
double denominator(Vec k, Vec q);
double denominator_diff(Vec k, Vec q);
void get_spacing_vec(vector<double> &spacing, double w, double a, double b, int pts);
double bound_chi_sum4(Vec q, double w, double T, int pts, double b, double a, double (*func)(Vec k, Vec q), double (*func_diff)(Vec k, Vec q));
void get_spacing_curve_consts(double w, double a, double b, double &A, double &upr, double &lwr);
double comparison_integral(Vec q, double w, double b, double a, int pts, double (*func)(Vec k, Vec q));
double close_to_zero(Vec q, double w, double T, int pts);
double num_states(double w, double T, int pts);

#endif
