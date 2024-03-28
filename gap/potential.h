#pragma once
#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <vector>
#include "vec.h"
#include "fermi_surface.h"

using namespace std;

struct indices_values {
    vector<int> indices;
    double value;
};

double potential_const(Vec k1, Vec k2);
double potential_test(Vec k1, Vec k2);
double phonon_coulomb(Vec q);
double potential_scal(Vec k1, Vec k2, double T);
double potential_scalapino_cube(Vec k1, Vec k2, double w, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_map);
double potential_scalapino_triplet(Vec k1, Vec k2, double w, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_map);

// Scalapino Potential Section
double get_k(double i, double n);
double f(double E, double T);
double ratio(Vec q, Vec k, double T, double mu, double w);
double ratio2(Vec q, Vec k, double T, double mu);
double chi_trapezoidal(Vec q, double T, double mu, int num_points);
double integrate_susceptibility(Vec q, double T, double mu, double w);
double trapezoidal_integration(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int num_points);
double trap_cube(auto &f, double x0, double x1, double y0, double y1, double z0, double z1);
double trap_8_cubes(auto &f, double x0, double x1, double y0, double y1, double z0, double z1);
double adaptive_trapezoidal(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int xdivs, int ydivs, int zdivs, double error_relative, int num_splits);
double iteratively_splitting_cubes(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, double error_total, double error_relative);
vector<vector<vector<double>>> chi_cube(double T, double mu, double DOS, double w);
double calculate_chi_from_cube(const vector<vector<vector<double>>> &chi_cube, Vec q);
Vec to_IBZ_2(const Vec k);
Vec to_IBZ_spherical(const Vec k);


#endif
