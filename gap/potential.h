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
double potential_scalapino_cube(Vec k1, Vec k2, double T, const vector<vector<vector<double>>> &chi_cube);
double potential_scalapino_triplet(Vec k1, Vec k2, double T, const vector<vector<vector<double>>> &chi_cube);

// Scalapino Potential Section
double get_k(double i, double n);
double f(double E, double T);
double ratio(Vec q, Vec k, double T, double mu);
double ratio2(Vec q, Vec k, double T, double mu);
double chi_trapezoidal(Vec q, double T, double mu, int num_points);
double integrate_susceptibility(Vec q, double T, double mu);
double trapezoidal_integration(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int num_points);
double trap_cube(auto &f, double x0, double x1, double y0, double y1, double z0, double z1);
double trap_8_cubes(auto &f, double x0, double x1, double y0, double y1, double z0, double z1);
double adaptive_trapezoidal(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int xdivs, int ydivs, int zdivs, double error_relative);
double iteratively_splitting_cubes(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, double error_total, double error_relative);
vector<vector<vector<double>>> chi_cube(double T, double mu, double DOS);
vector<vector<vector<double>>> chi_cube2(double T, double mu, double DOS);
double calculate_chi_from_cube(const vector<vector<vector<double>>> &chi_cube, Vec q);
Vec to_IBZ_2(const Vec k);
Vec to_IBZ_spherical(const Vec k);


#endif
