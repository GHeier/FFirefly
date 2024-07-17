#pragma once
#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <vector>
#include "vec.h"
#include "fermi_surface.h"

using namespace std;

struct indices_values {
    vector<int> indices;
    float value;
};

/**
 * @brief Constant potential
 * 
 * @param k1 
 * @param k2 
 * @return float 
 */
float potential_const(Vec k1, Vec k2);

/**
 * @brief Test potential
 * 
 * @param k1 
 * @param k2 
 * @return float 
 */
float potential_test(Vec k1, Vec k2);

/**
 * @brief Phonon Coulomb potential
 * 
 * @param q 
 * @return float 
 */
float phonon_coulomb(Vec q);

/**
 * @brief Scalapino singlet potential (FLEX) done via direct integration
 * 
 * @param k1 
 * @param k2 
 * @param T 
 * @return float 
 */
float potential_scal(Vec k1, Vec k2, float T);

/**
 * @brief Scalapino singlet potential (FLEX) done via chi cube interpolation
 * 
 * @param k1 
 * @param k2 
 * @param w 
 * @param T 
 * @param chi_map 
 * @return float 
 */
float potential_scalapino_cube(Vec k1, Vec k2, float w, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_map);

/**
 * @brief Scalapino triplet potential (FLEX) done via chi cube interpolation
 * 
 * @param k1 
 * @param k2 
 * @param w 
 * @param T 
 * @param chi_map 
 * @return float 
 */
float potential_scalapino_triplet(Vec k1, Vec k2, float w, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_map);

// Scalapino Potential Section
<<<<<<< HEAD
=======
double get_k(double i, double n);
double f(double E, double T);
double ratio(Vec q, Vec k, double T, double mu, double w);
double modified_ratio(Vec q, Vec k, double T, double mu, double w, double delta);
double imaginary_ratio(Vec q, Vec k, double T, double mu, double w, double eta);
double imaginary_integration(Vec q, double T, double mu, double w, int num_points, double eta);
double chi_trapezoidal(Vec q, double T, double mu, double w, int num_points);
double integrate_susceptibility(Vec q, double T, double mu, double w, int num_points);
double trapezoidal_integration(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int num_points);
double modified_integral_wrapper(Vec q, double T, double mu, double w, double delta, int num_points);
double modified_integral_wrapper_1D(double a, double b, double delta, int num_points);
double modified_integration(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int num_points, auto &surface, double delta, auto &avg);
double modified_integration_1D(auto &f, double x0, double x1, int num_points, auto &surface, double delta, auto &avg);
double trap_cube(auto &f, double x0, double x1, double y0, double y1, double z0, double z1);
double trap_8_cubes(auto &f, double x0, double x1, double y0, double y1, double z0, double z1);
double adaptive_trapezoidal(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, int xdivs, int ydivs, int zdivs, double error_relative, int num_splits);
double iteratively_splitting_cubes(auto &f, double x0, double x1, double y0, double y1, double z0, double z1, double error_total, double error_relative);
vector<vector<vector<double>>> chi_cube(double T, double mu, double DOS, double w);
double calculate_chi_from_cube(const vector<vector<vector<double>>> &chi_cube, Vec q);
Vec to_IBZ_2(const Vec k);
Vec to_IBZ_spherical(const Vec k);
>>>>>>> origin/main


#endif
