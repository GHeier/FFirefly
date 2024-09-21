#pragma once
#ifndef FREQUENCY_INCLUSION_H
#define FREQUENCY_INCLUSION_H

#include <complex>
#include "fermi_surface.h"
#include "calculations.h"

/**
 * @brief Calculates the size of the matrix from the surfaces defined, including those beyond
 * the Fermi Surface
 * 
 * @param freq_FS The frequencies and Fermi surface
 * @return int The size of the matrix
 */
int matrix_size_from_freq_FS(vector<vector<Vec>> &freq_FS);

/**
 * @brief Calculates the energy surfaces off of the Fermi Surface using the tetrahedron method
 * 
 * @param mu The chemical potential
 * @return vector<vector<Vec>> The surfaces
 */
vector<vector<Vec>> freq_tetrahedron_method(float mu);

/**
 * @brief Creates the P matrix for these surfaces
 *
 * @param P The matrix to be filled
 * @param k The surfaces 
 * @param T The temperature
 * @param chi_cube2 The susceptibility cube
 */
void create_P_freq(Matrix &P, vector<vector<Vec>> &k, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube2);

/**
 * @brief Creates the susceptibility cubes for each difference in energy between surfaces
 * 
 * @param T The temperature
 * @param mu The chemical potential
 * @return unordered_map<float, vector<vector<vector<float>>>>
 */
unordered_map <float, vector<vector<vector<float>>>> chi_cube_freq(float T, float mu);

vector<vector<vector<vector<complex<float>>>>> create_matsubara_cube(float T, float MU, int m_pts, int w_pts, float w_min, float w_max);

/**
 * @brief Is the numerator of the susceptibility integral
 *
 * @param k The momentum
 * @param q The momentum transfer
 * @param w The frequency
 * @param T The temperature
 *
 * @return float The numerator
 */
float integrand(Vec k, Vec q, float w, float T);

/**
 * @brief The denominator of the susceptibility integral
 *
 * This is used in order to define the peaks and zeros of the susceptibility, that way we can
 * evenly integrate around them without having to worry about the instability of the poles.
 *
 * @param k The momentum
 * @param q The momentum transfer
 *
 * @return float The denominator
 */
float denominator(Vec k, Vec q);

/**
 * @brief The derivative of the denominator
 *
 * This is used in order to integrate along the surfaces (ie surface space), rather than k-space
 * integral of dS*dk^2 rather than dk^3
 *
 * @param k The momentum
 * @param q The momentum transfer
 *
 * @return float The derivative
 */
float denominator_diff(Vec k, Vec q);

/**
 * @brief The spacing of the surfaces to be integrated over
 *
 * @param spacing The spacing vector
 * @param w The frequency
 * @param a The lower bound
 * @param b The upper bound
 * @param pts The number of points
 */
void get_spacing_vec(vector<float> &spacing, float w, float a, float b, int pts);

/**
 * @brief The constants for the spacing curve. It is defined as a linear y=mx+b curve
 *
 * @param w The frequency
 * @param a The lower bound
 * @param b The upper bound
 * @param A The constant for the curve
 * @param upr The upper bound
 * @param lwr The lower bound
 */
void get_spacing_curve_consts(float w, float a, float b, float &A, float &upr, float &lwr);

/**
 * @brief Interpolate chi from the map for a given frequency
 *
 * @param chi_cube_map The map of chi values
 * @param q The momentum transfer
 * @param w The frequency
 *
 * @return float The interpolated chi value
 */
float calculate_chi_from_cube_map(const unordered_map<float, vector<vector<vector<float>>>> &chi_cube_map, Vec q, float w);

/**
 * @brief The bounds of the surfaces
 *
 * The bounds are calculated by essentially random sampling, and finding the highest and lowest.
 * This value is then reduced by 1% in order to allow for an actual surface to exist there, rather
 * than a single point
 *
 * @param q The momentum transfer
 * @param upper The upper bound
 * @param lower The lower bound
 * @param func The function to be integrated
 */
void get_bounds(Vec q, float &upper, float &lower, float (*func)(Vec k, Vec q));

#endif
