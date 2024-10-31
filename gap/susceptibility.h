#pragma once
#ifndef SUSCEPTIBILITY_H_
#define SUSCEPTIBILITY_H_

#include <vector>
#include <complex>
#include <functional>
#include "vec.h"
#include "fermi_surface.h"
#include "potential.h"

using namespace std;

/**
 * @brief Fermi-Dirac distribution function
 * 
 * @param E Energy
 * @param T Temperature
 * @return float 
 */
float f(float E, float T);

/**
 * @brief Same as integrand, but without using .freq to define constant surfaces
 */
float ratio(Vec k, Vec q, float w, float T);

/**
 * @brief Integrand of the susceptibility
 *
 * @param k Momentum
 * @param q Momentum transfer
 * @param w Frequency
 * @param T Temperature
 *
 * @return float
 */
float integrand(Vec k, Vec q, float w, float T);

/**
 * @brief Integrate the susceptibility
 *
 * @param q Momentum transfer
 * @param T Temperature
 * @param mu Chemical potential
 * @param w Frequency
 * @param num_points Number of points
 *
 * @return float
 */
float integrate_susceptibility(Vec q, float T, float mu, float w, int num_points);

/**
 * @brief Integrate the complex susceptibility (uses iw instead of w for frequency)
 *
 * @param q Momentum transfer
 * @param T Temperature
 * @param mu Chemical potential
 * @param w Frequency
 * @param num_points Number of points
 *
 * @return complex<float>
 */
complex<float> complex_susceptibility_integration(Vec q, float T, float mu, complex<float> w, int num_points);

/**
 * @brief Integrate a complex function via the trapezoidal method
 *
 * @param q Momentum transfer
 * @param T Temperature
 * @param mu Chemical potential
 * @param w Frequency
 * @param num_points Number of points
 *
 * @return float
 */
complex<float> complex_trapezoidal_integration(const function<complex<float>(float, float, float)> &f, float x0, float x1, float y0, float y1, float z0, float z1, int num_points);

/**
 * @brief Integrate a function via the trapezoidal method
 *
 * @param q Momentum transfer
 * @param T Temperature
 * @param mu Chemical potential
 * @param w Frequency
 * @param num_points Number of points
 *
 * @return float
 */
float trapezoidal_integration(const function<double(float, float, float)> &f, float x0, float x1, float y0, float y1, float z0, float z1, int num_points);
/**
 * @brief Calculate the susceptibility cube
 *
 * This creates a cube in k-space, where at each point the susceptibility integral is evaluated.
 * Linear interpolation is used to calculate the values between the points, giving us the 
 * susceptibility function over all k-space.
 *
 * @param T Temperature
 * @param mu Chemical potential
 * @param w Frequency
 * @param message Message to print
 *
 * @return vector<vector<vector<float>>
 */
vector<vector<vector<float>>> chi_cube(float T, float mu, float w, string message);

/**
 * @brief Calculate the susceptibility from the cube
 *
 * This function takes the susceptibility cube and calculates the value at a given point in k-space.
 * Uses interpolation.
 *
 * @param chi_cube Susceptibility cube
 * @param q Momentum
 *
 * @return float
 */
float calculate_chi_from_cube(const vector<vector<vector<float>>> &chi_cube, Vec q);

/**
 * @brief Converts a point in the first Brillouin zone to the first IBZ
 *
 * This function takes a point in the first Brillouin zone and converts it to the first IBZ.
 * The benefit of this is not doing identical integrals multiple times.
 *
 * @param k Momentum
 * @return Vec
 */
Vec to_IBZ_2(const Vec k);

#endif
