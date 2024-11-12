#pragma once
#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <vector>
#include <complex>
#include <functional>

#include "vec.h"

extern int s_div;
extern int s_pts;
extern float *weights[5];
extern float *points[5];

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

void get_spacing_curve_consts(float w, float a, float b, float &A, float &upr, float &lwr);
void get_spacing_vec(vector<float> &spacing, float w, float a, float b, int pts);
void get_surface_transformed_bounds(float &upper, float &lower, function<float(Vec)> func);

float surface_transform_integral(function<float(Vec)> integrand,
        function<float(Vec)> func,
        function<float(Vec)> func_diff,
        vector<float> &svals);

#endif
