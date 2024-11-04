#pragma once
#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <vector>
#include <complex>
#include <functional>

#include "vec.h"

extern int s_div;
extern int s_pts;

void get_spacing_curve_consts2(float w, float a, float b, float &A, float &upr, float &lwr);
void get_spacing_vec2(vector<float> &spacing, float w, float a, float b, int pts);
void get_surface_transformed_bounds(float &upper, float &lower, function<float(Vec)> func);

float surface_transform_integral(function<float(Vec)> integrand,
        function<float(Vec)> func,
        function<float(Vec)> func_diff,
        vector<float> &svals);

#endif
