#pragma once
#ifndef PLOT_H
#define PLOT_H

#include <complex>

#include "../gap/vec.h"

void plot_interpolated_chi();
void plot_real_susceptibility_integration(float w);
void plot_complex_susceptibility_integration(complex<float> w);
void plot_complex_susceptibility_integration_v_w(Vec q);
void plot_real_trapezoidal_susceptibility_integration(float w);
void plot_matsubara_cube_v_q(complex<float> w);
void plot_matsubara_cube_v_w(Vec q);

void plot_surfaces(Vec q, float T, float w);

void plot_surfaces2(Vec q, float T, float w);

#endif
