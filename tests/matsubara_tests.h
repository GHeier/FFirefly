#pragma once
#ifndef MATSUBARA_TESTS_H
#define MATSUBARA_TESTS_H

#include <complex>

#include "../gap/vec.h"

using namespace std;

bool check_matcube_interpolation();
bool compare_real_vs_complex_susceptibility();
bool compare_real_vs_complex_susceptibility_integration(Vec q, float T, float mu, float w, float num_points);
bool check_nonzero_imaginary_part();
void complex_integration_convergence_test();
void analytical_integration_convergence_test();
bool test_real_integration();
bool test_complex_integration();

#endif
