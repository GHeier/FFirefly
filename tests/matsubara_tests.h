#pragma once
#ifndef MATSUBARA_TESTS_H
#define MATSUBARA_TESTS_H

#include <complex>

#include "../gap/vec.h"

using namespace std;

bool check_matcube_interpolation();
bool compare_real_vs_complex_susceptibility();
bool compare_real_vs_complex_susceptibility_integration(Vec q, float T, float mu, float w, float num_points);

#endif
