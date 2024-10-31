#pragma once
#ifndef MATSUBARA_TESTS_H
#define MATSUBARA_TESTS_H

#include <complex>

using namespace std;

void plot_matsubara_cube(complex<float> w);
bool compare_real_vs_complex_susceptibility();
bool test_real_integration();
bool test_complex_integration();
void plot_real_susceptibility();

#endif
