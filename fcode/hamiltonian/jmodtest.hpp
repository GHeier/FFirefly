#pragma once

#include <complex>
#include <iostream>
#include "../objects/cmf_types.hpp"  // Include CMF class

float epsilon_julia(int n, float kx, float ky = 0, float kz = 0);

// Create a new CMF instance and return a pointer
CMF_CS* create_CMF_CS();

// Load CMF from a file
CMF_CS* load_CMF_cs(const char* filename);

// Call operator() overload for `Vec`
void cmf_cs_call(CMF_CS* cmf, Vec point, float w, std::complex<float>* result);

// Call operator() overload for `float`
void cmf_cs_call2(CMF_CS* cmf, float w, std::complex<float>* result);

// Destroy CMF instance
void destroy_CMF(CMF_CS* cmf);

