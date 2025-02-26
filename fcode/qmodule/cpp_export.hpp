#pragma once

#include <complex>
#include <iostream>
#include "../objects/CondensedMatterField/fields.hpp"  // Include CMF class

float epsilon_julia(int n, float kx, float ky = 0, float kz = 0);

// Create a new CMF instance and return a pointer
Field_C* create_CMF_CS();

// Load CMF from a file
Field_C* load_CMF_cs(const char* filename);

// Call operator() overload for `Vec`
void cmf_cs_call(Field_C* cmf, Vec point, float w, std::complex<float>* result);

// Call operator() overload for `float`
void cmf_cs_call2(Field_C* cmf, float w, std::complex<float>* result);

void save_data(string filename, vector<Vec> &points, vector<complex<Vec>> &values, int dimension, bool with_w, bool is_complex, bool is_vector);

// Destroy CMF instance
void destroy_CMF(Field_C* cmf);

