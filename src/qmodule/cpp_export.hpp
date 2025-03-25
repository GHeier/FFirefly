#pragma once

#include <complex>
#include "../objects/CMField/fields.hpp"  // Include field class

float epsilon_julia(int n, float kx, float ky = 0, float kz = 0);

CMF* create_CMF();
// Load CMF from a file
CMF* load_CMF(const char* filename);

// Create a new field instance and return a pointer
Field_C* create_field_CS();

// Load field from a file
Field_C* load_field_cs(const char* filename);

// Call operator() overload for `Vec`
void field_cs_call(Field_C* field, Vec point, float w, std::complex<float>* result);

// Call operator() overload for `float`
void field_cs_call2(Field_C* field, float w, std::complex<float>* result);

void save_data(string filename, vector<Vec> &points, vector<complex<Vec>> &values, int dimension, bool with_w, bool is_complex, bool is_vector);

// Destroy field instance
void destroy_Field_C(Field_C* field);
void destroy_Field_R(Field_R* field);

CMF* create_cmf();
CMF* load_cmf(const char* filename);

void get_points(CMF* cmf, vector<Vec> &points);
void get_w_points(CMF* cmf, vector<float> &w_points);
void get_values(CMF* cmf, vector<complex<Vec>> &values);
void get_domain(CMF* cmf, vector<Vec> &domain);
void get_inv_domain(CMF* cmf, vector<Vec> &inv_domain);
void get_first(CMF* cmf, Vec &first);
void get_num_points(CMF* cmf, int &nx, int &ny, int &nz, int &nw);
void get_dimension(CMF* cmf, int &dimension);
void get_w_max_min(CMF* cmf, float &wmax, float &wmin);
void get_is_complex(CMF* cmf, bool &is_complex);
void get_is_vector(CMF* cmf, bool &is_vector);
void get_with_w(CMF* cmf, bool &with_w);
