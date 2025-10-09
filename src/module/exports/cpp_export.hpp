#pragma once

#include "../../objects/CMField/fields.hpp" // Include field class
#include "../../objects/CMField/bands.hpp"
#include <complex>

float epsilon_export(int n, float kx, float ky = 0, float kz = 0);
Bands *Bands_export0();
float Bands_operator_export0(Bands *obj, int n, const float *point, int len);
void Bands_operator_export0_numpy(Bands *obj, int n, const float *points, int num_points, int len, float *output);
float Bands_operator_export1(Bands *obj, int n, const float *point, int len);
void Bands_operator_export1_numpy(Bands *obj, int n, const float *points, int num_points, int len, float *output);
int Field_R_nbnd_export0(Field_R* a);
int Field_C_nbnd_export0(Field_C* a);
Field_R *Field_R_export0();
Field_R *Field_R_export2(const char *filename);

float Field_R_operator_export0(Field_R *obj, float w);
float Field_R_operator_export1(Field_R *obj, int n, float w);
float Field_R_operator_export2(Field_R *obj, const float *point, int len,
                               float w);
float Field_R_operator_export3(Field_R *obj, int n, const float *point, int len,
                               float w);
void Field_R_operator_export_list(Field_R *obj, const float *points, int num_points, int len,
                                   float w, float *output);

Field_C *Field_C_export();
Field_C *Field_C_export2(const char *filename);

void Field_C_operator_export0(Field_C *obj, float w, float *real_result,
                              float *imag_result);
void Field_C_operator_export1(Field_C *obj, int n, float w, float *real_result,
                              float *imag_result);
void Field_C_operator_export2(Field_C *obj, const float *point, int len,
                              float w, float *real_result, float *imag_result);
void Field_C_operator_export3(Field_C *obj, int n, const float *point, int len,
                              float w, float *real_result, float *imag_result);
void Field_C_operator_export_list(Field_C *obj, const float *points, int num_points, int len,
                                  float w, float *real_output, float *imag_output);

Field_RM *Field_RM_export0();
Field_RM *Field_RM_export2(const char *filename);
void Field_RM_operator_export0(Field_RM *obj, const float *point, int len,
                               float w, float **result, int *matrix_size);
void Field_RM_operator_export_list(Field_RM *obj, const float *points, int num_points, int len,
                                   float w, float *output, int *matrix_size);

Field_CM *Field_CM_export0();
Field_CM *Field_CM_export2(const char *filename);
void Field_CM_operator_export0(Field_CM *obj, const float *point, int len,
                               float w, float **real_result, float **imag_result, int *matrix_size);
void Field_CM_operator_export_list(Field_CM *obj, const float *points, int num_points, int len,
                                   float w, float *real_output, float *imag_output, int *matrix_size);


// Create a new field instance and return a pointer
Field_C *create_field_CS();

// Load field from a file
Field_C *load_field_cs(const char *filename);

// Call operator() overload for `Vec`
void field_cs_call(Field_C *field, Vec point, float w,
                   std::complex<float> *result);

// Call operator() overload for `float`
void field_cs_call2(Field_C *field, float w, std::complex<float> *result);

void save_data(string filename, vector<Vec> &points,
               vector<complex<Vec>> &values, int dimension, bool with_w,
               bool is_complex, bool is_vector);

// Destroy field instance
void destroy_Field_C(Field_C *field);
void destroy_Field_R(Field_R *field);
void destroy_Field_RM(Field_RM *field);
void destroy_Field_CM(Field_CM *field);
void destroy_Bands(Bands *field);



// Save data exports
void save_data_scalar_export0(const char *filename, const float *data_interleaved,
                               int total_size, bool is_complex,
                               const int *mesh, int mesh_size,
                               const float *domain_flat, int domain_rows, int domain_cols,
                               const float *w_points, int w_size);

void save_data_vector_export0(const char *filename, const float *data_interleaved,
                               int nk, int vec_len, bool is_complex,
                               const int *mesh, int mesh_size,
                               const float *domain_flat, int domain_rows, int domain_cols,
                               const float *w_points, int w_size);

void save_data_matrix_export0(const char *filename, const float *data_interleaved,
                               int num_matrices, int mat_dim, bool is_complex,
                               const int *mesh, int mesh_size,
                               const float *domain_flat, int domain_rows, int domain_cols,
                               const float *w_points, int w_size);
