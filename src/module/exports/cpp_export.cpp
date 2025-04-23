#include <string>

#include "../../config/load/c_config.h"
#include "../../config/load/cpp_config.hpp"
#include "../../hamiltonian/band_structure.hpp"
#include "../../objects/CMField/bands.hpp"   // Include CMF class
#include "../../objects/CMField/cmfield.hpp" // Include CMF class
#include "../../objects/CMField/fields.hpp"  // Include CMF class
#include "../../objects/CMField/vertex.hpp"  // Include CMF class
#include "../../objects/vec.hpp"

extern "C" float epsilon_export0(int n, float *k, int size) {
  Vec kvec;
  for (int i = 0; i < size; i++) {
    kvec(i) = k[i];
  }
  return epsilon(n, kvec);
}

extern "C" void load_config_export0(const char *filename) {
  read_c_config(filename);
  load_cpp_config();
}

extern "C" {

Bands *Bands_export0() { return new Bands(); }
Vertex *Vertex_export0() { return new Vertex(); }

Field_C *Field_C_export0() { return new Field_C(); }
Field_C *Field_C_export1(CMField cmf) { return new Field_C(cmf); }
Field_C *Field_C_export2(const char *filename) { return new Field_C(filename); }

Field_R *Field_R_export0() { return new Field_R(); }
Field_R *Field_R_export1(CMField cmf) { return new Field_R(cmf); }
Field_R *Field_R_export2(const char *filename) { return new Field_R(filename); }

void Field_C_operator_export0(Field_C *obj, float w, float *real_result,
                              float *imag_result) {
  complex<float> r = obj->operator()(w);
  *real_result = real(r);
  *imag_result = imag(r);
}
void Field_C_operator_export1(Field_C *obj, int n, float w, float *real_result,
                              float *imag_result) {
  complex<float> r = obj->operator()(n, w);
  *real_result = real(r);
  *imag_result = imag(r);
}
void Field_C_operator_export2(Field_C *obj, const float *point, int len,
                              float w, float *real_result, float *imag_result) {
  Vec v(point, len);
  complex<float> r = obj->operator()(v, w);
  *real_result = real(r);
  *imag_result = imag(r);
}
void Field_C_operator_export3(Field_C *obj, int n, const float *point, int len,
                              float w, float *real_result, float *imag_result) {
  Vec v(point, len);
  complex<float> r = obj->operator()(n, v, w);
  *real_result = real(r);
  *imag_result = imag(r);
}

float Field_R_operator_export0(Field_R *obj, float w) {
  return obj->operator()(w);
}
float Field_R_operator_export1(Field_R *obj, int n, float w) {
  return obj->operator()(n, w);
}
float Field_R_operator_export2(Field_R *obj, const float *point, int len,
                               float w) {
  Vec v(point, len);
  return obj->operator()(v, w);
}
float Field_R_operator_export3(Field_R *obj, int n, const float *point, int len,
                               float w) {
  Vec v(point, len);
  return obj->operator()(n, v, w);
}

// Create a new CMF instance and return a pointer
CMField *create_CMField() { return new CMField(); }

// Load CMField from a file
CMField *load_CMField(const char *filename) {
  CMField *cmf = new CMField();
  *cmf = *load_CMField(filename);
  return cmf;
}

void cmf_save(CMField *cmf, const char *filename) {
  save_CMField(filename, *cmf);
}

// Destroy CMField instance
void destroy_Field_C(Field_C *field) { delete field; }

void destroy_Field_R(Field_R *field) { delete field; }

void get_points(CMField *cmf, vector<Vec> &points) {
  points = cmf->data.points;
}
void get_w_points(CMField *cmf, vector<float> &w_points) {
  w_points = cmf->data.w_points;
}
void get_values(CMField *cmf, vector<complex<Vec>> &values) {
  values = cmf->data.values;
}
void get_domain(CMField *cmf, vector<Vec> &domain) { domain = cmf->domain; }
void get_inv_domain(CMField *cmf, vector<Vec> &inv_domain) {
  inv_domain = cmf->inv_domain;
}
void get_first(CMField *cmf, Vec &first) { first = cmf->first; }
void get_num_points(CMField *cmf, int &nx, int &ny, int &nz, int &nw) {
  nx = cmf->nx;
  ny = cmf->ny;
  nz = cmf->nz;
  nw = cmf->nw;
}
void get_dimension(CMField *cmf, int &dimension) {
  dimension = cmf->data.dimension;
}
void get_w_max_min(CMField *cmf, float &wmax, float &wmin) {
  wmax = cmf->wmax;
  wmin = cmf->wmin;
}
void get_is_complex(CMField *cmf, bool &is_complex) {
  is_complex = cmf->data.is_complex;
}
void get_is_vector(CMField *cmf, bool &is_vector) {
  is_vector = cmf->data.is_vector;
}
void get_with_w(CMField *cmf, bool &with_w) { with_w = cmf->data.with_w; }
}
