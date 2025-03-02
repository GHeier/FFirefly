#include <string>

#include "../objects/CondensedMatterField/fields.hpp"  // Include CMF class
#include "../objects/CondensedMatterField/cmf.hpp"  // Include CMF class
#include "../objects/vec.hpp"
#include "../config/load/c_config.h"
#include "../config/load/cpp_config.hpp"
#include "../hamiltonian/band_structure.hpp"

extern "C" float epsilon_julia(int n, float *k, int size) {
    Vec kvec;
    for (int i = 0; i < size; i++) {
        kvec(i) = k[i];
    }
    return epsilon(n, kvec);
}

extern "C" void load_config_julia(const char *filename) {
    read_c_config(filename);
    load_cpp_config();
}


extern "C" {

// Create a new CMF instance and return a pointer
CMF* create_CMF() {
    return new CMF();
}

// Load CMF from a file
CMF* load_CMF(const char* filename) {
    CMF *cmf = new CMF();
    *cmf = load_CMF_from_file(filename);
    return cmf;
}

Field_C* create_Field_C() {
    return new Field_C();
}

// Load CMF from a file
Field_C* load_Field_C(const char* filename) {
    return new Field_C(filename);
}

Field_R* create_Field_R() {
    return new Field_R();
}

Field_R* load_Field_R(const char* filename) {
    return new Field_R(filename);
}

// Call operator() overload for `Vec`
void Field_C_call(Field_C* field, const float *point, float w, int len, float *real_result, float *imag_result) {
    Vec q; q.dimension = len;
    for (int i = 0; i < len; i++) {
        q(i) = point[i];
    }
    complex<float> result = field->operator()(q, w);
    *real_result = real(result);
    *imag_result = imag(result);
}

// Call operator() overload for `float`
void Field_C_call_w(Field_C* field, float w, float *real_result, float *imag_result) {
    complex<float> result = field->operator()(w);
    *real_result = real(result);
    *imag_result = imag(result);
}

void Field_R_call(Field_R* field, const float *point, float w, int len, float *result) {
    Vec q; q.dimension = len;
    for (int i = 0; i < len; i++) {
        q(i) = point[i];
    }
    float temp = field->operator()(q, w);
    *result = temp;
}

void Field_R_call_w(Field_R* field, float w, float *result) {
    float temp = field->operator()(w);
    *result = temp;
}

void cmf_save(CMF* cmf, const char* filename) {
    save_CMF_to_file(filename, *cmf);
}

void save_data(string filename, vector<Vec> &points, vector<complex<Vec>> &values, int dimension, bool with_w, bool is_complex, bool is_vector) {
    save_to_file(filename, points, values, dimension, with_w, is_complex, is_vector);
}

// Destroy CMF instance
void destroy_Field_C(Field_C* field) {
    delete field;
}

void destroy_Field_R(Field_R* field) {
    delete field;
}

void get_points(CMF* cmf, vector<Vec> &points) {
    points = cmf->points;
}
void get_w_points(CMF* cmf, vector<float> &w_points) {
    w_points = cmf->w_points;
}
void get_values(CMF* cmf, vector<complex<Vec>> &values) {
    values = cmf->values;
}
void get_domain(CMF* cmf, vector<Vec> &domain) {
    domain = cmf->domain;
}
void get_inv_domain(CMF* cmf, vector<Vec> &inv_domain) {
    inv_domain = cmf->inv_domain;
}
void get_first(CMF* cmf, Vec &first) {
    first = cmf->first;
}
void get_num_points(CMF* cmf, int &nx, int &ny, int &nz, int &nw) {
    nx = cmf->nx;
    ny = cmf->ny;
    nz = cmf->nz;
    nw = cmf->nw;
}
void get_dimension(CMF* cmf, int &dimension) {
    dimension = cmf->dimension;
}
void get_w_max_min(CMF* cmf, float &wmax, float &wmin) {
    wmax = cmf->wmax;
    wmin = cmf->wmin;
}
void get_is_complex(CMF* cmf, bool &is_complex) {
    is_complex = cmf->is_complex;
}
void get_is_vector(CMF* cmf, bool &is_vector) {
    is_vector = cmf->is_vector;
}
void get_with_w(CMF* cmf, bool &with_w) {
    with_w = cmf->with_w;
}


}
