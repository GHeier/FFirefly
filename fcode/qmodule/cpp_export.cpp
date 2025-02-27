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
Field_C* create_CMF_CS() {
    return new Field_C();
}

// Load CMF from a file
Field_C* load_cmf_cs(const char* filename) {
    return new Field_C(filename);
}

Field_R* create_CMF_RS() {
    return new Field_R();
}

Field_R* load_cmf_rs(const char* filename) {
    return new Field_R(filename);
}

// Call operator() overload for `Vec`
void cmf_cs_call(Field_C* cmf, const float *point, float w, int len, float *real_result, float *imag_result) {
    Vec q; q.dimension = len;
    for (int i = 0; i < len; i++) {
        q(i) = point[i];
    }
    complex<float> result = cmf->operator()(q, w);
    *real_result = real(result);
    *imag_result = imag(result);
}

// Call operator() overload for `float`
void cmf_cs_call2(Field_C* cmf, float w, float *real_result, float *imag_result) {
    complex<float> result = cmf->operator()(w);
    *real_result = real(result);
    *imag_result = imag(result);
}

void cmf_rs_call(Field_R* cmf, const float *point, float w, int len, float *result) {
    Vec q; q.dimension = len;
    for (int i = 0; i < len; i++) {
        q(i) = point[i];
    }
    float temp = cmf->operator()(q, w);
    *result = temp;
}

void cmf_rs_call2(Field_R* cmf, float w, float *result) {
    float temp = cmf->operator()(w);
    *result = temp;
}

void cmf_save(Field_C* cmf, const char* filename) {
    save_CMF_to_file(filename, cmf->base);
}

void save_data(string filename, vector<Vec> &points, vector<complex<Vec>> &values, int dimension, bool with_w, bool is_complex, bool is_vector) {
    save_to_file(filename, points, values, dimension, with_w, is_complex, is_vector);
}

// Destroy CMF instance
void destroy_CMF(Field_C* cmf) {
    delete cmf;
}

}
