#include <string>

#include "../../config/load/c_config.h"
#include "../../config/load/cpp_config.hpp"
#include "../../objects/surfaces.hpp"
// Begin include
#include "../../hamiltonian/band_structure.hpp"
#include "../../objects/CMField/bands.hpp"
#include "../../objects/CMField/vertex.hpp"
#include "../../objects/CMField/fields.hpp"
// End include
#include "../../objects/vec.hpp"

extern "C" void load_config_export0(const char *filename) {
    read_c_config(filename);
    load_cpp_config();
}

static float (*user_callback)(Vec);

float call_func_adapter(Vec k) { return user_callback(k); }

extern "C" Surface *Surface_export0(float (*func)(Vec), float s_val) {
    user_callback = func;
    return new Surface(call_func_adapter, s_val);
}

extern "C" int Surface_num_faces_export0(Surface *a) { return a->faces.size(); }
extern "C" void Surface_var_faces_export0(Surface *a, float *b, int *c,
                                          int *d) {
    *d = a->faces.size();
    for (int i = 0; i < *d; i++) {
        c[i] = a->faces[i].dimension;
        for (int j = 0; j < c[i]; j++) {
            b[i * c[i] + j] = a->faces[i](j);
        }
    }
}

// Begin functions
extern "C" float epsilon_export0(int a, float* b, int c) {
    Vec v(b, c);
    return epsilon(a, v);
}
extern "C" Bands* Bands_export0() {
    return new Bands();
}
extern "C" float Bands_operator_export0(Bands* a, int b, float* c, int d) {
    Vec v(c, d);
    return a->operator()(b, v);
}
extern "C" Vertex* Vertex_export0() {
    return new Vertex();
}
extern "C" void Vertex_operator_export0(Vertex* a, float* b, int c, float d, char* e, char* f, float* g, float* h) {
    Vec v(b, c);
    complex<float> r = a->operator()(v, d, e, f);
    *g = real(r);
    *h = imag(r);
}
extern "C" Field_R* Field_R_export0() {
    return new Field_R();
}
extern "C" Field_R* Field_R_export1(char* a) {
    return new Field_R(a);
}
extern "C" float Field_R_operator_export0(Field_R* a, float b) {
    return a->operator()(b);
}
extern "C" float Field_R_operator_export1(Field_R* a, int b, float c) {
    return a->operator()(b, c);
}
extern "C" float Field_R_operator_export2(Field_R* a, float* b, int c, float d) {
    Vec v(b, c);
    return a->operator()(v, d);
}
extern "C" float Field_R_operator_export3(Field_R* a, int b, float* c, int d, float e) {
    Vec v(c, d);
    return a->operator()(b, v, e);
}
extern "C" Field_C* Field_C_export0() {
    return new Field_C();
}
extern "C" Field_C* Field_C_export1(char* a) {
    return new Field_C(a);
}
extern "C" void Field_C_operator_export0(Field_C* a, float b, float* c, float* d) {
    complex<float> r = a->operator()(b);
    *c = real(r);
    *d = imag(r);
}
extern "C" void Field_C_operator_export1(Field_C* a, int b, float c, float* d, float* e) {
    complex<float> r = a->operator()(b, c);
    *d = real(r);
    *e = imag(r);
}
extern "C" void Field_C_operator_export2(Field_C* a, float* b, int c, float d, float* e, float* f) {
    Vec v(b, c);
    complex<float> r = a->operator()(v, d);
    *e = real(r);
    *f = imag(r);
}
extern "C" void Field_C_operator_export3(Field_C* a, int b, float* c, int d, float e, float* f, float* g) {
    Vec v(c, d);
    complex<float> r = a->operator()(b, v, e);
    *f = real(r);
    *g = imag(r);
}
// End functions

extern "C" {
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
