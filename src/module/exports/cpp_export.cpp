#include <string>

#include "../../config/load/c_config.h"
#include "../../config/load/cpp_config.hpp"
#include "../../objects/CMField/cmfield.hpp"
#include "../../objects/CMField/fields.hpp"
#include "../../objects/CMField/vertex.hpp"
#include "../../objects/CMField/bands.hpp"
#include "../../objects/surfaces.hpp"
#include "../../hamiltonian/band_structure.hpp"
// Begin include
#include "../../objects/vec.hpp"
// End include

void vector_to_ptr(vector<float> r, float *a, int *b) {
    *b = r.size();
    for (int i = 0; i < *b; i++) {
        a[i] = r[i];
    }
}

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

typedef float (*callback_t)(Vec*);

static callback_t user_callback = nullptr;

float call_func_adapter(Vec* k) {
    return user_callback(k);
}
extern "C" Surface* Surface_export0(callback_t func, float s_val) {
    user_callback = func;
    return new Surface([](Vec k) {
        return call_func_adapter(&k);  // Pass pointer to callback
    }, s_val);
}
//extern "C" Surface* Surface_export0(float (*func)(Vec), float s_val) {
//    user_callback = func;
//    return new Surface(call_func_adapter, s_val);
//}

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
extern "C" void Surface_var_faces_export1(Surface *a, Vec *b) {
    int d = a->faces.size();
    for (int i = 0; i < d; i++) {
        b[i] = a->faces[i];
    }
}
// Begin Class gets
extern "C" float Vec_x_export0(Vec* a) {
    return a->x;
}
extern "C" float Vec_y_export0(Vec* a) {
    return a->y;
}
extern "C" float Vec_z_export0(Vec* a) {
    return a->z;
}
extern "C" float Vec_w_export0(Vec* a) {
    return a->w;
}
extern "C" float Vec_area_export0(Vec* a) {
    return a->area;
}
extern "C" int Vec_dimension_export0(Vec* a) {
    return a->dimension;
}
extern "C" int Vec_n_export0(Vec* a) {
    return a->n;
}
// End Class gets

// Begin Class functions
// End Class functions

// Begin functions
extern "C" Vec* string_to_vec_export0(Vec* a, char* b) {
    Vec* result = new Vec(string_to_vec(b));
    return result;
}
extern "C" void unpack_string_export0(Vec* a, char* b, float* c, int* d) {
    vector<float> result = unpack_string(b);
}
extern "C" char* vec_to_string_export0(Vec* a) {
    return strdup(vec_to_string(*a).c_str());
}
extern "C" Vec* Vec_export0() {
    return new Vec();
}
extern "C" Vec* Vec_export1(float a, float b, float c, float d, float e, int f, int g) {
    return new Vec(a, b, c, d, e, f, g);
}
extern "C" Vec* Vec_export2(vector<float> a) {
    return new Vec(a);
}
extern "C" Vec* Vec_export3(const float a, int b) {
    return new Vec(a, b);
}
extern "C" float operator_export0(Vec* a, int b) {
    return a->operator()(b);
}
extern "C" Vec* round_export0(Vec* a, int b) {
    Vec* result = new Vec(round(b));
    return result;
}
extern "C" float norm_export0(Vec* a) {
    return a->norm();
}
// End functions

extern "C" {

Bands *Bands_export0() { return new Bands(); }
float Bands_operator_export0(Bands *obj, int n, const float *point, int len) {
    Vec v(point, len);
    return obj->operator()(n, v);
}

Vec* vk_export0(int n, Vec* k, Bands* band) {
    Vec temp = vk(n, *k, *band);
    Vec *result = new Vec(temp);
    return result;
}

Vertex *Vertex_export0() { return new Vertex(); }

void Vertex_operator_export0(Vertex *obj, const float *point, int len, float w,
                             float *real_result, float *imag_result) {
    Vec v(point, len);
    complex<float> r = obj->operator()(v, w);
    *real_result = real(r);
    *imag_result = imag(r);
}

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
