#include <string>

#include "../../config/load/c_config.h"
#include "../../config/load/cpp_config.hpp"
#include "../../objects/CMField/cmfield.hpp"
#include "../../objects/CMField/fields.hpp"
#include "../../objects/CMField/vertex.hpp"
#include "../../objects/CMField/self_energy.hpp"
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


void data_save_export0(string filename, const float *points, const float *values, int num_points, int dimension,
               bool with_w, bool with_n, bool is_complex, bool is_vector) {
    vector<Vec> cpoints(num_points);
    vector<complex<Vec>> cvalues(num_points);
    int a = 0;
    int c = is_complex;
    int v = is_vector;
    for (int i = 0; i < num_points; i+=dimension) {
        cpoints[i] = Vec(vector<float>(points + i, points + i + dimension));
        Vec rv(vector<float>(values + a, values + a + 1 + 3*v));
        Vec cv;
        if (is_complex) 
            cv = Vec(vector<float>(values + a, values + a + 1 + 3*v));
        complex<Vec> val = complex<Vec>(rv, cv);
        cvalues[a] = val;
    }
    CMData data(cpoints, cvalues, dimension, with_w, with_n, is_complex, is_vector);
    data.save_hdf5(filename);
}

void field_save_export0(char* filename, float *domain_c, int* mesh_c, int dimension, int nbnd, float* w_points_c, int w_size, bool is_complex, bool is_vector, bool with_w, bool with_n, float* values_c) {
    vector<vector<float>> domain(dimension);
    vector<float> first(dimension);
    int num_vals = 1;
    for (int i = 0; i < dimension; i++) {
        vector<float> temp(dimension);
        for (int j = 0; j < dimension; j++) {
            temp[j] = domain_c[i * dimension + j];
            first[j] += temp[j] * (-0.5);
        }
        domain.push_back(temp);
        num_vals *= mesh_c[i];
    }
    if (with_w) num_vals *= w_size;
    vector<int> mesh(mesh_c, mesh_c + dimension);
    vector<float> w_points(w_points_c, w_points_c + w_size);

    int size = (3 * is_vector  + (1 - is_vector) ) * (1 + is_complex);
    vector<vector<vector<float>>> values(nbnd, vector<vector<float>>(num_vals, vector<float>(size)));
    printf("nbnd: %d\n", nbnd);
    printf("num_vals: %d\n", num_vals);
    printf("size: %d\n", size);
    for (int i = 0; i < nbnd; i++) {
        for (int j = 0; j < num_vals; j++) {
            for (int k = 0; k < size; k++) {
                int idx = i * nbnd * num_vals + j * size + k;
                values[i][j][k] = values_c[idx];
            }

        }
    }
    string str_file = filename;
    save_to_hdf5(str_file, domain, first, mesh, dimension, nbnd, w_points, is_complex, is_vector, with_w, with_n, values);
}

Bands *Bands_export0() { return new Bands(); }
float Bands_operator_export0(Bands *obj, int n, const float *point, int len) {
    Vec v(point, len);
    return obj->operator()(n, v);
}

void Bands_operator_export0_numpy(Bands *obj, int n, const float *points, int num_points, int len, float *output) {
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        Vec v(point_row, len);
        output[i] = obj->operator()(n, v);
    }
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

Self_Energy *Self_Energy_export0() { return new Self_Energy(); }

void Self_Energy_operator_export0(Self_Energy *obj, const float *point, int len, float w,
                             float *real_result, float *imag_result) {
    Vec v(point, len);
    complex<float> r = obj->operator()(v, w);
    *real_result = real(r);
    *imag_result = imag(r);
}

Field_C *Field_C_export0() { return new Field_C(); }
Field_C *Field_C_export1(CMField cmf) { return new Field_C(cmf); }
Field_C *Field_C_export2(const char *filename) { 
    return new Field_C(filename); 
}

Field_R *Field_R_export0() { return new Field_R(); }
Field_R *Field_R_export1(CMField cmf) { return new Field_R(cmf); }
Field_R *Field_R_export2(const char *filename) { 
    return new Field_R(filename); }

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
void destroy_Field_C(Field_C *a) { delete a; }

void destroy_Field_R(Field_R *a) { delete a; }

void destroy_Bands(Bands *a) { delete a; }

void destroy_Vertex(Vertex *a) { delete a; }

void destroy_Self_Energy(Self_Energy *a) { delete a; }

void destroy_Surface(Surface *a) { delete a; }

void destroy_Vec(Vec *a) { delete a; }

void CMF_points_export0(CMField *cmf, vector<Vec> &points) {
    points = cmf->data.points;
}
void CMF_w_points_export0(CMField *cmf, vector<float> &w_points) {
    w_points = cmf->data.w_points;
}
void CMF_values_export0(CMField *cmf, vector<complex<Vec>> &values) {
    values = cmf->data.values;
}
void CMF_domain_export0(CMField *cmf, vector<Vec> &domain) { domain = cmf->domain; }
void CMF_inv_domain_export0(CMField *cmf, vector<Vec> &inv_domain) {
    inv_domain = cmf->inv_domain;
}
void CMF_first_export0(CMField *cmf, Vec &first) { first = cmf->first; }
void CMF_num_points_export0(CMField *cmf, int &nx, int &ny, int &nz, int &nw) {
    nx = cmf->nx;
    ny = cmf->ny;
    nz = cmf->nz;
    nw = cmf->nw;
}
void CMF_dimension_export0(CMField *cmf, int &dimension) {
    dimension = cmf->data.dimension;
}
void CMF_w_max_min_export0(CMField *cmf, float &wmax, float &wmin) {
    wmax = cmf->wmax;
    wmin = cmf->wmin;
}
void CMF_is_complex_export0(CMField *cmf, bool &is_complex) {
    is_complex = cmf->data.is_complex;
}
void CMF_is_vector_export0(CMField *cmf, bool &is_vector) {
    is_vector = cmf->data.is_vector;
}
void CMF_with_w_export0(CMField *cmf, bool &with_w) { with_w = cmf->data.with_w; }
}
