#include <string>

#include "src/config/load/c_config.h"
#include "src/config/load/cpp_config.hpp"
#include "src/objects/CMField/fields.hpp"
#include "src/objects/CMField/vertex.hpp"
#include "src/objects/CMField/self_energy.hpp"
#include "src/objects/CMField/hamiltonian.hpp"
#include "src/objects/CMField/base_data.hpp"
#include "src/objects/CMData/cmdata.hpp"
#include "src/objects/CMField/bands.hpp"
#include "src/objects/surfaces.hpp"
#include "src/hamiltonian/models/band_structure.hpp"
// Begin include
#include "src/objects/vec.hpp"
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
    ensure_cpp_config_loaded();
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

extern "C" void Surface_faces_and_areas_export0(Surface *a, float *kpoints,
                                                  int *dims, float *areas, int *n) {
    *n = a->faces.size();
    for (int i = 0; i < *n; i++) {
        dims[i] = a->faces[i].dimension;
        areas[i] = a->faces[i].area;
        for (int j = 0; j < dims[i]; j++) {
            kpoints[i * dims[i] + j] = a->faces[i](j);
        }
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

//void field_save_export0(char* filename, float *domain_c, int* mesh_c, int dimension, int nbnd, float* w_points_c, int w_size, bool is_complex, bool is_vector, bool with_w, bool with_n, float* values_c) {
//    vector<vector<float>> domain(dimension);
//    vector<float> first(dimension);
//    int num_vals = 1;
//    for (int i = 0; i < dimension; i++) {
//        vector<float> temp(dimension);
//        for (int j = 0; j < dimension; j++) {
//            temp[j] = domain_c[i * dimension + j];
//            first[j] += temp[j] * (-0.5);
//        }
//        domain.push_back(temp);
//        num_vals *= mesh_c[i];
//    }
//    if (with_w) num_vals *= w_size;
//    vector<int> mesh(mesh_c, mesh_c + dimension);
//    vector<float> w_points(w_points_c, w_points_c + w_size);
//
//    int size = (3 * is_vector  + (1 - is_vector) ) * (1 + is_complex);
//    vector<vector<vector<float>>> values(nbnd, vector<vector<float>>(num_vals, vector<float>(size)));
//    for (int i = 0; i < nbnd; i++) {
//        for (int j = 0; j < num_vals; j++) {
//            for (int k = 0; k < size; k++) {
//                int idx = i * nbnd * num_vals + j * size + k;
//                values[i][j][k] = values_c[idx];
//            }
//
//        }
//    }
//    string str_file = filename;
//    save_to_hdf5(str_file, domain, first, mesh, dimension, nbnd, w_points, is_complex, is_vector, with_w, with_n, values);
//}

Bands *Bands_export0() {
    ensure_cpp_config_loaded();
    return new Bands();
}
float Bands_operator_export0(Bands *obj, int n, const float *point, int len) {
    Vec v(point, len);
    return obj->operator()(n, v);
}

float Bands_operator_export1(Bands *obj, const float *point, int len) {
    Vec v(point, len);
    return obj->operator()(v);
}

void Bands_operator_export0_numpy(Bands *obj, int n, const float *points, int num_points, int len, float *output) {
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        Vec v(point_row, len);
        output[i] = obj->operator()(n, v);
    }
}

void Bands_operator_export1_numpy(Bands *obj, const float *points, int num_points, int len, float *output) {
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        Vec v(point_row, len);
        output[i] = obj->operator()(v);
    }
}

Vertex *Vertex_export0() {
    ensure_cpp_config_loaded();
    return new Vertex();
}

void Vertex_operator_export0(Vertex *obj, const float *point, int len, float w,
                             float *real_result, float *imag_result) {
    Vec v(point, len);
    complex<float> r = obj->operator()(v, w);
    *real_result = real(r);
    *imag_result = imag(r);
}

Self_Energy *Self_Energy_export0() {
    ensure_cpp_config_loaded();
    return new Self_Energy();
}

void Self_Energy_operator_export0(Self_Energy *obj, const float *point, int len, float w,
                             float *real_result, float *imag_result) {
    Vec v(point, len);
    complex<float> r = obj->operator()(v, w);
    *real_result = real(r);
    *imag_result = imag(r);
}

Hamiltonian *Hamiltonian_export0() {
    ensure_cpp_config_loaded();
    return new Hamiltonian();
}

void Hamiltonian_operator_export0(Hamiltonian *obj, const float *point, int len,
                                   float *real_result, float *imag_result, int *matrix_size) {
    Vec v(point, len);
    vector<vector<complex<float>>> H = obj->operator()(v);

    if (H.empty()) {
        *matrix_size = 0;
        return;
    }

    int n = H.size();
    *matrix_size = n;

    // Flatten matrix to 1D array: real[0,0], imag[0,0], real[0,1], imag[0,1], ...
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i * n + j;
            real_result[idx] = real(H[i][j]);
            imag_result[idx] = imag(H[i][j]);
        }
    }
}

extern "C" int Field_R_nbnd_export0(Field_R* a) {
    // nbnd is always 1 now (no multi-band support)
    int temp = 1;
    return temp;
}

extern "C" int Field_C_nbnd_export0(Field_C* a) {
    // nbnd is always 1 now (no multi-band support)
    return 1;
}

Field_C *Field_C_export0() {
    ensure_cpp_config_loaded();
    return new Field_C();
}
Field_C *Field_C_export2(const char *filename) {
    ensure_cpp_config_loaded();
    return new Field_C(filename);
}

Field_R *Field_R_export0() {
    ensure_cpp_config_loaded();
    return new Field_R();
}
Field_R *Field_R_export2(const char *filename) {
    ensure_cpp_config_loaded();
    return new Field_R(filename);
}

void Field_C_operator_export0(Field_C *obj, float w, float *real_result,
                              float *imag_result) {
    complex<float> r = obj->operator()(w);
    *real_result = real(r);
    *imag_result = imag(r);
}
//void Field_C_operator_export1(Field_C *obj, int n, float w, float *real_result,
//                              float *imag_result) {
//    complex<float> r = obj->operator()(n, w);
//    *real_result = real(r);
//    *imag_result = imag(r);
//}
void Field_C_operator_export2(Field_C *obj, const float *point, int len,
                              float w, float *real_result, float *imag_result) {
    Vec v(point, len);
    complex<float> r = obj->operator()(v, w);
    *real_result = real(r);
    *imag_result = imag(r);
}
//void Field_C_operator_export3(Field_C *obj, int n, const float *point, int len,
//                              float w, float *real_result, float *imag_result) {
//    Vec v(point, len);
//    complex<float> r = obj->operator()(n, v, w);
//    *real_result = real(r);
//    *imag_result = imag(r);
//}

void Field_C_operator_export_list(Field_C *obj, const float *points, int num_points, int len,
                                  float w, float *real_output, float *imag_output) {
    vector<Vec> vec_points;
    vec_points.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        vec_points.emplace_back(point_row, len);
    }
    vector<complex<float>> results = obj->operator()(vec_points, w);
    for (int i = 0; i < num_points; ++i) {
        real_output[i] = real(results[i]);
        imag_output[i] = imag(results[i]);
    }
}

void Field_C_operator_export_w_list(Field_C *obj, const float *w_points, int num_w, float *real_output, float *imag_output) {
    vector<float> w_vec(w_points, w_points + num_w);
    vector<complex<float>> results = obj->operator()(w_vec);
    for (int i = 0; i < num_w; ++i) {
        real_output[i] = real(results[i]);
        imag_output[i] = imag(results[i]);
    }
}

float Field_R_operator_export0(Field_R *obj, float w) {
    return obj->operator()(w);
}
//float Field_R_operator_export1(Field_R *obj, int n, float w) {
//    return obj->operator()(n, w);
//}
float Field_R_operator_export2(Field_R *obj, const float *point, int len,
                               float w) {
    Vec v(point, len);
    return obj->operator()(v, w);
}
//float Field_R_operator_export3(Field_R *obj, int n, const float *point, int len,
//                               float w) {
//    Vec v(point, len);
//    return obj->operator()(n, v, w);
//}

void Field_R_operator_export_list(Field_R *obj, const float *points, int num_points, int len,
                                   float w, float *output) {
    vector<Vec> vec_points;
    vec_points.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        vec_points.emplace_back(point_row, len);
    }
    vector<float> results = obj->operator()(vec_points, w);
    for (int i = 0; i < num_points; ++i) {
        output[i] = results[i];
    }
}

void Field_R_operator_export_w_list(Field_R *obj, const float *w_points, int num_w, float *output) {
    vector<float> w_vec(w_points, w_points + num_w);
    vector<float> results = obj->operator()(w_vec);
    for (int i = 0; i < num_w; ++i) {
        output[i] = results[i];
    }
}

Field_RM *Field_RM_export0() {
    ensure_cpp_config_loaded();
    return new Field_RM();
}
Field_RM *Field_RM_export2(const char *filename) {
    ensure_cpp_config_loaded();
    return new Field_RM(filename);
}

void Field_RM_operator_export0(Field_RM *obj, const float *point, int len,
                               float w, float *result, int *matrix_size) {
    Vec v(point, len);
    vector<vector<float>> mat = obj->operator()(v, w);

    if (mat.empty()) {
        *matrix_size = 0;
        return;
    }

    int n = mat.size();
    *matrix_size = n;

    // Flatten matrix to 1D array
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i * n + j;
            result[idx] = mat[i][j];
        }
    }
}

void Field_RM_operator_export_list(Field_RM *obj, const float *points, int num_points, int len,
                                   float w, float *output, int *matrix_size) {
    vector<Vec> vec_points;
    vec_points.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        vec_points.emplace_back(point_row, len);
    }
    vector<vector<vector<float>>> results = obj->operator()(vec_points, w);

    if (results.empty() || results[0].empty()) {
        *matrix_size = 0;
        return;
    }

    int n = results[0].size();
    *matrix_size = n;

    // Flatten all matrices to output array
    for (int p = 0; p < num_points; ++p) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int idx = p * n * n + i * n + j;
                output[idx] = results[p][i][j];
            }
        }
    }
}

void Field_RM_operator_export_w(Field_RM *obj, float w, float *result, int *matrix_size) {
    vector<vector<float>> mat = obj->operator()(w);

    if (mat.empty()) {
        *matrix_size = 0;
        return;
    }

    int n = mat.size();
    *matrix_size = n;

    // Flatten matrix to 1D array
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i * n + j;
            result[idx] = mat[i][j];
        }
    }
}

void Field_RM_operator_export_w_list(Field_RM *obj, const float *w_points, int num_w, float *output, int *matrix_size) {
    vector<float> w_vec(w_points, w_points + num_w);
    vector<vector<vector<float>>> results = obj->operator()(w_vec);

    if (results.empty() || results[0].empty()) {
        *matrix_size = 0;
        return;
    }

    int n = results[0].size();
    *matrix_size = n;

    // Flatten all matrices to output array
    for (int w_idx = 0; w_idx < num_w; ++w_idx) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int idx = w_idx * n * n + i * n + j;
                output[idx] = results[w_idx][i][j];
            }
        }
    }
}

Field_CM *Field_CM_export0() {
    ensure_cpp_config_loaded();
    return new Field_CM();
}
Field_CM *Field_CM_export2(const char *filename) {
    ensure_cpp_config_loaded();
    return new Field_CM(filename);
}

void Field_CM_operator_export0(Field_CM *obj, const float *point, int len,
                               float w, float *real_result, float *imag_result, int *matrix_size) {
    Vec v(point, len);
    vector<vector<complex<float>>> mat = obj->operator()(v, w);

    if (mat.empty()) {
        *matrix_size = 0;
        return;
    }

    int n = mat.size();
    *matrix_size = n;

    // Flatten matrix to 1D arrays (real and imaginary parts)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i * n + j;
            real_result[idx] = real(mat[i][j]);
            imag_result[idx] = imag(mat[i][j]);
        }
    }
}

void Field_CM_operator_export_list(Field_CM *obj, const float *points, int num_points, int len,
                                   float w, float *real_output, float *imag_output, int *matrix_size) {
    vector<Vec> vec_points;
    vec_points.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        vec_points.emplace_back(point_row, len);
    }
    vector<vector<vector<complex<float>>>> results = obj->operator()(vec_points, w);

    if (results.empty() || results[0].empty()) {
        *matrix_size = 0;
        return;
    }

    int n = results[0].size();
    *matrix_size = n;

    // Flatten all matrices to output arrays
    for (int p = 0; p < num_points; ++p) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int idx = p * n * n + i * n + j;
                real_output[idx] = real(results[p][i][j]);
                imag_output[idx] = imag(results[p][i][j]);
            }
        }
    }
}

void Field_CM_operator_export_w(Field_CM *obj, float w, float *real_result, float *imag_result, int *matrix_size) {
    vector<vector<complex<float>>> mat = obj->operator()(w);

    if (mat.empty()) {
        *matrix_size = 0;
        return;
    }

    int n = mat.size();
    *matrix_size = n;

    // Flatten matrix to 1D arrays (real and imaginary parts)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i * n + j;
            real_result[idx] = real(mat[i][j]);
            imag_result[idx] = imag(mat[i][j]);
        }
    }
}

void Field_CM_operator_export_w_list(Field_CM *obj, const float *w_points, int num_w, float *real_output, float *imag_output, int *matrix_size) {
    vector<float> w_vec(w_points, w_points + num_w);
    vector<vector<vector<complex<float>>>> results = obj->operator()(w_vec);

    if (results.empty() || results[0].empty()) {
        *matrix_size = 0;
        return;
    }

    int n = results[0].size();
    *matrix_size = n;

    // Flatten all matrices to output arrays
    for (int w_idx = 0; w_idx < num_w; ++w_idx) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int idx = w_idx * n * n + i * n + j;
                real_output[idx] = real(results[w_idx][i][j]);
                imag_output[idx] = imag(results[w_idx][i][j]);
            }
        }
    }
}

// Field metadata exports - Field_R
int Field_R_get_mesh_size(Field_R *obj) {
    return obj->cmf.data.mesh.size();
}

void Field_R_get_mesh(Field_R *obj, int *mesh_out) {
    for (size_t i = 0; i < obj->cmf.data.mesh.size(); i++) {
        mesh_out[i] = obj->cmf.data.mesh[i];
    }
}

int Field_R_get_domain_rows(Field_R *obj) {
    return obj->cmf.data.domain.size();
}

int Field_R_get_domain_cols(Field_R *obj) {
    if (obj->cmf.data.domain.empty()) return 0;
    return obj->cmf.data.domain[0].size();
}

void Field_R_get_domain(Field_R *obj, float *domain_out) {
    int idx = 0;
    for (size_t i = 0; i < obj->cmf.data.domain.size(); i++) {
        for (size_t j = 0; j < obj->cmf.data.domain[i].size(); j++) {
            domain_out[idx++] = obj->cmf.data.domain[i][j];
        }
    }
}

int Field_R_get_w_points_size(Field_R *obj) {
    return obj->cmf.data.w_points.size();
}

void Field_R_get_w_points(Field_R *obj, float *w_points_out) {
    for (size_t i = 0; i < obj->cmf.data.w_points.size(); i++) {
        w_points_out[i] = obj->cmf.data.w_points[i];
    }
}

int Field_R_get_dimension(Field_R *obj) {
    return obj->cmf.data.dimension;
}

// Field metadata exports - Field_C
int Field_C_get_mesh_size(Field_C *obj) {
    return obj->cmf.data.mesh.size();
}

void Field_C_get_mesh(Field_C *obj, int *mesh_out) {
    for (size_t i = 0; i < obj->cmf.data.mesh.size(); i++) {
        mesh_out[i] = obj->cmf.data.mesh[i];
    }
}

int Field_C_get_domain_rows(Field_C *obj) {
    return obj->cmf.data.domain.size();
}

int Field_C_get_domain_cols(Field_C *obj) {
    if (obj->cmf.data.domain.empty()) return 0;
    return obj->cmf.data.domain[0].size();
}

void Field_C_get_domain(Field_C *obj, float *domain_out) {
    int idx = 0;
    for (size_t i = 0; i < obj->cmf.data.domain.size(); i++) {
        for (size_t j = 0; j < obj->cmf.data.domain[i].size(); j++) {
            domain_out[idx++] = obj->cmf.data.domain[i][j];
        }
    }
}

int Field_C_get_w_points_size(Field_C *obj) {
    return obj->cmf.data.w_points.size();
}

void Field_C_get_w_points(Field_C *obj, float *w_points_out) {
    for (size_t i = 0; i < obj->cmf.data.w_points.size(); i++) {
        w_points_out[i] = obj->cmf.data.w_points[i];
    }
}

int Field_C_get_dimension(Field_C *obj) {
    return obj->cmf.data.dimension;
}

// Field metadata exports - Field_RM
int Field_RM_get_mesh_size(Field_RM *obj) {
    return obj->cmf.data.mesh.size();
}

void Field_RM_get_mesh(Field_RM *obj, int *mesh_out) {
    for (size_t i = 0; i < obj->cmf.data.mesh.size(); i++) {
        mesh_out[i] = obj->cmf.data.mesh[i];
    }
}

int Field_RM_get_domain_rows(Field_RM *obj) {
    return obj->cmf.data.domain.size();
}

int Field_RM_get_domain_cols(Field_RM *obj) {
    if (obj->cmf.data.domain.empty()) return 0;
    return obj->cmf.data.domain[0].size();
}

void Field_RM_get_domain(Field_RM *obj, float *domain_out) {
    int idx = 0;
    for (size_t i = 0; i < obj->cmf.data.domain.size(); i++) {
        for (size_t j = 0; j < obj->cmf.data.domain[i].size(); j++) {
            domain_out[idx++] = obj->cmf.data.domain[i][j];
        }
    }
}

int Field_RM_get_w_points_size(Field_RM *obj) {
    return obj->cmf.data.w_points.size();
}

void Field_RM_get_w_points(Field_RM *obj, float *w_points_out) {
    for (size_t i = 0; i < obj->cmf.data.w_points.size(); i++) {
        w_points_out[i] = obj->cmf.data.w_points[i];
    }
}

int Field_RM_get_dimension(Field_RM *obj) {
    return obj->cmf.data.dimension;
}

// Field metadata exports - Field_CM
int Field_CM_get_mesh_size(Field_CM *obj) {
    return obj->cmf.data.mesh.size();
}

void Field_CM_get_mesh(Field_CM *obj, int *mesh_out) {
    for (size_t i = 0; i < obj->cmf.data.mesh.size(); i++) {
        mesh_out[i] = obj->cmf.data.mesh[i];
    }
}

int Field_CM_get_domain_rows(Field_CM *obj) {
    return obj->cmf.data.domain.size();
}

int Field_CM_get_domain_cols(Field_CM *obj) {
    if (obj->cmf.data.domain.empty()) return 0;
    return obj->cmf.data.domain[0].size();
}

void Field_CM_get_domain(Field_CM *obj, float *domain_out) {
    int idx = 0;
    for (size_t i = 0; i < obj->cmf.data.domain.size(); i++) {
        for (size_t j = 0; j < obj->cmf.data.domain[i].size(); j++) {
            domain_out[idx++] = obj->cmf.data.domain[i][j];
        }
    }
}

int Field_CM_get_w_points_size(Field_CM *obj) {
    return obj->cmf.data.w_points.size();
}

void Field_CM_get_w_points(Field_CM *obj, float *w_points_out) {
    for (size_t i = 0; i < obj->cmf.data.w_points.size(); i++) {
        w_points_out[i] = obj->cmf.data.w_points[i];
    }
}

int Field_CM_get_dimension(Field_CM *obj) {
    return obj->cmf.data.dimension;
}

// Destroy CMField instance
void destroy_Field_C(Field_C *a) { delete a; }

void destroy_Field_R(Field_R *a) { delete a; }

void destroy_Field_RM(Field_RM *a) { delete a; }

void destroy_Field_CM(Field_CM *a) { delete a; }

void destroy_Bands(Bands *a) { delete a; }

void destroy_Vertex(Vertex *a) { delete a; }

void destroy_Self_Energy(Self_Energy *a) { delete a; }

void destroy_Hamiltonian(Hamiltonian *a) { delete a; }

void destroy_Surface(Surface *a) { delete a; }

void destroy_Vec(Vec *a) { delete a; }


// Save data exports - handles BaseData::DataVariant conversion
// For scalar fields (rank = 0, inds = {})
// data_interleaved format: [real0, imag0, real1, imag1, ...]
void save_data_scalar_export0(const char *filename, const float *data_interleaved,
                               int total_size, bool is_complex,
                               const int *mesh, int mesh_size,
                               const float *domain_flat, int domain_rows, int domain_cols,
                               const float *w_points, int w_size) {
    // Convert flat arrays to C++ types
    vector<int> mesh_vec(mesh, mesh + mesh_size);
    vector<vector<float>> domain_vec(domain_rows, vector<float>(domain_cols));
    for (int i = 0; i < domain_rows; i++) {
        for (int j = 0; j < domain_cols; j++) {
            domain_vec[i][j] = domain_flat[i * domain_cols + j];
        }
    }
    vector<float> w_vec(w_points, w_points + w_size);

    // Convert interleaved data to complex vector (DataVariant type 0)
    vector<cfloat> data_vec(total_size);
    if (is_complex) {
        for (int i = 0; i < total_size; i++) {
            data_vec[i] = cfloat(data_interleaved[2*i], data_interleaved[2*i + 1]);
        }
    } else {
        for (int i = 0; i < total_size; i++) {
            data_vec[i] = cfloat(data_interleaved[i], 0.0f);
        }
    }

    // Create DataVariant and call save_data
    BaseData::DataVariant data = data_vec;
    vector<int> inds = {};  // Scalar field has empty inds
    save_data(filename, data, is_complex, mesh_vec, domain_vec, w_vec, inds);
}

// For vector fields (rank = 0, but is_vector = true, inds = {})
void save_data_vector_export0(const char *filename, const float *data_interleaved,
                               int nk, int vec_len, bool is_complex,
                               const int *mesh, int mesh_size,
                               const float *domain_flat, int domain_rows, int domain_cols,
                               const float *w_points, int w_size) {
    // Convert flat arrays to C++ types
    vector<int> mesh_vec(mesh, mesh + mesh_size);
    vector<vector<float>> domain_vec(domain_rows, vector<float>(domain_cols));
    for (int i = 0; i < domain_rows; i++) {
        for (int j = 0; j < domain_cols; j++) {
            domain_vec[i][j] = domain_flat[i * domain_cols + j];
        }
    }
    vector<float> w_vec(w_points, w_points + w_size);

    // Convert interleaved data to 2D complex vector (DataVariant type 1)
    vector<vector<cfloat>> data_vec(nk, vector<cfloat>(vec_len));
    int idx = 0;
    if (is_complex) {
        for (int i = 0; i < nk; i++) {
            for (int j = 0; j < vec_len; j++) {
                data_vec[i][j] = cfloat(data_interleaved[idx], data_interleaved[idx + 1]);
                idx += 2;
            }
        }
    } else {
        for (int i = 0; i < nk; i++) {
            for (int j = 0; j < vec_len; j++) {
                data_vec[i][j] = cfloat(data_interleaved[idx], 0.0f);
                idx++;
            }
        }
    }

    // Create DataVariant and call save_data
    BaseData::DataVariant data = data_vec;
    vector<int> inds = {};  // Vector field has empty inds (vec_len is handled separately)
    save_data(filename, data, is_complex, mesh_vec, domain_vec, w_vec, inds);
}

// For matrix fields (rank = 2, inds = {mat_dim, mat_dim})
void save_data_matrix_export0(const char *filename, const float *data_interleaved,
                               int num_matrices, int mat_dim, bool is_complex,
                               const int *mesh, int mesh_size,
                               const float *domain_flat, int domain_rows, int domain_cols,
                               const float *w_points, int w_size) {
    // Convert flat arrays to C++ types
    vector<int> mesh_vec(mesh, mesh + mesh_size);
    vector<vector<float>> domain_vec(domain_rows, vector<float>(domain_cols));
    for (int i = 0; i < domain_rows; i++) {
        for (int j = 0; j < domain_cols; j++) {
            domain_vec[i][j] = domain_flat[i * domain_cols + j];
        }
    }
    vector<float> w_vec(w_points, w_points + w_size);

    // Convert interleaved data to 3D complex vector (DataVariant type 2)
    vector<vector<vector<cfloat>>> data_vec(num_matrices,
        vector<vector<cfloat>>(mat_dim, vector<cfloat>(mat_dim)));

    int idx = 0;
    if (is_complex) {
        for (int m = 0; m < num_matrices; m++) {
            for (int i = 0; i < mat_dim; i++) {
                for (int j = 0; j < mat_dim; j++) {
                    data_vec[m][i][j] = cfloat(data_interleaved[idx], data_interleaved[idx + 1]);
                    idx += 2;
                }
            }
        }
    } else {
        for (int m = 0; m < num_matrices; m++) {
            for (int i = 0; i < mat_dim; i++) {
                for (int j = 0; j < mat_dim; j++) {
                    data_vec[m][i][j] = cfloat(data_interleaved[idx], 0.0f);
                    idx++;
                }
            }
        }
    }

    // Create DataVariant and call save_data
    BaseData::DataVariant data = data_vec;
    vector<int> inds = {mat_dim, mat_dim};  // 2D matrix
    save_data(filename, data, is_complex, mesh_vec, domain_vec, w_vec, inds);
}

// For 3D tensor fields (rank = 3, inds = {ten_dim, ten_dim, ten_dim})
void save_data_tensor3_export0(const char *filename, const float *data_interleaved,
                               int num_tensors, int ten_dim, bool is_complex,
                               const int *mesh, int mesh_size,
                               const float *domain_flat, int domain_rows, int domain_cols,
                               const float *w_points, int w_size) {
    // Convert flat arrays to C++ types
    vector<int> mesh_vec(mesh, mesh + mesh_size);
    vector<vector<float>> domain_vec(domain_rows, vector<float>(domain_cols));
    for (int i = 0; i < domain_rows; i++) {
        for (int j = 0; j < domain_cols; j++) {
            domain_vec[i][j] = domain_flat[i * domain_cols + j];
        }
    }
    vector<float> w_vec(w_points, w_points + w_size);

    // Convert interleaved data to 4D complex vector (DataVariant type 3)
    vector<vector<vector<vector<cfloat>>>> data_vec(num_tensors,
        vector<vector<vector<cfloat>>>(ten_dim,
            vector<vector<cfloat>>(ten_dim,
                vector<cfloat>(ten_dim))));

    int idx = 0;
    if (is_complex) {
        for (int t = 0; t < num_tensors; t++) {
            for (int i = 0; i < ten_dim; i++) {
                for (int j = 0; j < ten_dim; j++) {
                    for (int k = 0; k < ten_dim; k++) {
                        data_vec[t][i][j][k] = cfloat(data_interleaved[idx], data_interleaved[idx + 1]);
                        idx += 2;
                    }
                }
            }
        }
    } else {
        for (int t = 0; t < num_tensors; t++) {
            for (int i = 0; i < ten_dim; i++) {
                for (int j = 0; j < ten_dim; j++) {
                    for (int k = 0; k < ten_dim; k++) {
                        data_vec[t][i][j][k] = cfloat(data_interleaved[idx], 0.0f);
                        idx++;
                    }
                }
            }
        }
    }

    // Create DataVariant and call save_data
    BaseData::DataVariant data = data_vec;
    vector<int> inds = {ten_dim, ten_dim, ten_dim};  // 3D tensor
    save_data(filename, data, is_complex, mesh_vec, domain_vec, w_vec, inds);
}

// For 4D tensor fields (rank = 4, inds = {ten_dim, ten_dim, ten_dim, ten_dim})
void save_data_tensor4_export0(const char *filename, const float *data_interleaved,
                               int num_tensors, int ten_dim, bool is_complex,
                               const int *mesh, int mesh_size,
                               const float *domain_flat, int domain_rows, int domain_cols,
                               const float *w_points, int w_size) {
    // Convert flat arrays to C++ types
    vector<int> mesh_vec(mesh, mesh + mesh_size);
    vector<vector<float>> domain_vec(domain_rows, vector<float>(domain_cols));
    for (int i = 0; i < domain_rows; i++) {
        for (int j = 0; j < domain_cols; j++) {
            domain_vec[i][j] = domain_flat[i * domain_cols + j];
        }
    }
    vector<float> w_vec(w_points, w_points + w_size);

    // Convert interleaved data to 5D complex vector (DataVariant type 4)
    vector<vector<vector<vector<vector<cfloat>>>>> data_vec(num_tensors,
        vector<vector<vector<vector<cfloat>>>>(ten_dim,
            vector<vector<vector<cfloat>>>(ten_dim,
                vector<vector<cfloat>>(ten_dim,
                    vector<cfloat>(ten_dim)))));

    int idx = 0;
    if (is_complex) {
        for (int t = 0; t < num_tensors; t++) {
            for (int i = 0; i < ten_dim; i++) {
                for (int j = 0; j < ten_dim; j++) {
                    for (int k = 0; k < ten_dim; k++) {
                        for (int l = 0; l < ten_dim; l++) {
                            data_vec[t][i][j][k][l] = cfloat(data_interleaved[idx], data_interleaved[idx + 1]);
                            idx += 2;
                        }
                    }
                }
            }
        }
    } else {
        for (int t = 0; t < num_tensors; t++) {
            for (int i = 0; i < ten_dim; i++) {
                for (int j = 0; j < ten_dim; j++) {
                    for (int k = 0; k < ten_dim; k++) {
                        for (int l = 0; l < ten_dim; l++) {
                            data_vec[t][i][j][k][l] = cfloat(data_interleaved[idx], 0.0f);
                            idx++;
                        }
                    }
                }
            }
        }
    }

    // Create DataVariant and call save_data
    BaseData::DataVariant data = data_vec;
    vector<int> inds = {ten_dim, ten_dim, ten_dim, ten_dim};  // 4D tensor
    save_data(filename, data, is_complex, mesh_vec, domain_vec, w_vec, inds);
}

// BaseData exports
extern "C" BaseData* BaseData_load(const char *filename) {
    return new BaseData(load_data_from_hdf5(filename));
}

extern "C" BaseData* BaseData_load_with_ordering(const char *filename, const char *ordering) {
    return new BaseData(load_data_from_hdf5(filename, ordering));
}

extern "C" void BaseData_save(BaseData *data, const char *filename) {
    save_data_to_hdf5(*data, filename);
}

extern "C" void BaseData_save_with_ordering(BaseData *data, const char *filename, const char *ordering) {
    save_data_to_hdf5(*data, filename, ordering);
}

extern "C" void destroy_BaseData(BaseData *data) {
    delete data;
}

// BaseData metadata getters
extern "C" int BaseData_get_is_complex(BaseData *data) { return data->is_complex; }
extern "C" int BaseData_get_is_vector(BaseData *data) { return data->is_vector; }
extern "C" int BaseData_get_is_matrix(BaseData *data) { return data->is_matrix; }
extern "C" int BaseData_get_with_k(BaseData *data) { return data->with_k; }
extern "C" int BaseData_get_with_w(BaseData *data) { return data->with_w; }
extern "C" int BaseData_get_as_mesh(BaseData *data) { return data->as_mesh; }
extern "C" int BaseData_get_rank(BaseData *data) { return data->rank(); }
extern "C" int BaseData_get_inds_size(BaseData *data) { return data->inds.size(); }
extern "C" void BaseData_get_inds(BaseData *data, int *inds_out) {
    for (size_t i = 0; i < data->inds.size(); i++) {
        inds_out[i] = data->inds[i];
    }
}
extern "C" int BaseData_get_dimension(BaseData *data) { return data->dimension; }
extern "C" int BaseData_get_nk(BaseData *data) { return data->nk(); }
extern "C" int BaseData_get_nw(BaseData *data) { return data->nw(); }

// BaseData array getters
extern "C" int BaseData_get_mesh_size(BaseData *data) {
    return data->mesh.size();
}

extern "C" void BaseData_get_mesh(BaseData *data, int *mesh_out) {
    for (size_t i = 0; i < data->mesh.size(); i++) {
        mesh_out[i] = data->mesh[i];
    }
}

extern "C" int BaseData_get_domain_rows(BaseData *data) {
    return data->domain.size();
}

extern "C" int BaseData_get_domain_cols(BaseData *data) {
    if (data->domain.empty()) return 0;
    return data->domain[0].size();
}

extern "C" void BaseData_get_domain(BaseData *data, float *domain_out) {
    int idx = 0;
    for (size_t i = 0; i < data->domain.size(); i++) {
        for (size_t j = 0; j < data->domain[i].size(); j++) {
            domain_out[idx++] = data->domain[i][j];
        }
    }
}

extern "C" int BaseData_get_w_points_size(BaseData *data) {
    return data->w_points.size();
}

extern "C" void BaseData_get_w_points(BaseData *data, float *w_points_out) {
    for (size_t i = 0; i < data->w_points.size(); i++) {
        w_points_out[i] = data->w_points[i];
    }
}

// BaseData data extraction
extern "C" void BaseData_get_data_scalar(BaseData *data, float *real_out, float *imag_out) {
    auto& flat = data->get<std::vector<cfloat>>();
    for (size_t i = 0; i < flat.size(); i++) {
        real_out[i] = flat[i].real();
        if (data->is_complex) {
            imag_out[i] = flat[i].imag();
        }
    }
}

extern "C" void BaseData_get_data_matrix(BaseData *data, float *real_out, float *imag_out) {
    auto& matrices = data->get<std::vector<std::vector<std::vector<cfloat>>>>();
    int idx = 0;
    for (const auto& mat : matrices) {
        for (const auto& row : mat) {
            for (const auto& val : row) {
                real_out[idx] = val.real();
                if (data->is_complex) {
                    imag_out[idx] = val.imag();
                }
                idx++;
            }
        }
    }
}

extern "C" void BaseData_get_data_tensor3(BaseData *data, float *real_out, float *imag_out) {
    auto& tensors = data->get<std::vector<std::vector<std::vector<std::vector<cfloat>>>>>();
    int idx = 0;
    for (const auto& tensor : tensors) {
        for (const auto& slice : tensor) {
            for (const auto& row : slice) {
                for (const auto& val : row) {
                    real_out[idx] = val.real();
                    if (data->is_complex) {
                        imag_out[idx] = val.imag();
                    }
                    idx++;
                }
            }
        }
    }
}

extern "C" void BaseData_get_data_tensor4(BaseData *data, float *real_out, float *imag_out) {
    auto& tensors = data->get<std::vector<std::vector<std::vector<std::vector<std::vector<cfloat>>>>>>();
    int idx = 0;
    for (const auto& tensor : tensors) {
        for (const auto& cube : tensor) {
            for (const auto& slice : cube) {
                for (const auto& row : slice) {
                    for (const auto& val : row) {
                        real_out[idx] = val.real();
                        if (data->is_complex) {
                            imag_out[idx] = val.imag();
                        }
                        idx++;
                    }
                }
            }
        }
    }
}

// Field get_data exports
extern "C" BaseData* Field_R_get_data(Field_R *obj) {
    return obj->get_data();
}

extern "C" BaseData* Field_C_get_data(Field_C *obj) {
    return obj->get_data();
}

extern "C" BaseData* Field_RM_get_data(Field_RM *obj) {
    return obj->get_data();
}

extern "C" BaseData* Field_CM_get_data(Field_CM *obj) {
    return obj->get_data();
}

// Field exports
extern "C" Field* Field_create() {
    return new Field();
}

extern "C" Field* Field_from_file(const char* filename) {
    return new Field(string(filename));
}

extern "C" void Field_destroy(Field* obj) {
    delete obj;
}

extern "C" void Field_save(Field* obj, const char* filename) {
    obj->save(string(filename));
}

extern "C" bool Field_is_complex(Field* obj) {
    return obj->is_complex;
}

extern "C" bool Field_is_vector(Field* obj) {
    return obj->is_vector;
}

extern "C" bool Field_is_matrix(Field* obj) {
    return obj->is_matrix;
}

extern "C" const char* Field_get_default_plot_type(Field* obj) {
    return obj->default_plot_type.c_str();
}

extern "C" const char* Field_get_title(Field* obj) {
    return obj->title.c_str();
}

extern "C" const char* Field_get_x_label(Field* obj) {
    return obj->x_label.c_str();
}

extern "C" const char* Field_get_y_label(Field* obj) {
    return obj->y_label.c_str();
}

// Scalar complex operators
extern "C" void Field_call_scalar_complex_kw(Field* obj, const float* point, int len, float w, float* real_out, float* imag_out) {
    Vec vec;
    vec.dimension = len;
    for (int i = 0; i < len; i++) vec(i) = point[i];
    auto result = obj->operator_scalar_complex(vec, w);
    *real_out = result.real();
    *imag_out = result.imag();
}

extern "C" void Field_call_scalar_complex_w(Field* obj, float w, float* real_out, float* imag_out) {
    auto result = obj->operator_scalar_complex(w);
    *real_out = result.real();
    *imag_out = result.imag();
}

extern "C" void Field_call_scalar_complex_list(Field* obj, const float* points, int num_points, int len, float w, float* real_out, float* imag_out) {
    vector<Vec> vecs(num_points);
    for (int i = 0; i < num_points; i++) {
        vecs[i].dimension = len;
        for (int j = 0; j < len; j++) {
            vecs[i](j) = points[i * len + j];
        }
    }
    auto results = obj->operator_scalar_complex(vecs, w);
    for (size_t i = 0; i < results.size(); i++) {
        real_out[i] = results[i].real();
        imag_out[i] = results[i].imag();
    }
}

// Scalar real operators
extern "C" float Field_call_scalar_real_kw(Field* obj, const float* point, int len, float w) {
    Vec vec;
    vec.dimension = len;
    for (int i = 0; i < len; i++) vec(i) = point[i];
    return obj->operator_scalar_real(vec, w);
}

extern "C" float Field_call_scalar_real_w(Field* obj, float w) {
    return obj->operator_scalar_real(w);
}

extern "C" void Field_call_scalar_real_list(Field* obj, const float* points, int num_points, int len, float w, float* out) {
    vector<Vec> vecs(num_points);
    for (int i = 0; i < num_points; i++) {
        vecs[i].dimension = len;
        for (int j = 0; j < len; j++) {
            vecs[i](j) = points[i * len + j];
        }
    }
    auto results = obj->operator_scalar_real(vecs, w);
    for (size_t i = 0; i < results.size(); i++) {
        out[i] = results[i];
    }
}

// Matrix complex operators
extern "C" void Field_call_matrix_complex(Field* obj, const float* point, int len, float w, float* real_out, float* imag_out, int* size_out) {
    Vec vec;
    vec.dimension = len;
    for (int i = 0; i < len; i++) vec(i) = point[i];
    auto result = obj->operator_matrix_complex(vec, w);
    *size_out = result.size();
    for (size_t i = 0; i < result.size(); i++) {
        for (size_t j = 0; j < result[i].size(); j++) {
            int idx = i * result[i].size() + j;
            real_out[idx] = result[i][j].real();
            imag_out[idx] = result[i][j].imag();
        }
    }
}

// Matrix real operators
extern "C" void Field_call_matrix_real(Field* obj, const float* point, int len, float w, float* out, int* size_out) {
    Vec vec;
    vec.dimension = len;
    for (int i = 0; i < len; i++) vec(i) = point[i];
    auto result = obj->operator_matrix_real(vec, w);
    *size_out = result.size();
    for (size_t i = 0; i < result.size(); i++) {
        for (size_t j = 0; j < result[i].size(); j++) {
            out[i * result[i].size() + j] = result[i][j];
        }
    }
}

extern "C" BaseData* Field_get_data(Field *obj) {
    return obj->get_data();
}

// ============================================================================
// Additional Hamiltonian exports
// ============================================================================

void Hamiltonian_operator_export_list(Hamiltonian *obj, const float *points, int num_points, int len, float *real_output, float *imag_output, int *matrix_size) {
    vector<Vec> vec_points;
    vec_points.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        vec_points.emplace_back(point_row, len);
    }

    // For Hamiltonian, we need to call operator() for each point individually
    // since it doesn't have a batch interface like Field_CM
    if (vec_points.empty()) {
        *matrix_size = 0;
        return;
    }

    auto first_result = obj->operator()(vec_points[0]);
    if (first_result.empty()) {
        *matrix_size = 0;
        return;
    }

    int n = first_result.size();
    *matrix_size = n;

    // Store first result
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int idx = i * n + j;
            real_output[idx] = real(first_result[i][j]);
            imag_output[idx] = imag(first_result[i][j]);
        }
    }

    // Process remaining points
    for (int k = 1; k < num_points; ++k) {
        auto mat = obj->operator()(vec_points[k]);
        int offset = k * n * n;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int idx = offset + i * n + j;
                real_output[idx] = real(mat[i][j]);
                imag_output[idx] = imag(mat[i][j]);
            }
        }
    }
}

bool Hamiltonian_file_found(Hamiltonian *obj) {
    return obj->file_found;
}

void Hamiltonian_get_bands_export0(Hamiltonian *obj, const float *point, int len,
                                    float *eigenvalues_out, int *num_bands) {
    Vec v(point, len);
    vector<float> eigs = obj->get_bands(v);

    *num_bands = eigs.size();
    for (size_t i = 0; i < eigs.size(); i++) {
        eigenvalues_out[i] = eigs[i];
    }
}

void Hamiltonian_get_wavefunctions_export0(Hamiltonian *obj, const float *point, int len,
                                            float *eigenvalues_out, float *eigvecs_real,
                                            float *eigvecs_imag, int *num_bands) {
    Vec v(point, len);
    vector<eigvec> result = obj->get_wavefunctions(v);

    if (result.empty()) {
        *num_bands = 0;
        return;
    }

    int n = result.size();
    *num_bands = n;

    // Extract eigenvalues and eigenvectors
    for (int i = 0; i < n; i++) {
        eigenvalues_out[i] = result[i].eigenvalue;

        // Eigenvectors: store i-th eigenvector in i-th column
        for (int j = 0; j < n; j++) {
            int idx = j * n + i;  // Column-major: eigvec i is in column i
            eigvecs_real[idx] = result[i].eigenvector[j].real();
            eigvecs_imag[idx] = result[i].eigenvector[j].imag();
        }
    }
}

void Hamiltonian_get_bands_export_list(Hamiltonian *obj, const float *points, int num_points, int len,
                                        float *eigenvalues_out, int *num_bands) {
    vector<Vec> vec_points;
    vec_points.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        vec_points.emplace_back(point_row, len);
    }

    vector<vector<float>> results = obj->get_bands(vec_points);

    if (results.empty() || results[0].empty()) {
        *num_bands = 0;
        return;
    }

    int n = results[0].size();
    *num_bands = n;

    // Flatten all eigenvalue arrays to output
    for (int p = 0; p < num_points; ++p) {
        for (int i = 0; i < n; i++) {
            eigenvalues_out[p * n + i] = results[p][i];
        }
    }
}

void Hamiltonian_get_wavefunctions_export_list(Hamiltonian *obj, const float *points, int num_points, int len,
                                                 float *eigenvalues_out, float *eigvecs_real,
                                                 float *eigvecs_imag, int *num_bands) {
    vector<Vec> vec_points;
    vec_points.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        vec_points.emplace_back(point_row, len);
    }

    vector<vector<eigvec>> results = obj->get_wavefunctions(vec_points);

    if (results.empty() || results[0].empty()) {
        *num_bands = 0;
        return;
    }

    int n = results[0].size();
    *num_bands = n;

    // Extract eigenvalues and eigenvectors for all k-points
    for (int p = 0; p < num_points; ++p) {
        for (int i = 0; i < n; i++) {
            // Store eigenvalue
            eigenvalues_out[p * n + i] = results[p][i].eigenvalue;

            // Store eigenvectors: for k-point p, eigenvector i, component j
            for (int j = 0; j < n; j++) {
                int idx = p * n * n + j * n + i;  // Column-major: eigvec i is in column i
                eigvecs_real[idx] = results[p][i].eigenvector[j].real();
                eigvecs_imag[idx] = results[p][i].eigenvector[j].imag();
            }
        }
    }
}

extern "C" void Hamiltonian_get_fermi_velocity_export0(Hamiltonian *obj, const float *point, int len,
                                                        float *velocities_out, int *num_bands) {
    Vec v(point, len);
    vector<Vec> vels = obj->get_fermi_velocity(v);

    *num_bands = vels.size();
    for (size_t i = 0; i < vels.size(); i++) {
        velocities_out[i * 3 + 0] = vels[i].x;
        velocities_out[i * 3 + 1] = vels[i].y;
        velocities_out[i * 3 + 2] = vels[i].z;
    }
}

extern "C" void Hamiltonian_get_fermi_velocity_export_list(Hamiltonian *obj, const float *points,
                                                             int num_points, int len,
                                                             float *velocities_out, int *num_bands) {
    vector<Vec> vec_points;
    vec_points.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        const float* point_row = points + i * len;
        vec_points.emplace_back(point_row, len);
    }

    vector<vector<Vec>> results = obj->get_fermi_velocity(vec_points);

    if (results.empty() || results[0].empty()) {
        *num_bands = 0;
        return;
    }

    int n = results[0].size();
    *num_bands = n;

    // Flatten all velocity vectors to output
    for (int p = 0; p < num_points; ++p) {
        for (int i = 0; i < n; i++) {
            int idx = (p * n + i) * 3;
            velocities_out[idx + 0] = results[p][i].x;
            velocities_out[idx + 1] = results[p][i].y;
            velocities_out[idx + 2] = results[p][i].z;
        }
    }
}

}
