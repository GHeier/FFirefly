#include <string>

#include "../../config/load/c_config.h"
#include "../../config/load/cpp_config.hpp"
#include "../../objects/CMField/fields.hpp"
#include "../../objects/CMField/vertex.hpp"
#include "../../objects/CMField/self_energy.hpp"
#include "../../objects/CMField/hamiltonian.hpp"
#include "../../objects/CMField/base_data.hpp"
#include "../../objects/CMData/cmdata.hpp"
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

//extern "C" float epsilon_export0(int n, float *k, int size) {
//    Vec kvec;
//    for (int i = 0; i < size; i++) {
//        kvec(i) = k[i];
//    }
//    return epsilon(n, kvec);
//}

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

void Hamiltonian_operator_export0(Hamiltonian *obj, const float *point, int len, float w,
                                   float *real_result, float *imag_result, int *matrix_size) {
    Vec v(point, len);
    vector<vector<complex<float>>> H = obj->operator()(v, w);

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
    printf("cnbnd = %d\n", temp);
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

Field_RM *Field_RM_export0() {
    ensure_cpp_config_loaded();
    return new Field_RM();
}
Field_RM *Field_RM_export2(const char *filename) {
    ensure_cpp_config_loaded();
    return new Field_RM(filename);
}

void Field_RM_operator_export0(Field_RM *obj, const float *point, int len,
                               float w, float **result, int *matrix_size) {
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
            (*result)[idx] = mat[i][j];
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

Field_CM *Field_CM_export0() {
    ensure_cpp_config_loaded();
    return new Field_CM();
}
Field_CM *Field_CM_export2(const char *filename) {
    ensure_cpp_config_loaded();
    return new Field_CM(filename);
}

void Field_CM_operator_export0(Field_CM *obj, const float *point, int len,
                               float w, float **real_result, float **imag_result, int *matrix_size) {
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
            (*real_result)[idx] = real(mat[i][j]);
            (*imag_result)[idx] = imag(mat[i][j]);
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
// For scalar fields (n_indices = 0)
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
    save_data(filename, data, is_complex, mesh_vec, domain_vec, w_vec, 0, 0);
}

// For vector fields (n_indices = 0, but is_vector = true)
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
    save_data(filename, data, is_complex, mesh_vec, domain_vec, w_vec, 0, 0);
}

// For matrix fields (n_indices = 2)
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
    save_data(filename, data, is_complex, mesh_vec, domain_vec, w_vec, 2, mat_dim);
}

}
