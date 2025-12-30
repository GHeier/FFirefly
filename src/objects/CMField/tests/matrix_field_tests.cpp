#include "src/objects/CMField/fields.hpp"
#include <filesystem>
#include <iostream>

#include "src/config/load/c_config.h"

using namespace std;

namespace matrix_field_tests_ns {

static int mpts = 3;
static int mat_dim = 3;  // 3x3 matrices

// Helper to create matrix at a given spatial point
static vector<vector<cfloat>> create_matrix(Vec p, int dim, int mat_size, float w_val = 0.0) {
    vector<vector<cfloat>> matrix(mat_size, vector<cfloat>(mat_size));

    // Create a simple pattern: M_ij = (sum of spatial coords) + i + j + w
    float spatial_val = 0.0;
    for (int i = 0; i < dim; i++) {
        spatial_val += p(i);
    }

    for (int i = 0; i < mat_size; i++) {
        for (int j = 0; j < mat_size; j++) {
            float val = spatial_val + i + j + w_val;
            matrix[i][j] = cfloat(val, val / 10.0);
        }
    }
    return matrix;
}

static float get_pnt(int i, int pnts) {
    return 1.0 * i / (pnts - 1);
}

static Vec get_vec(int i, int j, int k, int pnts) {
    return Vec(
            get_pnt(i, pnts),
            get_pnt(j, pnts),
            get_pnt(k, pnts)
            );
}

// Create matrix data for all spatial points (w-k ordering)
static vector<vector<vector<cfloat>>> create_matrix_data(int dim, int pnts, int mat_size, vector<float> w_points = {}) {
    vector<vector<vector<cfloat>>> data;
    int w_pts = w_points.empty() ? 1 : w_points.size();

    // w-k ordering: frequency varies slowest, then spatial indices
    if (dim == 1) {
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
                float w_val = w_points.empty() ? 0.0 : w_points[w];
                Vec point = get_vec(i, 0, 0, pnts);
                data.push_back(create_matrix(point, dim, mat_size, w_val));
            }
        }
    }
    if (dim == 2) {
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
                for (int j = 0; j < pnts; j++) {
                    float w_val = w_points.empty() ? 0.0 : w_points[w];
                    Vec point = get_vec(i, j, 0, pnts);
                    data.push_back(create_matrix(point, dim, mat_size, w_val));
                }
            }
        }
    }
    if (dim == 3) {
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
                for (int j = 0; j < pnts; j++) {
                    for (int k = 0; k < pnts; k++) {
                        float w_val = w_points.empty() ? 0.0 : w_points[w];
                        Vec point = get_vec(i, j, k, pnts);
                        data.push_back(create_matrix(point, dim, mat_size, w_val));
                    }
                }
            }
        }
    }
    return data;
}

// Test Field_CM (complex matrix) - 1D spatial
static bool field_cm_1d_k() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_matrix_data(1, mpts, mat_dim);

    Field_CM field(data, {mat_dim, mat_dim}, mesh, domain);

    Vec v(0.0);  // Centered at origin, corresponds to x=0.5 in [0,1]
    auto result = field(v);

    // At x=0.5: spatial_val = 0.5, M_00 = 0.5 + 0 + 0 = 0.5
    cfloat expected(0.5, 0.05);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected) < 1e-6;
}

// Test Field_RM (real matrix) - 1D spatial
static bool field_rm_1d_k() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_matrix_data(1, mpts, mat_dim);

    Field_RM field(data, {mat_dim, mat_dim}, mesh, domain);

    Vec v(0.0);
    auto result = field(v);

    // At x=0.5: spatial_val = 0.5, M_00 = 0.5 + 0 + 0 = 0.5
    float expected = 0.5;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected) < 1e-6;
}

// Test Field_CM - 2D spatial
static bool field_cm_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    auto data = create_matrix_data(2, mpts, mat_dim);

    Field_CM field(data, {mat_dim, mat_dim}, mesh, domain);

    Vec v(0.0, 0.0);  // Center corresponds to (0.5, 0.5)
    auto result = field(v);

    // At (0.5, 0.5): spatial_val = 1.0, M_12 = 1.0 + 1 + 2 = 4.0
    cfloat expected(4.0, 0.4);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[1][2] - expected) < 1e-6;
}

// Test Field_RM - 2D spatial
static bool field_rm_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    auto data = create_matrix_data(2, mpts, mat_dim);

    Field_RM field(data, {mat_dim, mat_dim}, mesh, domain);

    Vec v(0.0, 0.0);
    auto result = field(v);

    // At (0.5, 0.5): spatial_val = 1.0, M_12 = 1.0 + 1 + 2 = 4.0
    float expected = 4.0;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[1][2] - expected) < 1e-6;
}

// Test Field_CM - 3D spatial
static bool field_cm_3d_k() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    auto data = create_matrix_data(3, mpts, mat_dim);

    Field_CM field(data, {mat_dim, mat_dim}, mesh, domain);

    Vec v(0.0, 0.0, 0.0);  // Center corresponds to (0.5, 0.5, 0.5)
    auto result = field(v);

    // At (0.5, 0.5, 0.5): spatial_val = 1.5, M_22 = 1.5 + 2 + 2 = 5.5
    cfloat expected(5.5, 0.55);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[2][2] - expected) < 1e-6;
}

// Test Field_RM - 3D spatial
static bool field_rm_3d_k() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    auto data = create_matrix_data(3, mpts, mat_dim);

    Field_RM field(data, {mat_dim, mat_dim}, mesh, domain);

    Vec v(0.0, 0.0, 0.0);
    auto result = field(v);

    // At (0.5, 0.5, 0.5): spatial_val = 1.5, M_22 = 1.5 + 2 + 2 = 5.5
    float expected = 5.5;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[2][2] - expected) < 1e-6;
}

// Test Field_CM with frequency dependence
static bool field_cm_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(2, mpts, mat_dim, w_points);

    Field_CM field(data, {mat_dim, mat_dim}, mesh, domain, w_points);

    Vec v(0.1, 0.1);  // Slightly offset from center
    float w = 1.5;
    auto result = field(v, w);

    // At centered (0.6, 0.6), w=1.5: spatial_val = 1.2, w_val = 1.5
    // M_01 = 1.2 + 0 + 1 + 1.5 = 3.7
    cfloat expected(3.7, 0.37);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][1] - expected) < 1e-5;  // Slightly looser tolerance due to interpolation
}

// Test Field_RM with frequency dependence
static bool field_rm_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(2, mpts, mat_dim, w_points);

    Field_RM field(data, {mat_dim, mat_dim}, mesh, domain, w_points);

    Vec v(0.1, 0.1);
    float w = 1.5;
    auto result = field(v, w);

    // At centered (0.6, 0.6), w=1.5: spatial_val = 1.2, w_val = 1.5
    // M_01 = 1.2 + 0 + 1 + 1.5 = 3.7
    float expected = 3.7;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][1] - expected) < 1e-5;
}

// Test save/load for Field_CM
static bool field_cm_save_load() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_matrix_data(1, mpts, mat_dim);

    Field_CM field1(data, {mat_dim, mat_dim}, mesh, domain);

    string fname = "test_matrix_field.h5";
    field1.save(fname);

    Field_CM field2(fname);

    Vec v(0.0);
    auto result1 = field1(v);
    auto result2 = field2(v);

    bool passed = true;
    for (int i = 0; i < mat_dim; i++) {
        for (int j = 0; j < mat_dim; j++) {
            if (fabs(result1[i][j] - result2[i][j]) > 1e-6) {
                passed = false;
            }
        }
    }

    filesystem::remove(fname);
    return passed;
}

// Test save/load for Field_RM
static bool field_rm_save_load() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_matrix_data(1, mpts, mat_dim);

    Field_RM field1(data, {mat_dim, mat_dim}, mesh, domain);

    string fname = "test_matrix_field.h5";
    field1.save(fname);

    Field_RM field2(fname);

    Vec v(0.0);
    auto result1 = field1(v);
    auto result2 = field2(v);

    bool passed = true;
    for (int i = 0; i < mat_dim; i++) {
        for (int j = 0; j < mat_dim; j++) {
            if (fabs(result1[i][j] - result2[i][j]) > 1e-6) {
                passed = false;
            }
        }
    }

    filesystem::remove(fname);
    return passed;
}

// Test copy assignment for Field_CM
static bool field_cm_copy() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_matrix_data(1, mpts, mat_dim);

    Field_CM field1(data, {mat_dim, mat_dim}, mesh, domain);
    Field_CM field2;
    field2 = field1;

    Vec v(0.0);
    auto result1 = field1(v);
    auto result2 = field2(v);

    bool passed = true;
    for (int i = 0; i < mat_dim; i++) {
        for (int j = 0; j < mat_dim; j++) {
            if (fabs(result1[i][j] - result2[i][j]) > 1e-6) {
                passed = false;
            }
        }
    }

    return passed;
}

// Test copy assignment for Field_RM
static bool field_rm_copy() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_matrix_data(1, mpts, mat_dim);

    Field_RM field1(data, {mat_dim, mat_dim}, mesh, domain);
    Field_RM field2;
    field2 = field1;

    Vec v(0.0);
    auto result1 = field1(v);
    auto result2 = field2(v);

    bool passed = true;
    for (int i = 0; i < mat_dim; i++) {
        for (int j = 0; j < mat_dim; j++) {
            if (fabs(result1[i][j] - result2[i][j]) > 1e-6) {
                passed = false;
            }
        }
    }

    return passed;
}

} // namespace matrix_field_tests_ns

bool matrix_field_tests() {
    using namespace matrix_field_tests_ns;
    int num_tests = 12;
    bool all_tests[num_tests] = {
        field_cm_1d_k(),
        field_rm_1d_k(),
        field_cm_2d_k(),
        field_rm_2d_k(),
        field_cm_3d_k(),
        field_rm_3d_k(),
        field_cm_2d_w(),
        field_rm_2d_w(),
        field_cm_save_load(),
        field_rm_save_load(),
        field_cm_copy(),
        field_rm_copy(),
    };
    return print_test_results(all_tests, num_tests, "Matrix Field tests");
}
