#include "../fields.hpp"
#include <filesystem>
#include <iostream>

#include "../../../config/load/c_config.h"

using namespace std;

namespace field_tests_ns {

static int mpts = 3;

static cfloat func_linear(Vec p, int dim) {
    float val = 0.0;
    for (int i = 0; i < dim; i++)
        val += p(i);
    return cfloat(val, val / 10);
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

static vector<cfloat> create_data(int dim, int pnts, vector<float> w_points = {}) {
    vector<cfloat> values;
    int w_pts = w_points.empty() ? 1 : w_points.size();

    if (dim == 3) {
        for (int i = 0; i < pnts; i++) {
            for (int j = 0; j < pnts; j++) {
                for (int k = 0; k < pnts; k++) {
                    for (int w = 0; w < w_pts; w++) {
                        float w_val = w_points.empty() ? 0.0 : w_points[w];
                        Vec point = get_vec(i, j, k, pnts);
                        cfloat base = func_linear(point, dim);
                        cfloat value = base + cfloat(w_val, w_val / 10);
                        values.push_back(value);
                    }
                }
            }
        }
    }
    if (dim == 2) {
        for (int i = 0; i < pnts; i++) {
            for (int j = 0; j < pnts; j++) {
                for (int w = 0; w < w_pts; w++) {
                    float w_val = w_points.empty() ? 0.0 : w_points[w];
                    Vec point = get_vec(i, j, 0, pnts);
                    cfloat base = func_linear(point, dim);
                    cfloat value = base + cfloat(w_val, w_val / 10);
                    values.push_back(value);
                }
            }
        }
    }
    if (dim == 1) {
        for (int i = 0; i < pnts; i++) {
            for (int w = 0; w < w_pts; w++) {
                float w_val = w_points.empty() ? 0.0 : w_points[w];
                Vec point = get_vec(i, 0, 0, pnts);
                cfloat base = func_linear(point, dim);
                cfloat value = base + cfloat(w_val, w_val / 10);
                values.push_back(value);
            }
        }
    }
    return values;
}

static bool field_r_1d_k() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<cfloat> data = create_data(1, mpts);

    Field_R field(data, mesh, domain);

    Vec v(0.0);  // Centered at origin, corresponds to x=0.5 in [0,1]
    float result = field(v);

    return fabs(result - 0.5) < 1e-6;
}

static bool field_c_1d_k() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<cfloat> data = create_data(1, mpts);

    Field_C field(data, mesh, domain);

    Vec v(0.0);
    cfloat result = field(v);

    return fabs(result - cfloat(0.5, 0.05)) < 1e-6;
}

static bool field_r_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<cfloat> data = create_data(2, mpts);

    Field_R field(data, mesh, domain);

    Vec v(0.0, 0.0);  // Center corresponds to (0.5, 0.5)
    float result = field(v);

    return fabs(result - 1.0) < 1e-6;
}

static bool field_c_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<cfloat> data = create_data(2, mpts);

    Field_C field(data, mesh, domain);

    Vec v(0.0, 0.0);
    cfloat result = field(v);

    return fabs(result - cfloat(1.0, 0.1)) < 1e-6;
}

static bool field_r_3d_k() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<cfloat> data = create_data(3, mpts);

    Field_R field(data, mesh, domain);

    Vec v(0.0, 0.0, 0.0);  // Center corresponds to (0.5, 0.5, 0.5)
    float result = field(v);

    return fabs(result - 1.5) < 1e-6;
}

static bool field_c_3d_k() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<cfloat> data = create_data(3, mpts);

    Field_C field(data, mesh, domain);

    Vec v(0.0, 0.0, 0.0);
    cfloat result = field(v);

    return fabs(result - cfloat(1.5, 0.15)) < 1e-6;
}

static bool field_r_1d_w() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(1, mpts, w_points);

    Field_R field(data, mesh, domain, w_points);

    Vec v(0.0);
    float result = field(v, 1.5);

    // At centered v=0 (original 0.5), w=1.5: val = 0.5 + 1.5 = 2.0
    return fabs(result - 2.0) < 1e-6;
}

static bool field_c_1d_w() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(1, mpts, w_points);

    Field_C field(data, mesh, domain, w_points);

    Vec v(0.0);
    cfloat result = field(v, 1.5);

    return fabs(result - cfloat(2.0, 0.20)) < 1e-6;
}

static bool field_r_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(2, mpts, w_points);

    Field_R field(data, mesh, domain, w_points);

    Vec v(0.1, 0.1);  // Centered coords
    float result = field(v, 1.1);

    // Original point at (0.6, 0.6): spatial = 1.2, w = 1.1, total = 2.3
    return fabs(result - 2.3) < 1e-6;
}

static bool field_c_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(2, mpts, w_points);

    Field_C field(data, mesh, domain, w_points);

    Vec v(0.1, 0.1);
    cfloat result = field(v, 1.1);

    return fabs(result - cfloat(2.3, 0.23)) < 1e-6;
}

static bool field_r_3d_w() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(3, mpts, w_points);

    Field_R field(data, mesh, domain, w_points);

    Vec v(-0.25, -0.25, -0.25);  // Centered coords -> original (0.25, 0.25, 0.25)
    float result = field(v, 1.5);

    // At (0.25, 0.25, 0.25), w=1.5: spatial = 0.75, w = 1.5, total = 2.25
    return fabs(result - 2.25) < 1e-6;
}

static bool field_c_3d_w() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(3, mpts, w_points);

    Field_C field(data, mesh, domain, w_points);

    Vec v(-0.25, -0.25, -0.25);
    cfloat result = field(v, 1.5);

    return fabs(result - cfloat(2.25, 0.225)) < 1e-6;
}

static bool create_destroy() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<cfloat> data = create_data(1, mpts);

    Field_C field(data, mesh, domain);

    string fname = "testfield.h5";
    field.save(fname);

    Field_C loaded(fname);

    Vec v(0.0);
    cfloat result1 = field(v);
    cfloat result2 = loaded(v);

    return fabs(result1 - result2) < 1e-6;
}

// ============= Matrix Field Tests =============

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

// Create matrix data for all spatial points
static vector<vector<vector<cfloat>>> create_matrix_data(int dim, int pnts, int mat_size, vector<float> w_points = {}) {
    vector<vector<vector<cfloat>>> data;
    int w_pts = w_points.empty() ? 1 : w_points.size();

    if (dim == 1) {
        for (int i = 0; i < pnts; i++) {
            for (int w = 0; w < w_pts; w++) {
                float w_val = w_points.empty() ? 0.0 : w_points[w];
                Vec point = get_vec(i, 0, 0, pnts);
                data.push_back(create_matrix(point, dim, mat_size, w_val));
            }
        }
    }
    if (dim == 2) {
        for (int i = 0; i < pnts; i++) {
            for (int j = 0; j < pnts; j++) {
                for (int w = 0; w < w_pts; w++) {
                    float w_val = w_points.empty() ? 0.0 : w_points[w];
                    Vec point = get_vec(i, j, 0, pnts);
                    data.push_back(create_matrix(point, dim, mat_size, w_val));
                }
            }
        }
    }
    if (dim == 3) {
        for (int i = 0; i < pnts; i++) {
            for (int j = 0; j < pnts; j++) {
                for (int k = 0; k < pnts; k++) {
                    for (int w = 0; w < w_pts; w++) {
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

static bool field_cm_1d_k() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_matrix_data(1, mpts, mat_dim);

    Field_CM field(data, mat_dim, mesh, domain);

    Vec v(0.0);  // Centered at origin, corresponds to x=0.5 in [0,1]
    auto result = field(v);

    // At x=0.5: spatial_val = 0.5, M_00 = 0.5 + 0 + 0 = 0.5
    cfloat expected(0.5, 0.05);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_1d_k() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_matrix_data(1, mpts, mat_dim);

    Field_RM field(data, mat_dim, mesh, domain);

    Vec v(0.0);  // Centered at origin, corresponds to x=0.5 in [0,1]
    auto result = field(v);

    // At x=0.5: spatial_val = 0.5, M_00 = 0.5, M_11 = 1.5, M_22 = 2.5
    float expected_00 = 0.5;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected_00) < 1e-6;
}

static bool field_cm_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    auto data = create_matrix_data(2, mpts, mat_dim);

    Field_CM field(data, mat_dim, mesh, domain);

    Vec v(0.0, 0.0);  // Center corresponds to (0.5, 0.5)
    auto result = field(v);

    // At (0.5, 0.5): spatial_val = 1.0, M_00 = 1.0
    cfloat expected(1.0, 0.1);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    auto data = create_matrix_data(2, mpts, mat_dim);

    Field_RM field(data, mat_dim, mesh, domain);

    Vec v(0.0, 0.0);  // Center corresponds to (0.5, 0.5)
    auto result = field(v);

    // At (0.5, 0.5): spatial_val = 1.0, M_00 = 1.0
    float expected_00 = 1.0;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected_00) < 1e-6;
}

static bool field_cm_1d_w() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(1, mpts, mat_dim, w_points);

    Field_CM field(data, mat_dim, mesh, domain, w_points);

    Vec v(0.0);  // Centered at x=0.5
    auto result = field(v, 1.5);

    // At x=0.5, w=1.5: spatial_val = 0.5, w_val = 1.5, M_00 = 0.5 + 0 + 0 + 1.5 = 2.0
    cfloat expected(2.0, 0.20);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_1d_w() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(1, mpts, mat_dim, w_points);

    Field_RM field(data, mat_dim, mesh, domain, w_points);

    Vec v(0.0);  // Centered at x=0.5
    auto result = field(v, 1.5);

    // At x=0.5, w=1.5: spatial_val = 0.5, w_val = 1.5, M_00 = 2.0
    float expected_00 = 2.0;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected_00) < 1e-6;
}

static bool field_cm_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(2, mpts, mat_dim, w_points);

    Field_CM field(data, mat_dim, mesh, domain, w_points);

    Vec v(0.1, 0.1);  // Centered coords -> (0.6, 0.6)
    auto result = field(v, 1.1);

    // At (0.6, 0.6), w=1.1: spatial_val = 1.2, w_val = 1.1, M_00 = 1.2 + 0 + 0 + 1.1 = 2.3
    cfloat expected(2.3, 0.23);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(2, mpts, mat_dim, w_points);

    Field_RM field(data, mat_dim, mesh, domain, w_points);

    Vec v(0.1, 0.1);  // Centered coords -> (0.6, 0.6)
    auto result = field(v, 1.1);

    // At (0.6, 0.6), w=1.1: spatial_val = 1.2, w_val = 1.1, M_00 = 2.3
    float expected_00 = 2.3;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected_00) < 1e-6;
}

static bool field_cm_3d_w() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(3, mpts, mat_dim, w_points);

    Field_CM field(data, mat_dim, mesh, domain, w_points);

    Vec v(-0.25, -0.25, -0.25);  // Centered coords -> (0.25, 0.25, 0.25)
    auto result = field(v, 1.5);

    // At (0.25, 0.25, 0.25), w=1.5: spatial_val = 0.75, w_val = 1.5, M_00 = 0.75 + 0 + 0 + 1.5 = 2.25
    cfloat expected(2.25, 0.225);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_3d_w() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(3, mpts, mat_dim, w_points);

    Field_RM field(data, mat_dim, mesh, domain, w_points);

    Vec v(-0.25, -0.25, -0.25);  // Centered coords -> (0.25, 0.25, 0.25)
    auto result = field(v, 1.5);

    // At (0.25, 0.25, 0.25), w=1.5: spatial_val = 0.75, w_val = 1.5, M_00 = 2.25
    float expected_00 = 2.25;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return fabs(result[0][0] - expected_00) < 1e-6;
}

} // namespace field_tests_ns

bool field_tests() {
    using namespace field_tests_ns;
    int num_tests = 23;
    bool all_tests[num_tests] = {
        create_destroy(),
        field_r_1d_k(),
        field_c_1d_k(),
        field_r_2d_k(),
        field_c_2d_k(),
        field_r_3d_k(),
        field_c_3d_k(),
        field_r_1d_w(),
        field_c_1d_w(),
        field_r_2d_w(),
        field_c_2d_w(),
        field_r_3d_w(),
        field_c_3d_w(),
        field_cm_1d_k(),
        field_rm_1d_k(),
        field_cm_2d_k(),
        field_rm_2d_k(),
        field_cm_1d_w(),
        field_rm_1d_w(),
        field_cm_2d_w(),
        field_rm_2d_w(),
        field_cm_3d_w(),
        field_rm_3d_w(),
    };
    filesystem::remove("testfield.h5");
    return print_test_results(all_tests, num_tests, "Field tests");
}
