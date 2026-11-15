#include "src/objects/CMField/fields.hpp"
#include "src/objects/CMField/field.hpp"
#include <filesystem>
#include <H5Cpp.h>
#include <iostream>
#include <cmath>

#include "src/config/load/c_config.h"

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
    return 1.0 * i / pnts;  // periodic grid: sample i at i/n
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

    // w-k ordering: frequency varies slowest, then spatial indices
    if (dim == 3) {
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
                for (int j = 0; j < pnts; j++) {
                    for (int k = 0; k < pnts; k++) {
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
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
                for (int j = 0; j < pnts; j++) {
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
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
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

    return abs(result - 0.5) < 1e-6;
}

static bool field_c_1d_k() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<cfloat> data = create_data(1, mpts);

    Field_C field(data, mesh, domain);

    Vec v(0.0);
    cfloat result = field(v);

    return abs(result - cfloat(0.5, 0.05)) < 1e-6;
}

static bool field_r_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<cfloat> data = create_data(2, mpts);

    Field_R field(data, mesh, domain);

    Vec v(0.0, 0.0);  // Center corresponds to (0.5, 0.5)
    float result = field(v);

    return abs(result - 1.0) < 1e-6;
}

static bool field_c_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<cfloat> data = create_data(2, mpts);

    Field_C field(data, mesh, domain);

    Vec v(0.0, 0.0);
    cfloat result = field(v);

    return abs(result - cfloat(1.0, 0.1)) < 1e-6;
}

static bool field_r_3d_k() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<cfloat> data = create_data(3, mpts);

    Field_R field(data, mesh, domain);

    Vec v(0.0, 0.0, 0.0);  // Center corresponds to (0.5, 0.5, 0.5)
    float result = field(v);

    return abs(result - 1.5) < 1e-6;
}

static bool field_c_3d_k() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<cfloat> data = create_data(3, mpts);

    Field_C field(data, mesh, domain);

    Vec v(0.0, 0.0, 0.0);
    cfloat result = field(v);

    return abs(result - cfloat(1.5, 0.15)) < 1e-6;
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
    return abs(result - 2.0) < 1e-6;
}

static bool field_c_1d_w() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(1, mpts, w_points);

    Field_C field(data, mesh, domain, w_points, false);

    Vec v(0.1);
    cfloat result = field(v, 1.5);

    return abs(result - cfloat(1.6, 0.16)) < 1e-6;
}

static bool field_r_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(2, mpts, w_points);
    Field_R field(data, mesh, domain, w_points, false);

    Vec v1(0.1, 0.3);  // Centered coords
    float result = field(v1, 1.5);

    // Original point at (0.6, 0.6): spatial = 1.2, w = 1.1, total = 2.3
    return abs(result - 1.9) < 1e-6;
}

static bool field_c_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(2, mpts, w_points);

    Field_C field(data, mesh, domain, w_points);

    Vec v(0.1, 0.1);
    cfloat result = field(v, 1.1);

    return abs(result - cfloat(2.3, 0.23)) < 1e-6;
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
    return abs(result - 2.25) < 1e-6;
}

static bool field_c_3d_w() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    vector<cfloat> data = create_data(3, mpts, w_points);

    Field_C field(data, mesh, domain, w_points);

    Vec v(-0.25, -0.25, -0.25);
    cfloat result = field(v, 1.5);

    return abs(result - cfloat(2.25, 0.225)) < 1e-6;
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

    return abs(result1 - result2) < 1e-6;
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

    return abs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_1d_k() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_matrix_data(1, mpts, mat_dim);

    Field_RM field(data, {mat_dim, mat_dim}, mesh, domain);

    Vec v(0.0);  // Centered at origin, corresponds to x=0.5 in [0,1]
    auto result = field(v);

    // At x=0.5: spatial_val = 0.5, M_00 = 0.5, M_11 = 1.5, M_22 = 2.5
    float expected_00 = 0.5;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return abs(result[0][0] - expected_00) < 1e-6;
}

static bool field_cm_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    auto data = create_matrix_data(2, mpts, mat_dim);

    Field_CM field(data, {mat_dim, mat_dim}, mesh, domain);

    Vec v(0.0, 0.0);  // Center corresponds to (0.5, 0.5)
    auto result = field(v);

    // At (0.5, 0.5): spatial_val = 1.0, M_00 = 1.0
    cfloat expected(1.0, 0.1);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return abs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_2d_k() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    auto data = create_matrix_data(2, mpts, mat_dim);

    Field_RM field(data, {mat_dim, mat_dim}, mesh, domain);

    Vec v(0.0, 0.0);  // Center corresponds to (0.5, 0.5)
    auto result = field(v);

    // At (0.5, 0.5): spatial_val = 1.0, M_00 = 1.0
    float expected_00 = 1.0;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return abs(result[0][0] - expected_00) < 1e-6;
}

static bool field_cm_1d_w() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(1, mpts, mat_dim, w_points);

    Field_CM field(data, {mat_dim, mat_dim}, mesh, domain, w_points);

    Vec v(0.0);  // Centered at x=0.5
    auto result = field(v, 1.5);

    // At x=0.5, w=1.5: spatial_val = 0.5, w_val = 1.5, M_00 = 0.5 + 0 + 0 + 1.5 = 2.0
    cfloat expected(2.0, 0.20);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return abs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_1d_w() {
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(1, mpts, mat_dim, w_points);

    Field_RM field(data, {mat_dim, mat_dim}, mesh, domain, w_points);

    Vec v(0.0);  // Centered at x=0.5
    auto result = field(v, 1.5);

    // At x=0.5, w=1.5: spatial_val = 0.5, w_val = 1.5, M_00 = 2.0
    float expected_00 = 2.0;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return abs(result[0][0] - expected_00) < 1e-6;
}

static bool field_cm_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(2, mpts, mat_dim, w_points);

    Field_CM field(data, {mat_dim, mat_dim}, mesh, domain, w_points);

    Vec v(0.1, 0.1);  // Centered coords -> (0.6, 0.6)
    auto result = field(v, 1.1);

    // At (0.6, 0.6), w=1.1: spatial_val = 1.2, w_val = 1.1, M_00 = 1.2 + 0 + 0 + 1.1 = 2.3
    cfloat expected(2.3, 0.23);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return abs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_2d_w() {
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(2, mpts, mat_dim, w_points);

    Field_RM field(data, {mat_dim, mat_dim}, mesh, domain, w_points);

    Vec v(0.1, 0.1);  // Centered coords -> (0.6, 0.6)
    auto result = field(v, 1.1);

    // At (0.6, 0.6), w=1.1: spatial_val = 1.2, w_val = 1.1, M_00 = 2.3
    float expected_00 = 2.3;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return abs(result[0][0] - expected_00) < 1e-6;
}

static bool field_cm_3d_w() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(3, mpts, mat_dim, w_points);

    Field_CM field(data, {mat_dim, mat_dim}, mesh, domain, w_points);

    Vec v(-0.25, -0.25, -0.25);  // Centered coords -> (0.25, 0.25, 0.25)
    auto result = field(v, 1.5);

    // At (0.25, 0.25, 0.25), w=1.5: spatial_val = 0.75, w_val = 1.5, M_00 = 0.75 + 0 + 0 + 1.5 = 2.25
    cfloat expected(2.25, 0.225);

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return abs(result[0][0] - expected) < 1e-6;
}

static bool field_rm_3d_w() {
    vector<int> mesh = {mpts, mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    vector<float> w_points = {1.0, 2.0, 3.0};
    auto data = create_matrix_data(3, mpts, mat_dim, w_points);

    Field_RM field(data, {mat_dim, mat_dim}, mesh, domain, w_points);

    Vec v(-0.25, -0.25, -0.25);  // Centered coords -> (0.25, 0.25, 0.25)
    auto result = field(v, 1.5);

    // At (0.25, 0.25, 0.25), w=1.5: spatial_val = 0.75, w_val = 1.5, M_00 = 2.25
    float expected_00 = 2.25;

    if (result.size() != mat_dim || result[0].size() != mat_dim) {
        return false;
    }

    return abs(result[0][0] - expected_00) < 1e-6;
}


// ============= Regression Tests =============

// Non-cubic 3D mesh: flat index must use (i*ny + j)*nz + k
static bool regression_noncubic_3d() {
    int nx = 4, ny = 3, nz = 2;
    vector<cfloat> data;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                data.push_back(cfloat(1.0f * i / nx + 10.0f * j / ny, 0));
    vector<vector<float>> domain = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    Field_R field(data, {nx, ny, nz}, domain, {}, false);
    float result = field(Vec(0.5, 0.2, 0.7));
    return abs(result - 2.5f) < 1e-5;
}

// Periodic grid: between the last sample and the next period's first
static bool regression_periodic_wrap() {
    int n = 4;
    vector<cfloat> data;
    for (int i = 0; i < n; i++) data.push_back(cfloat(cos(2 * M_PI * i / n), 0));
    Field_R field(data, {n}, {{1.0}}, {}, false);
    // x = 0.875 lies halfway between x=3/4 (f=0) and x=1~0 (f=1)
    bool wrap = abs(field(Vec(0.875)) - 0.5f) < 1e-5;
    // Grid points are reproduced exactly, including outside the first zone
    bool exact = abs(field(Vec(0.25)) - 0.0f) < 1e-5 && abs(field(Vec(-1.0)) - 1.0f) < 1e-5;
    return wrap && exact;
}

// Non-orthogonal (hexagonal) basis: fractional coords are (B^T)^-1 k
static bool regression_hex_domain() {
    int n = 4;
    vector<vector<float>> B = {{1.0f, 0.0f}, {-0.5f, (float)(sqrt(3.0) / 2)}};
    vector<cfloat> data;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            data.push_back(cfloat(1.0f * i / n, 1.0f * j / n));
    Field_C uncentered(data, {n, n}, B, {}, false);
    Field_C centered(data, {n, n}, B, {}, true);
    Vec f(0.3, 0.6);
    f.dimension = 2;
    Vec k = frac_to_cart(B, f);
    // Centered: sample i sits at fractional coordinate i/n - 1/2
    Vec fc(0.3 - 0.5, 0.6 - 0.5);
    fc.dimension = 2;
    Vec kc = frac_to_cart(B, fc);
    return abs(uncentered(k) - cfloat(0.3, 0.6)) < 1e-5 && abs(centered(kc) - cfloat(0.3, 0.6)) < 1e-5;
}

// Matrix field with a single frequency point used to index w_idx = -1
static bool regression_single_w_matrix() {
    vector<vector<vector<cfloat>>> data(4, vector<vector<cfloat>>(2, vector<cfloat>(2, cfloat(1.5, 0))));
    Field_CM field(data, {2, 2}, {4}, {{1.0}}, {0.0});
    auto m = field(Vec(0.1), 0.0);
    return m.size() == 2 && abs(m[1][1] - cfloat(1.5, 0)) < 1e-6;
}

// Frequency-only real matrix field: used to crash, and Field_RM(w) always
// returned an empty matrix
static bool regression_w_only_real_matrix() {
    vector<float> w_points = {-1.0, -0.5, 0.0, 0.5, 1.0};
    vector<vector<vector<cfloat>>> data;
    for (float w : w_points) data.push_back(vector<vector<cfloat>>(2, vector<cfloat>(2, cfloat(w, 0))));
    Field_RM field(data, {2, 2}, {}, {}, w_points);
    auto m = field(0.2f);
    auto mk = field(Vec(0.3), 0.2f);
    return m.size() == 2 && abs(m[0][1] - 0.2f) < 1e-6 && mk.size() == 2 && abs(mk[1][0] - 0.2f) < 1e-6;
}

// Frequency-only evaluation of a k-dependent field averages over k
static bool regression_k_average() {
    int n = 4;
    vector<float> w_points = {0.0, 1.0};
    vector<cfloat> data;
    for (float w : w_points)
        for (int i = 0; i < n; i++) data.push_back(cfloat(i + w, 0));
    Field_C field(data, {n}, {{1.0}}, w_points);
    // mean_k(i) = 1.5, plus w = 0.5
    return abs(field(0.5f) - cfloat(2.0, 0)) < 1e-5;
}

// Rank-1 (band vector) fields round-trip through HDF5 and load as Field_CM
static bool regression_rank1_save_load() {
    int n = 4;
    vector<vector<cfloat>> data;
    for (int i = 0; i < n; i++) data.push_back({cfloat(i, 0), cfloat(10 + i, 0)});
    Field_CM field(data, {2}, {n}, {{1.0}}, {}, false);
    field.save("testfield_rank1.h5");
    Field_CM loaded("testfield_rank1.h5");
    filesystem::remove("testfield_rank1.h5");
    auto m = loaded(Vec(0.5));
    return m.size() == 2 && m[0].size() == 1 && abs(m[0][0] - cfloat(2, 0)) < 1e-5
        && abs(m[1][0] - cfloat(12, 0)) < 1e-5;
}

// save_data's default domain {{}} used to be written as a 1x0 dataset that
// crashed on reload
static bool regression_default_domain_roundtrip() {
    int n = 4;
    vector<cfloat> data;
    for (int i = 0; i < n; i++) data.push_back(cfloat(1.0f * i / n, 0));
    save_data("testfield_domain.h5", data, {}, {n});
    Field_C loaded("testfield_domain.h5");
    filesystem::remove("testfield_domain.h5");
    // centered by default: physical x = 0 is grid x = 1/2
    return abs(loaded(Vec(0.0)) - cfloat(0.5, 0)) < 1e-5;
}

// Indexed evaluation returns the requested tensor element
static bool regression_indexed_element() {
    vector<vector<vector<cfloat>>> data(3, {{cfloat(1, 0), cfloat(2, 0)}, {cfloat(3, 0), cfloat(4, 0)}});
    FieldImpl field(data, true, false, {3}, {{1.0}}, {}, {2, 2});
    auto r = field(Vec(0.0), {1, 0});
    auto* v = std::get_if<vector<cfloat>>(&r);
    return v && v->size() == 1 && abs((*v)[0] - cfloat(3, 0)) < 1e-6;
}

// Vec constructed from fewer than 4 values zero-fills the rest
static bool regression_vec_init() {
    Vec a(vector<float>{1.0f});
    float pts[2] = {1.0f, 2.0f};
    Vec b(pts, 2);
    return a.y == 0 && a.z == 0 && a.w == 0 && b.z == 0 && b.w == 0;
}


// A field loaded from HDF5 keeps its values when saved again
static bool regression_load_then_save() {
    int n = 4;
    vector<cfloat> data;
    for (int i = 0; i < n; i++) data.push_back(cfloat(i, -i));
    save_data("testfield_resave_a.h5", data, {}, {n}, {{1.0}}, {}, ep(), false);
    Field_C loaded("testfield_resave_a.h5");
    loaded.save("testfield_resave_b.h5");
    Field_C again("testfield_resave_b.h5");
    filesystem::remove("testfield_resave_a.h5");
    filesystem::remove("testfield_resave_b.h5");
    return abs(again(Vec(0.5)) - cfloat(2, -2)) < 1e-5;
}

// Vector fields read from HDF5 keep their components per grid point
static bool regression_vector_field_from_file() {
    int n = 4;
    BaseData f;
    f.is_complex = true;
    f.is_vector = true;
    f.with_k = f.as_mesh = true;
    f.centered = false;
    f.mesh = {n, n};
    f.dimension = 2;
    vector<vector<cfloat>> data;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) data.push_back({cfloat(i, 0), cfloat(j, 1)});
    f.data = data;
    save_data_to_hdf5(f, "testfield_vector.h5");
    FieldImpl loaded("testfield_vector.h5");
    filesystem::remove("testfield_vector.h5");
    auto r = loaded(Vec(0.25, 0.5));
    auto* v = std::get_if<complex<Vec>>(&r);
    return v && abs(v->real().x - 1) < 1e-5 && abs(v->real().y - 2) < 1e-5
        && abs(v->imag().y - 1) < 1e-5;
}

// Fields built in code expose their values through get_data()
static bool regression_programmatic_get_data() {
    vector<cfloat> data = {cfloat(1, 0), cfloat(2, 0), cfloat(3, 0)};
    Field_C field(data, {3}, {{1.0}});
    BaseData* d = field.get_data();
    return d->real_values.size() == 3 && abs(d->real_values[2] - 3) < 1e-6;
}

// Attributes and the format version round-trip; files without a version read as 1
static bool regression_attributes() {
    BaseData f;
    f.is_complex = false;
    f.real_values = {1, 2, 3};
    f.attrs = {{"quantity", "Sigma"}, {"units", "eV"}};
    save_data_to_hdf5(f, "testfield_attrs.h5");
    write_field_attributes("testfield_attrs.h5", {{"temperature", "0.01"}});
    BaseData loaded = load_data_from_hdf5("testfield_attrs.h5");
    auto attrs = read_field_attributes("testfield_attrs.h5");
    filesystem::remove("testfield_attrs.h5");
    bool ok = loaded.format_version == BaseData::FORMAT_VERSION
        && loaded.attr("quantity") == "Sigma" && loaded.attr("units") == "eV"
        && loaded.attr("temperature") == "0.01" && attrs.size() == 3
        && loaded.attr("missing", "none") == "none";

    // Legacy file: same layout, no format_version attribute
    {
        H5::H5File file("testfield_legacy.h5", H5F_ACC_TRUNC);
        auto scalar = [&](const char* name, int v) {
            H5::DataSet ds = file.createDataSet(name, H5::PredType::NATIVE_INT, H5::DataSpace(H5S_SCALAR));
            ds.write(&v, H5::PredType::NATIVE_INT);
        };
        scalar("/is_complex", 0); scalar("/is_vector", 0); scalar("/with_k", 0);
        scalar("/with_w", 0); scalar("/dimension", 0); scalar("/as_mesh", 0);
        hsize_t zero[1] = {0}, three[1] = {3};
        file.createDataSet("/inds", H5::PredType::NATIVE_INT, H5::DataSpace(1, zero));
        H5::Group g = file.createGroup("/values");
        float vals[3] = {1, 2, 3};
        g.createDataSet("real", H5::PredType::NATIVE_FLOAT, H5::DataSpace(1, three))
            .write(vals, H5::PredType::NATIVE_FLOAT);
    }
    BaseData legacy = load_data_from_hdf5("testfield_legacy.h5");
    filesystem::remove("testfield_legacy.h5");
    return ok && legacy.format_version == 1 && legacy.attrs.empty();
}


// Plot labels use the stored quantity, units and frequency axis
static bool regression_attribute_labels() {
    vector<cfloat> data = {cfloat(1, 0), cfloat(2, 0)};
    save_data("testfield_sigma_w.h5", data, {}, {}, {{}}, {0.0f, 1.0f});
    write_field_attributes("testfield_sigma_w.h5",
                           {{"quantity", "\\Sigma"}, {"units", "eV"}, {"frequency_axis", "matsubara"}});
    Field field("testfield_sigma_w.h5");
    filesystem::remove("testfield_sigma_w.h5");
    return field.y_label == "$\\Sigma(\\omega)$ [eV]" && field.x_label == "$i\\omega_n$";
}

} // namespace field_tests_ns

bool field_tests() {
    using namespace field_tests_ns;
    int num_tests = 38;
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
        regression_noncubic_3d(),
        regression_periodic_wrap(),
        regression_hex_domain(),
        regression_single_w_matrix(),
        regression_w_only_real_matrix(),
        regression_k_average(),
        regression_rank1_save_load(),
        regression_default_domain_roundtrip(),
        regression_indexed_element(),
        regression_vec_init(),
        regression_load_then_save(),
        regression_vector_field_from_file(),
        regression_programmatic_get_data(),
        regression_attributes(),
        regression_attribute_labels(),
    };
    filesystem::remove("testfield.h5");
    return print_test_results(all_tests, num_tests, "Field tests");
}
