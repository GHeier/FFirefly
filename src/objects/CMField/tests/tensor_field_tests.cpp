#include "src/objects/CMField/fields.hpp"
#include <filesystem>
#include <iostream>

#include "src/config/load/c_config.h"

using namespace std;

namespace tensor_field_tests_ns {

static int mpts = 3;
static int ten3_dim = 2;  // 2x2x2 tensors (3 indices)
static int ten4_dim = 2;  // 2x2x2x2 tensors (4 indices)

// Helper to create 3D tensor at a given spatial point
static vector<vector<vector<cfloat>>> create_tensor3(Vec p, int dim, int ten_size, float w_val = 0.0) {
    vector<vector<vector<cfloat>>> tensor(ten_size, vector<vector<cfloat>>(ten_size, vector<cfloat>(ten_size)));

    float spatial_val = 0.0;
    for (int i = 0; i < dim; i++) {
        spatial_val += p(i);
    }

    for (int i = 0; i < ten_size; i++) {
        for (int j = 0; j < ten_size; j++) {
            for (int k = 0; k < ten_size; k++) {
                float val = spatial_val + i + j + k + w_val;
                tensor[i][j][k] = cfloat(val, val / 10.0);
            }
        }
    }
    return tensor;
}

// Helper to create 4D tensor at a given spatial point
static vector<vector<vector<vector<cfloat>>>> create_tensor4(Vec p, int dim, int ten_size, float w_val = 0.0) {
    vector<vector<vector<vector<cfloat>>>> tensor(ten_size,
        vector<vector<vector<cfloat>>>(ten_size,
            vector<vector<cfloat>>(ten_size,
                vector<cfloat>(ten_size))));

    float spatial_val = 0.0;
    for (int i = 0; i < dim; i++) {
        spatial_val += p(i);
    }

    for (int i = 0; i < ten_size; i++) {
        for (int j = 0; j < ten_size; j++) {
            for (int k = 0; k < ten_size; k++) {
                for (int l = 0; l < ten_size; l++) {
                    float val = spatial_val + i + j + k + l + w_val;
                    tensor[i][j][k][l] = cfloat(val, val / 10.0);
                }
            }
        }
    }
    return tensor;
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

// Create 3D tensor data for all spatial points (w-k ordering)
static vector<vector<vector<vector<cfloat>>>> create_tensor3_data(int dim, int pnts, int ten_size, vector<float> w_points = {}) {
    vector<vector<vector<vector<cfloat>>>> data;
    int w_pts = w_points.empty() ? 1 : w_points.size();

    if (dim == 1) {
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
                float w_val = w_points.empty() ? 0.0 : w_points[w];
                Vec point = get_vec(i, 0, 0, pnts);
                data.push_back(create_tensor3(point, dim, ten_size, w_val));
            }
        }
    }
    if (dim == 2) {
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
                for (int j = 0; j < pnts; j++) {
                    float w_val = w_points.empty() ? 0.0 : w_points[w];
                    Vec point = get_vec(i, j, 0, pnts);
                    data.push_back(create_tensor3(point, dim, ten_size, w_val));
                }
            }
        }
    }
    return data;
}

// Create 4D tensor data for all spatial points (w-k ordering)
static vector<vector<vector<vector<vector<cfloat>>>>> create_tensor4_data(int dim, int pnts, int ten_size, vector<float> w_points = {}) {
    vector<vector<vector<vector<vector<cfloat>>>>> data;
    int w_pts = w_points.empty() ? 1 : w_points.size();

    if (dim == 1) {
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
                float w_val = w_points.empty() ? 0.0 : w_points[w];
                Vec point = get_vec(i, 0, 0, pnts);
                data.push_back(create_tensor4(point, dim, ten_size, w_val));
            }
        }
    }
    if (dim == 2) {
        for (int w = 0; w < w_pts; w++) {
            for (int i = 0; i < pnts; i++) {
                for (int j = 0; j < pnts; j++) {
                    float w_val = w_points.empty() ? 0.0 : w_points[w];
                    Vec point = get_vec(i, j, 0, pnts);
                    data.push_back(create_tensor4(point, dim, ten_size, w_val));
                }
            }
        }
    }
    return data;
}

// Test 3D tensor field - 1D spatial
static bool tensor3_field_1d() {
    string filename = "/tmp/test_tensor3_1d.h5";
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_tensor3_data(1, mpts, ten3_dim);

    // Create Field_CM - for n_indices=3, dim_indices represents tensor dimensions
    Field_CM field(data, {ten3_dim, ten3_dim, ten3_dim}, mesh, domain);

    // Test evaluation at x=0.5 (center)
    Vec v(0.0);
    auto result = field(v);

    // Result should be flattened from [2][2][2] to [4][2]
    // At x=0.5: spatial_val = 0.5, T_000 = 0.5 + 0 + 0 + 0 = 0.5
    if (result.size() != ten3_dim * ten3_dim) {
        cout << "Test tensor3_field_1d failed: wrong outer dimension " << result.size() << endl;
        return false;
    }
    if (result[0].size() != ten3_dim) {
        cout << "Test tensor3_field_1d failed: wrong inner dimension " << result[0].size() << endl;
        return false;
    }

    cfloat expected(0.5, 0.05);
    float tol = 0.01;
    if (abs(result[0][0] - expected) > tol) {
        cout << "Test tensor3_field_1d failed: result[0][0]=" << result[0][0] << " expected=" << expected << endl;
        return false;
    }

    // Test save/load
    field.save(filename);
    Field_CM loaded_field(filename);
    auto loaded_result = loaded_field(v);
    if (abs(loaded_result[0][0] - expected) > tol) {
        cout << "Test tensor3_field_1d (loaded) failed: result[0][0]=" << loaded_result[0][0] << " expected=" << expected << endl;
        return false;
    }

    filesystem::remove(filename);
    return true;
}

// Test 3D tensor field - 2D spatial
static bool tensor3_field_2d() {
    string filename = "/tmp/test_tensor3_2d.h5";
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    auto data = create_tensor3_data(2, mpts, ten3_dim);

    Field_CM field(data, {ten3_dim, ten3_dim, ten3_dim}, mesh, domain);

    Vec v(0.0, 0.0);
    auto result = field(v);

    // At (0.5, 0.5): spatial_val = 1.0, T_000 = 1.0 + 0 + 0 + 0 = 1.0
    cfloat expected(1.0, 0.1);
    float tol = 0.01;
    if (abs(result[0][0] - expected) > tol) {
        cout << "Test tensor3_field_2d failed: result[0][0]=" << result[0][0] << " expected=" << expected << endl;
        return false;
    }

    // Test save/load
    field.save(filename);
    Field_CM loaded_field(filename);
    auto loaded_result = loaded_field(v);
    if (abs(loaded_result[0][0] - expected) > tol) {
        cout << "Test tensor3_field_2d (loaded) failed: result[0][0]=" << loaded_result[0][0] << " expected=" << expected << endl;
        return false;
    }

    filesystem::remove(filename);
    return true;
}

// Test 4D tensor field - 1D spatial
static bool tensor4_field_1d() {
    string filename = "/tmp/test_tensor4_1d.h5";
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_tensor4_data(1, mpts, ten4_dim);

    Field_CM field(data, {ten4_dim, ten4_dim, ten4_dim, ten4_dim}, mesh, domain);

    Vec v(0.0);
    auto result = field(v);

    // Result should be flattened from [2][2][2][2] to [4][4]
    // At x=0.5: spatial_val = 0.5, T_0000 = 0.5 + 0 + 0 + 0 + 0 = 0.5
    if (result.size() != ten4_dim * ten4_dim) {
        cout << "Test tensor4_field_1d failed: wrong outer dimension " << result.size() << endl;
        return false;
    }
    if (result[0].size() != ten4_dim * ten4_dim) {
        cout << "Test tensor4_field_1d failed: wrong inner dimension " << result[0].size() << endl;
        return false;
    }

    cfloat expected(0.5, 0.05);
    float tol = 0.01;
    if (abs(result[0][0] - expected) > tol) {
        cout << "Test tensor4_field_1d failed: result[0][0]=" << result[0][0] << " expected=" << expected << endl;
        return false;
    }

    // Test save/load
    field.save(filename);
    Field_CM loaded_field(filename);
    auto loaded_result = loaded_field(v);
    if (abs(loaded_result[0][0] - expected) > tol) {
        cout << "Test tensor4_field_1d (loaded) failed: result[0][0]=" << loaded_result[0][0] << " expected=" << expected << endl;
        return false;
    }

    filesystem::remove(filename);
    return true;
}

// Test 4D tensor field - 2D spatial
static bool tensor4_field_2d() {
    string filename = "/tmp/test_tensor4_2d.h5";
    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    auto data = create_tensor4_data(2, mpts, ten4_dim);

    Field_CM field(data, {ten4_dim, ten4_dim, ten4_dim, ten4_dim}, mesh, domain);

    Vec v(0.0, 0.0);
    auto result = field(v);

    // At (0.5, 0.5): spatial_val = 1.0, T_0000 = 1.0 + 0 + 0 + 0 + 0 = 1.0
    cfloat expected(1.0, 0.1);
    float tol = 0.01;
    if (abs(result[0][0] - expected) > tol) {
        cout << "Test tensor4_field_2d failed: result[0][0]=" << result[0][0] << " expected=" << expected << endl;
        return false;
    }

    // Test save/load
    field.save(filename);
    Field_CM loaded_field(filename);
    auto loaded_result = loaded_field(v);
    if (abs(loaded_result[0][0] - expected) > tol) {
        cout << "Test tensor4_field_2d (loaded) failed: result[0][0]=" << loaded_result[0][0] << " expected=" << expected << endl;
        return false;
    }

    filesystem::remove(filename);
    return true;
}

// Test 4D tensor with frequency dimension
static bool tensor4_field_w() {
    string filename = "/tmp/test_tensor4_w.h5";
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {0.0, 1.0, 2.0};
    auto data = create_tensor4_data(1, mpts, ten4_dim, w_points);

    Field_CM field(data, {ten4_dim, ten4_dim, ten4_dim, ten4_dim}, mesh, domain, w_points);

    // Test at w=0.0
    Vec v(0.0);
    auto result = field(v, 0.0);
    cfloat expected0(0.5, 0.05);  // spatial_val=0.5, w=0.0
    float tol = 0.01;
    if (abs(result[0][0] - expected0) > tol) {
        cout << "Test tensor4_field_w (w=0) failed: result[0][0]=" << result[0][0] << " expected=" << expected0 << endl;
        return false;
    }

    // Test at w=2.0
    result = field(v, 2.0);
    cfloat expected2(2.5, 0.25);  // spatial_val=0.5, w=2.0
    if (abs(result[0][0] - expected2) > tol) {
        cout << "Test tensor4_field_w (w=2) failed: result[0][0]=" << result[0][0] << " expected=" << expected2 << endl;
        return false;
    }

    // Test at w=1.0 (should interpolate between w=0 and w=2)
    result = field(v, 1.0);
    cfloat expected1(1.5, 0.15);  // spatial_val=0.5, w=1.0
    if (abs(result[0][0] - expected1) > tol) {
        cout << "Test tensor4_field_w (w=1) failed: result[0][0]=" << result[0][0] << " expected=" << expected1 << endl;
        return false;
    }

    // Test save/load
    field.save(filename);
    Field_CM loaded_field(filename);
    auto loaded_result = loaded_field(v, 1.0);
    if (abs(loaded_result[0][0] - expected1) > tol) {
        cout << "Test tensor4_field_w (loaded) failed: result[0][0]=" << loaded_result[0][0] << " expected=" << expected1 << endl;
        return false;
    }

    filesystem::remove(filename);
    return true;
}

// Test saving and loading 4D tensor field
static bool tensor4_save_load() {
    string filename = "/tmp/test_tensor4.h5";
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    auto data = create_tensor4_data(1, mpts, ten4_dim);

    // Create and save
    {
        Field_CM field(data, {ten4_dim, ten4_dim, ten4_dim, ten4_dim}, mesh, domain);
        field.save(filename);
    }

    // Load and test
    {
        Field_CM field(filename);
        Vec v(0.0);
        auto result = field(v);

        cfloat expected(0.5, 0.05);
        float tol = 0.01;
        if (abs(result[0][0] - expected) > tol) {
            cout << "Test tensor4_save_load failed: result[0][0]=" << result[0][0] << " expected=" << expected << endl;
            return false;
        }
    }

    // Cleanup
    filesystem::remove(filename);
    return true;
}

} // namespace tensor_field_tests_ns

bool tensor_field_tests() {
    using namespace tensor_field_tests_ns;

    int passed = 0;
    int total = 6;

    if (tensor3_field_1d()) passed++; else cout << "tensor3_field_1d FAILED\n";
    if (tensor3_field_2d()) passed++; else cout << "tensor3_field_2d FAILED\n";
    if (tensor4_field_1d()) passed++; else cout << "tensor4_field_1d FAILED\n";
    if (tensor4_field_2d()) passed++; else cout << "tensor4_field_2d FAILED\n";
    if (tensor4_field_w()) passed++; else cout << "tensor4_field_w FAILED\n";
    if (tensor4_save_load()) passed++; else cout << "tensor4_save_load FAILED\n";

    if (passed == total) {
        cout << "\033[1;32mAll " << total << " Tensor Field tests passed!\n\033[0m";
        return true;
    } else {
        cout << "\033[1;31m" << (total - passed) << " Tensor Field tests failed!\n\033[0m";
        return false;
    }
}
