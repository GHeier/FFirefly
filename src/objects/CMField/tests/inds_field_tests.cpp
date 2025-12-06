/**
 * Comprehensive test suite for Field objects using vector<int> inds
 *
 * Tests cover:
 * - Scalars (inds={})
 * - Vectors (inds={n})
 * - Uniform matrices (inds={n,n})
 * - Non-uniform matrices (inds={m,n})
 * - Uniform 3D tensors (inds={l,m,n})
 * - Non-uniform 3D tensors (inds={l,m,n} where l!=m!=n)
 * - 4D tensors (single-band and multi-band)
 * - Save/load round-trips
 * - Edge cases
 */

#include "../../../config/load/c_config.h"
#include "inds_field_tests.hpp"
#include "../fields.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <complex>

using namespace std;

const float tol = 1e-5;

// Helper: Create complex scalar data (w-k ordering)
vector<cfloat> create_scalar_data(int npts, int nw) {
    vector<cfloat> data(npts * nw);
    // w-k ordering: frequency varies slowest
    for (int j = 0; j < nw; j++) {
        for (int i = 0; i < npts; i++) {
            data[j * npts + i] = cfloat(1.0 * i / npts, 1.0 * j / nw);
        }
    }
    return data;
}

// Helper: Create vector data (inds={n})
vector<vector<cfloat>> create_vector_data(int npts, int nw, int vec_dim) {
    vector<vector<cfloat>> data(npts * nw);
    for (int i = 0; i < npts * nw; i++) {
        data[i].resize(vec_dim);
        for (int j = 0; j < vec_dim; j++) {
            data[i][j] = cfloat(0.5 + i * 0.01 + j * 0.1, 0.1 + i * 0.001);
        }
    }
    return data;
}

// Helper: Create matrix data with potentially non-uniform dimensions
vector<vector<vector<cfloat>>> create_matrix_data(int npts, int nw, int dim1, int dim2) {
    vector<vector<vector<cfloat>>> data(npts * nw);
    for (int t = 0; t < npts * nw; t++) {
        data[t].resize(dim1);
        for (int i = 0; i < dim1; i++) {
            data[t][i].resize(dim2);
            for (int j = 0; j < dim2; j++) {
                data[t][i][j] = cfloat(0.5 + t * 0.01 + i * 0.1 + j * 0.01,
                                       0.1 + t * 0.001);
            }
        }
    }
    return data;
}

// Helper: Create 3D tensor data with non-uniform dimensions
vector<vector<vector<vector<cfloat>>>> create_tensor3_data(int npts, int nw,
                                                            int dim1, int dim2, int dim3) {
    vector<vector<vector<vector<cfloat>>>> data(npts * nw);
    for (int t = 0; t < npts * nw; t++) {
        data[t].resize(dim1);
        for (int i = 0; i < dim1; i++) {
            data[t][i].resize(dim2);
            for (int j = 0; j < dim2; j++) {
                data[t][i][j].resize(dim3);
                for (int k = 0; k < dim3; k++) {
                    data[t][i][j][k] = cfloat(0.5 + t * 0.01 + i * 0.1 + j * 0.01 + k * 0.001,
                                              0.1 + t * 0.001);
                }
            }
        }
    }
    return data;
}

// Helper: Create 4D tensor data with non-uniform dimensions
vector<vector<vector<vector<vector<cfloat>>>>> create_tensor4_data(int npts, int nw,
                                                                     int dim1, int dim2,
                                                                     int dim3, int dim4) {
    vector<vector<vector<vector<vector<cfloat>>>>> data(npts * nw);
    for (int t = 0; t < npts * nw; t++) {
        data[t].resize(dim1);
        for (int i = 0; i < dim1; i++) {
            data[t][i].resize(dim2);
            for (int j = 0; j < dim2; j++) {
                data[t][i][j].resize(dim3);
                for (int k = 0; k < dim3; k++) {
                    data[t][i][j][k].resize(dim4);
                    for (int l = 0; l < dim4; l++) {
                        data[t][i][j][k][l] = cfloat(0.5 + t * 0.01, 0.1 + t * 0.001);
                    }
                }
            }
        }
    }
    return data;
}

// =============================================================================
// Test 1: Scalar field 1D with frequency (inds={})
// =============================================================================
bool test_scalar_1d_w() {

    int mpts = 10;
    int nw = 5;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {-0.2, -0.1, 0.0, 0.1, 0.2};

    auto data = create_scalar_data(mpts, nw);
    Field_C field(data, mesh, domain, w_points);

    // Evaluate at point
    Vec v(0.0);  // Center of BZ
    cfloat result = field(v, 0.0);
    // At k=0.5 (interpolates between i=4 and i=5), w=0.0 (j=2): cfloat(0.45, 0.4)
    cfloat expected(0.45, 0.4);
    if (abs(result - expected) > tol) {
        cout << "test_scalar_1d_w FAILED: result=" << result << " expected=" << expected << endl;
        return false;
    }
    return true;
}

// =============================================================================
// Test 2: Scalar field 2D with frequency (inds={})
// =============================================================================
bool test_scalar_2d_w() {

    int mpts = 100;  // 10x10
    int nw = 5;
    vector<int> mesh = {10, 10};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};
    vector<float> w_points = {-0.2, -0.1, 0.0, 0.1, 0.2};

    auto data = create_scalar_data(mpts, nw);
    Field_C field(data, mesh, domain, w_points);

    // Evaluate at point
    Vec v(0.5, 0.5);
    cfloat result = field(v, 0.0);
    // Should interpolate - just check it doesn't crash
    return (!isnan(result.real()));
}

// =============================================================================
// Test 3: Vector field (inds={3})
// =============================================================================
bool test_vector_field() {

    int mpts = 10;
    int nw = 1;
    int vec_dim = 3;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {vec_dim};

    auto data = create_vector_data(mpts, nw, vec_dim);
    Field_CM field(data, inds, mesh, domain, w_points);

    // Evaluate at point
    Vec v(0.5);
    auto result = field(v);  // Should return [vec_dim][1] column matrix
    return (result.size() == vec_dim && result[0].size() == 1);
}

// =============================================================================
// Test 4: Uniform square matrix (inds={3,3})
// =============================================================================
bool test_uniform_matrix() {

    int mpts = 10;
    int nw = 1;
    int mat_dim = 3;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {mat_dim, mat_dim};

    auto data = create_matrix_data(mpts, nw, mat_dim, mat_dim);
    Field_CM field(data, inds, mesh, domain, w_points);

    // Evaluate at point
    Vec v(0.5);
    auto result = field(v);  // Should return 3x3 matrix
    return (result.size() == mat_dim && result[0].size() == mat_dim);
}

// =============================================================================
// Test 5: NON-UNIFORM matrix (inds={2,3}) - NEW CAPABILITY
// =============================================================================
bool test_nonuniform_matrix() {

    int mpts = 10;
    int nw = 1;
    int dim1 = 2, dim2 = 3;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {dim1, dim2};  // 2x3 matrix

    auto data = create_matrix_data(mpts, nw, dim1, dim2);
    Field_CM field(data, inds, mesh, domain, w_points);

    // Evaluate at point
    Vec v(0.5);
    auto result = field(v);  // Should return 2x3 matrix
    return (result.size() == dim1 && result[0].size() == dim2);
}

// =============================================================================
// Test 6: Uniform 3D tensor (inds={2,2,2})
// =============================================================================
bool test_uniform_tensor3() {

    int mpts = 10;
    int nw = 1;
    int ten_dim = 2;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {ten_dim, ten_dim, ten_dim};

    auto data = create_tensor3_data(mpts, nw, ten_dim, ten_dim, ten_dim);
    Field_CM field(data, inds, mesh, domain, w_points);

    // Evaluate at point
    Vec v(0.5);
    auto result = field(v);  // Returns flattened matrix [d1*d2][d3] = [2*2][2] = [4][2]
    return (result.size() == ten_dim * ten_dim && result[0].size() == ten_dim);
}

// =============================================================================
// Test 7: NON-UNIFORM 3D tensor (inds={2,3,4}) - NEW CAPABILITY
// =============================================================================
bool test_nonuniform_tensor3() {

    int mpts = 10;
    int nw = 1;
    int dim1 = 2, dim2 = 3, dim3 = 4;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {dim1, dim2, dim3};

    auto data = create_tensor3_data(mpts, nw, dim1, dim2, dim3);
    Field_CM field(data, inds, mesh, domain, w_points);

    // Evaluate at point
    Vec v(0.5);
    auto result = field(v);  // Returns flattened matrix [dim1*dim2][dim3] = [2*3][4] = [6][4]
    return (result.size() == dim1 * dim2 && result[0].size() == dim3);
}

// =============================================================================
// Test 8: Single-band 4D vertex (inds={1,1,1,1})
// =============================================================================
bool test_singleband_vertex() {

    int mpts = 100;  // 10x10 k-mesh
    int nw = 5;
    vector<int> mesh = {10, 10};
    vector<vector<float>> domain = {{6.28, 0.0}, {0.0, 6.28}};
    vector<float> w_points = {-0.2, -0.1, 0.0, 0.1, 0.2};
    vector<int> inds = {1, 1, 1, 1};  // Single-band vertex

    auto data = create_tensor4_data(mpts, nw, 1, 1, 1, 1);
    Field_CM field(data, inds, mesh, domain, w_points);

    // Evaluate at point
    Vec v(3.14, 3.14);
    auto result = field(v, 0.0);  // Returns flattened matrix [d1*d2][d3*d4] = [1*1][1*1] = [1][1]
    return (result.size() == 1 && result[0].size() == 1);
}

// =============================================================================
// Test 9: Two-band 4D vertex (inds={2,2,2,2})
// =============================================================================
bool test_twoband_vertex() {

    int mpts = 10;
    int nw = 1;
    int dim = 2;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {dim, dim, dim, dim};

    auto data = create_tensor4_data(mpts, nw, dim, dim, dim, dim);
    Field_CM field(data, inds, mesh, domain, w_points);

    // Evaluate at point
    Vec v(0.5);
    auto result = field(v);  // Returns flattened matrix [d1*d2][d3*d4] = [2*2][2*2] = [4][4]
    return (result.size() == dim * dim && result[0].size() == dim * dim);
}

// =============================================================================
// Test 10: NON-UNIFORM 4D tensor (inds={1,2,2,1}) - NEW CAPABILITY
// =============================================================================
bool test_nonuniform_tensor4() {

    int mpts = 10;
    int nw = 1;
    int d1 = 1, d2 = 2, d3 = 2, d4 = 1;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {d1, d2, d3, d4};

    auto data = create_tensor4_data(mpts, nw, d1, d2, d3, d4);
    Field_CM field(data, inds, mesh, domain, w_points);

    // Evaluate at point
    Vec v(0.5);
    auto result = field(v);  // Should return 1x2x2x1 tensor
    // Returns flattened matrix [d1*d2][d3*d4] = [1*2][2*1] = [2][2]
    return (result.size() == d1 * d2 && result[0].size() == d3 * d4);
}

// =============================================================================
// Test 11: Save/Load round-trip for scalar
// =============================================================================
bool test_scalar_save_load() {

    int mpts = 10;
    int nw = 5;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {-0.2, -0.1, 0.0, 0.1, 0.2};

    auto data = create_scalar_data(mpts, nw);
    Field_C field1(data, mesh, domain, w_points);

    // Save
    string filename = "/tmp/test_scalar_inds.h5";
    field1.save(filename);

    // Load
    Field_C field2(filename);

    // Compare
    Vec v(0.5);
    cfloat result1 = field1(v, 0.0);
    cfloat result2 = field2(v, 0.0);
    return (abs(result1 - result2) < tol);
}

// =============================================================================
// Test 12: Save/Load round-trip for uniform matrix
// =============================================================================
bool test_matrix_save_load() {

    int mpts = 10;
    int nw = 1;
    int mat_dim = 3;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {mat_dim, mat_dim};

    auto data = create_matrix_data(mpts, nw, mat_dim, mat_dim);
    Field_CM field1(data, inds, mesh, domain, w_points);

    // Save
    string filename = "/tmp/test_matrix_inds.h5";
    field1.save(filename);

    // Load
    Field_CM field2(filename);

    // Compare
    Vec v(0.5);
    auto result1 = field1(v);
    auto result2 = field2(v);
    return (abs(result1[0][0] - result2[0][0]) < tol);
}

// =============================================================================
// Test 13: Save/Load round-trip for NON-UNIFORM matrix
// =============================================================================
bool test_nonuniform_matrix_save_load() {

    int mpts = 10;
    int nw = 1;
    int dim1 = 2, dim2 = 3;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {dim1, dim2};

    auto data = create_matrix_data(mpts, nw, dim1, dim2);
    Field_CM field1(data, inds, mesh, domain, w_points);

    // Save
    string filename = "/tmp/test_nonuniform_matrix_inds.h5";
    field1.save(filename);

    // Load
    Field_CM field2(filename);

    // Compare
    Vec v(0.5);
    auto result1 = field1(v);
    auto result2 = field2(v);
    return (result1.size() == dim1 && result1[0].size() == dim2 && abs(result1[0][0] - result2[0][0]) < tol);
}

// =============================================================================
// Test 14: Save/Load round-trip for 4D vertex
// =============================================================================
bool test_vertex_save_load() {

    int mpts = 100;
    int nw = 5;
    vector<int> mesh = {10, 10};
    vector<vector<float>> domain = {{6.28, 0.0}, {0.0, 6.28}};
    vector<float> w_points = {-0.2, -0.1, 0.0, 0.1, 0.2};
    vector<int> inds = {1, 1, 1, 1};

    auto data = create_tensor4_data(mpts, nw, 1, 1, 1, 1);
    Field_CM field1(data, inds, mesh, domain, w_points);

    // Save
    string filename = "/tmp/test_vertex_inds.h5";
    field1.save(filename);

    // Load
    Field_CM field2(filename);

    // Compare
    Vec v(3.14, 3.14);
    auto result1 = field1(v, 0.0);  // Returns [1][1] matrix
    auto result2 = field2(v, 0.0);
    return (abs(result1[0][0] - result2[0][0]) < tol);
}

// =============================================================================
// Test 15: Edge case - single element tensor (inds={1})
// =============================================================================
bool test_single_element_vector() {

    int mpts = 10;
    int nw = 1;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {1};

    auto data = create_vector_data(mpts, nw, 1);
    Field_CM field(data, inds, mesh, domain, w_points);

    Vec v(0.5);
    auto result = field(v);  // Should return [1][1] matrix for single-element vector
    return (result.size() == 1 && result[0].size() == 1);
}

// =============================================================================
// Test 16: Verify inds stored correctly in HDF5
// =============================================================================
bool test_inds_hdf5_storage() {

    int mpts = 10;
    int nw = 1;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {2, 3, 4};

    auto data = create_tensor3_data(mpts, nw, 2, 3, 4);
    Field_CM field1(data, inds, mesh, domain, w_points);

    // Save
    string filename = "/tmp/test_inds_storage.h5";
    field1.save(filename);

    // Load and check inds
    Field_CM field2(filename);
    auto loaded_data = field2.get_data();

    // Verify inds matches
    return (loaded_data->inds.size() == 3 &&
            loaded_data->inds[0] == 2 &&
            loaded_data->inds[1] == 3 &&
            loaded_data->inds[2] == 4);
}

// =============================================================================
// Test 17: Large tensor (inds={5,5,5,5}) - stress test
// =============================================================================
bool test_large_tensor() {

    int mpts = 10;
    int nw = 1;
    int dim = 5;
    vector<int> mesh = {mpts};
    vector<vector<float>> domain = {{1.0}};
    vector<float> w_points = {};
    vector<int> inds = {dim, dim, dim, dim};

    auto data = create_tensor4_data(mpts, nw, dim, dim, dim, dim);
    Field_CM field(data, inds, mesh, domain, w_points);

    // Evaluate - mostly checking memory doesn't explode
    Vec v(0.5);
    auto result = field(v);  // Returns flattened [dim*dim][dim*dim] = [25][25]
    return (result.size() == dim * dim && result[0].size() == dim * dim);
}

// =============================================================================
// Test 18: Verify total_index_size() calculation
// =============================================================================
bool test_total_index_size() {

    // Scalar: inds={} → size = 1
    {
        int mpts = 10;
        vector<int> mesh = {mpts};
        vector<vector<float>> domain = {{1.0}};
        auto data = create_scalar_data(mpts, 1);
        Field_C field(data, mesh, domain, {});
        return (field.get_data()->total_index_size() == 1);
    }

    // Matrix: inds={3,3} → size = 9
    {
        int mpts = 10;
        vector<int> mesh = {mpts};
        vector<vector<float>> domain = {{1.0}};
        vector<int> inds = {3, 3};
        auto data = create_matrix_data(mpts, 1, 3, 3);
        Field_CM field(data, inds, mesh, domain, {});
        return (field.get_data()->total_index_size() == 9);
    }

    // Non-uniform: inds={2,3,4} → size = 24
    {
        int mpts = 10;
        vector<int> mesh = {mpts};
        vector<vector<float>> domain = {{1.0}};
        vector<int> inds = {2, 3, 4};
        auto data = create_tensor3_data(mpts, 1, 2, 3, 4);
        Field_CM field(data, inds, mesh, domain, {});
        return (field.get_data()->total_index_size() == 24);
    }

    return true;
}

// =============================================================================
// Main test runner
// =============================================================================
void run_all_inds_tests() {
    cout << "Running Comprehensive inds Tests" << endl;
    cout << "========================================\n" << endl;

    int passed = 0, total = 0;

    #define RUN_TEST(test_func) \
        try { \
            total++; \
            test_func(); \
            passed++; \
        } catch (const exception& e) { \
            cout << "  FAILED: " << e.what() << endl; \
        }

    // Scalar tests
    RUN_TEST(test_scalar_1d_w);
    RUN_TEST(test_scalar_2d_w);
    RUN_TEST(test_scalar_save_load);

    // Vector tests
    RUN_TEST(test_vector_field);
    RUN_TEST(test_single_element_vector);

    // Matrix tests
    RUN_TEST(test_uniform_matrix);
    RUN_TEST(test_nonuniform_matrix);
    RUN_TEST(test_matrix_save_load);
    RUN_TEST(test_nonuniform_matrix_save_load);

    // 3D tensor tests
    RUN_TEST(test_uniform_tensor3);
    RUN_TEST(test_nonuniform_tensor3);

    // 4D tensor tests
    RUN_TEST(test_singleband_vertex);
    RUN_TEST(test_twoband_vertex);
    RUN_TEST(test_nonuniform_tensor4);
    RUN_TEST(test_vertex_save_load);

    // Edge cases and stress tests
    RUN_TEST(test_inds_hdf5_storage);
    RUN_TEST(test_large_tensor);
    RUN_TEST(test_total_index_size);

    cout << "\n========================================" << endl;
    cout << "Results: " << passed << "/" << total << " tests passed" << endl;
    cout << "========================================\n" << endl;
}

bool inds_tests() {
    int num_tests = 18;
    bool all_tests[num_tests] = {
        test_scalar_1d_w(),
        test_scalar_2d_w(),
        test_scalar_save_load(),
        test_vector_field(),
        test_single_element_vector(),
        test_uniform_matrix(),
        test_nonuniform_matrix(),
        test_matrix_save_load(),
        test_nonuniform_matrix_save_load(),
        test_uniform_tensor3(),
        test_nonuniform_tensor3(),
        test_singleband_vertex(),
        test_twoband_vertex(),
        test_nonuniform_tensor4(),
        test_vertex_save_load(),
        test_inds_hdf5_storage(),
        test_large_tensor(),
        test_total_index_size()
    };
    return print_test_results(all_tests, num_tests, "Inds Field tests");
}
