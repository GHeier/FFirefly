#include "../base_data.hpp"
#include "../data_evaluator.hpp"
#include <filesystem>
#include <iostream>

#include "../../../config/load/c_config.h"

using namespace std;

int mpts = 3;

cfloat func_linear(Vec p, int dim) {
    float val = 0.0;
    for (int i = 0; i < dim; i++) 
        val += p(i);
    return cfloat(val, val / 10);
}

float get_pnt(int i, int pnts) {
    return 1.0 * i / (pnts - 1);
}

Vec get_vec(int i, int j, int k, int pnts) {
    return Vec(
            get_pnt(i, pnts),
            get_pnt(j, pnts),
            get_pnt(k, pnts)
            );
}

vector<cfloat> create_data(int dim, int pnts, vector<float> w_points = {}) {
    vector<cfloat> values;
    int w_pts = w_points.empty() ? 1 : w_points.size();
    int idx = 0;

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

vector<cfloat> create_data_tensor(int dim, int pnts, vector<float> w_points = {}, int mat_dim = 1, int n_inds = 0) {
    vector<cfloat> values;
    int w_pts = w_points.empty() ? 1 : w_points.size();
    int idx = 0;

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
                        for (int n = 0; n < n_inds; n++) {
                            values.push_back(value);
                        }
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

bool evaluate_real_scalar_1d_w() {
    // 1. Create a simple BaseData
    BaseData data;
    data.is_complex = false;
    data.is_vector = false;
    data.with_k = false;
    data.with_w = true;
    data.as_mesh = false;
    data.dimension = 1; // 3-component vector
    vector<float> w_points = {1, 2, 3};
    vector<vector<float>> domain = {{1.0}};

    // Fill with 2x3 complex values
    vector<cfloat> vecs = {
        cfloat(1.0, 0.1), cfloat(2.0, 0.2), cfloat(3.0, 0.3)
    };

    // Store as matrix of scalars (variant type)
    data.w_points = w_points;
    data.domain = domain;
    data.data = vecs;
    DataEvaluator field(data);
    auto result = field(1.5);
    if (auto* s = std::get_if<float>(&result)) {
            return fabs(*s - 1.5) < 1e-6;
    }
    return false;
}

bool evaluate_complex_scalar_1d_w() {
    // 1. Create a simple BaseData
    BaseData data;
    data.is_complex = true;
    data.is_vector = false;
    data.with_k = false;
    data.with_w = true;
    data.as_mesh = false;
    data.dimension = 1; // 3-component vector
    vector<float> w_points = {1, 2, 3};
    vector<vector<float>> domain = {{1.0}};

    // Fill with 2x3 complex values
    vector<cfloat> vecs = {
        cfloat(1.0, 0.1), cfloat(2.0, 0.2), cfloat(3.0, 0.3)
    };

    // Store as matrix of scalars (variant type)
    data.w_points = w_points;
    data.domain = domain;
    data.data = vecs;
    DataEvaluator field(data);
    auto result = field(1.5);
    if (auto* s = std::get_if<cfloat>(&result)) {
            cfloat ans(1.5, 0.15);
            return fabs(*s - ans) < 1e-6;
    }
    return false;
}

bool evaluate_real_scalar_1d_k() {
    // 1. Create a simple BaseData
    BaseData data;
    data.is_complex = false;
    data.is_vector = false;
    data.with_k = true;
    data.with_w = false;
    data.as_mesh = true;
    data.mesh = {mpts};   // 3 k-points
    data.dimension = 1; // 3-component vector
    vector<float> w_points;
    vector<vector<float>> domain = {{1.0}};

    vector<cfloat> vecs = create_data(data.dimension, mpts);

    // Store as matrix of scalars (variant type)
    data.w_points = w_points;
    data.domain = domain;
    data.data = vecs;
    DataEvaluator field(data);
    Vec v(0.5);
    auto result = field(v);
    if (auto* s = std::get_if<float>(&result)) {
            return fabs(*s - 0.5) < 1e-6;
    }
    return false;
}

bool evaluate_complex_scalar_1d_k() {
    BaseData data;
    data.is_complex = true;
    data.is_vector = false;
    data.with_k = true;
    data.with_w = false;
    data.as_mesh = true;
    data.mesh = {mpts};
    data.dimension = 1;
    vector<float> w_points;
    vector<vector<float>> domain = {{1.0}};

    vector<cfloat> vecs = create_data(data.dimension, mpts);

    data.w_points = w_points;
    data.domain = domain;
    data.data = vecs;
    DataEvaluator field(data);
    Vec v(0.5);
    auto result = field(v);
    if (auto* s = std::get_if<cfloat>(&result)) {
        return fabs(*s - cfloat(0.5, 0.05)) < 1e-6;
    }
    return false;
}

bool evaluate_real_scalar_2d_k() {
    BaseData data;
    data.is_complex = false;
    data.is_vector = false;
    data.with_k = true;
    data.with_w = false;
    data.as_mesh = true;
    data.dimension = 2;
    data.mesh = {mpts, mpts};   // 3 k-points

    data.domain = {
        {1.0, 0.0}, {0.0, 1.0}
    };
    vector<cfloat> temp = create_data(data.dimension, mpts);
    data.data = temp;

    DataEvaluator field(data);
    Vec v(0.5, 0.5);
    auto result = field(v);
    if (auto* s = std::get_if<float>(&result)) {
        return fabs(*s - 1.0) < 1e-6;
    }
    return false;
}

bool evaluate_complex_scalar_2d_k() {
    BaseData data;
    data.is_complex = true;
    data.is_vector = false;
    data.with_k = true;
    data.with_w = false;
    data.as_mesh = true;
    data.dimension = 2;
    data.mesh = {mpts, mpts};

    data.domain = {
        {1.0, 0.0}, {0.0, 1.0}
    };
    vector<cfloat> temp = create_data(data.dimension, mpts);
    data.data = temp;

    DataEvaluator field(data);
    Vec v(0.5, 0.5);
    auto result = field(v);
    if (auto* s = std::get_if<cfloat>(&result)) {
        return fabs(*s - cfloat(1.0, 0.1)) < 1e-6;
    }
    return false;
}

bool evaluate_real_scalar_2d_w() {
    BaseData data;
    data.is_complex = false;
    data.is_vector = false;
    data.with_k = false;
    data.with_w = true;
    data.as_mesh = true;
    data.mesh = {mpts, mpts};
    data.dimension = 2;

    data.domain = {
        {1.0, 0.0}, {0.0, 1.0}
    };
    data.w_points = {1.0, 2.0, 3.0};
    data.data = create_data(data.dimension, mpts, data.w_points);

    DataEvaluator field(data);
    Vec p(0.6, 0.6);
    auto result = field(p, 1.1);
    if (auto* s = std::get_if<float>(&result)) {
        return fabs(*s - 2.3) < 1e-6;
    }
    return false;
}

bool evaluate_complex_scalar_2d_w() {
    BaseData data;
    data.is_complex = true;
    data.is_vector = false;
    data.with_k = false;
    data.with_w = true;
    data.as_mesh = true;
    data.mesh = {mpts, mpts};   // 3 k-points
    data.dimension = 2;

    data.domain = {
        {1.0, 0.0}, {0.0, 1.0}
    };
    data.w_points = {1.0, 2.0, 3.0};
    data.data = create_data(data.dimension, mpts, data.w_points);


    DataEvaluator field(data);
    Vec p(0.6, 0.6);
    auto result = field(p, 1.1);
    if (auto* s = std::get_if<cfloat>(&result)) {
        return fabs(*s - cfloat(2.3, 0.23)) < 1e-6;
    }
    return false;
}

bool evaluate_real_scalar_3d_k() {
    BaseData data;
    data.is_complex = false;
    data.is_vector = false;
    data.with_k = true;
    data.with_w = false;
    data.as_mesh = true;
    data.dimension = 3;
    data.mesh = {mpts, mpts, mpts};

    data.domain = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    };
    vector<cfloat> temp = create_data(data.dimension, mpts);
    data.data = temp;

    DataEvaluator field(data);
    Vec v(0.5, 0.5, 0.5);
    auto result = field(v);
    if (auto* s = std::get_if<float>(&result)) {
        return fabs(*s - 1.5) < 1e-6;
    }
    return false;
}

bool evaluate_complex_scalar_3d_k() {
    BaseData data;
    data.is_complex = true;
    data.is_vector = false;
    data.with_k = true;
    data.with_w = false;
    data.as_mesh = true;
    data.dimension = 3;
    data.mesh = {mpts, mpts, mpts};

    data.domain = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    };
    vector<cfloat> temp = create_data(data.dimension, mpts);
    data.data = temp;

    DataEvaluator field(data);
    Vec v(0.5, 0.5, 0.5);
    auto result = field(v);
    if (auto* s = std::get_if<cfloat>(&result)) {
        return fabs(*s - cfloat(1.5, 0.15)) < 1e-6;
    }
    return false;
}

bool evaluate_real_scalar_3d_w() {
    BaseData data;
    data.is_complex = false;
    data.is_vector = false;
    data.with_k = false;
    data.with_w = true;
    data.as_mesh = true;
    data.mesh = {mpts, mpts, mpts};
    data.dimension = 3;

    data.domain = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    };
    data.w_points = {1.0, 2.0, 3.0};
    data.data = create_data(data.dimension, mpts, data.w_points);

    DataEvaluator field(data);
    Vec p(0.25, 0.25, 0.25);
    auto result = field(p, 1.5);
    if (auto* s = std::get_if<float>(&result)) {
        return fabs(*s - 2.25) < 1e-6;
    }
    return false;
}

bool evaluate_complex_scalar_3d_w() {
    BaseData data;
    data.is_complex = true;
    data.is_vector = false;
    data.with_k = false;
    data.with_w = true;
    data.as_mesh = true;
    data.mesh = {mpts, mpts, mpts};
    data.dimension = 3;

    data.domain = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    };
    data.w_points = {1.0, 2.0, 3.0};
    data.data = create_data(data.dimension, mpts, data.w_points);

    DataEvaluator field(data);
    Vec p(0.25, 0.25, 0.25);
    auto result = field(p, 1.5);
    if (auto* s = std::get_if<cfloat>(&result)) {
        return fabs(*s - cfloat(2.25, 0.225)) < 1e-6;
    }
    return false;
}

bool create_destroy() {
    // 1. Create a simple BaseData
    BaseData field;
    field.is_complex = true;
    field.is_vector = true;
    field.with_k = true;
    field.with_w = false;
    field.as_mesh = true;

    field.mesh = {2};   // 2 k-points
    field.dimension = 3; // 3-component vector
    field.domain = {{1.0}};
    field.inds = {3};   // 3-component vector (rank-1)

    // Fill with 2x3 complex values
    vector<vector<cfloat>> vecs = {
        {cfloat(1.0, 0.1), cfloat(2.0, 0.2), cfloat(3.0, 0.3)},
        {cfloat(4.0, 0.4), cfloat(5.0, 0.5), cfloat(6.0, 0.6)}
    };

    // Store as matrix of scalars (variant type)
    field.data = vecs;
    auto& mat1 = field.get<vector<vector<cfloat>>>();
    if (mat1.size() != vecs.size()) return false;
    if (mat1[0].size() != vecs[0].size()) return false;

    for (size_t i = 0; i < mat1.size(); ++i) {
        for (size_t j = 0; j < mat1[i].size(); ++j) {
            if (abs(mat1[i][j] - vecs[i][j]) > 1e-6f) return false;
        }
    }

    // 2. Save to file
    string fname = "testfield.h5";
    save_data_to_hdf5(field, fname);

    // 3. Reload
    BaseData loaded = load_data_from_hdf5(fname);

    // 4. Verify contents
    auto& mat = loaded.get<vector<vector<cfloat>>>();
    if (mat.size() != vecs.size()) return false;
    if (mat[0].size() != vecs[0].size()) return false;

    for (size_t i = 0; i < mat.size(); ++i) {
        for (size_t j = 0; j < mat[i].size(); ++j) {
            if (abs(mat[i][j] - vecs[i][j]) > 1e-6f) return false;
        }
    }

    return true;
}

bool point_storage_scalar() {
    // Test point storage with as_mesh = false
    BaseData data;
    data.is_complex = true;
    data.is_vector = false;
    data.with_k = true;
    data.with_w = false;
    data.as_mesh = false;  // Use points instead of mesh
    data.dimension = 3;

    // Define arbitrary k-points
    data.points = {
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.0, 0.5, 0.0}
    };

    // Define domain (not used for interpolation when as_mesh = false, but kept for compatibility)
    data.domain = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    // Fill with complex values - one per k-point
    vector<cfloat> values = {
        cfloat(1.0, 0.1),
        cfloat(2.0, 0.2),
        cfloat(3.0, 0.3),
        cfloat(4.0, 0.4)
    };
    data.data = values;

    // Verify nk() works correctly
    if (data.nk() != 4) return false;

    // Save to file
    string fname = "testpoints.h5";
    save_data_to_hdf5(data, fname);

    // Reload
    BaseData loaded = load_data_from_hdf5(fname);

    // Verify metadata
    if (loaded.as_mesh != false) return false;
    if (loaded.with_k != true) return false;
    if (loaded.nk() != 4) return false;

    // Verify points
    if (loaded.points.size() != 4) return false;
    for (size_t i = 0; i < loaded.points.size(); ++i) {
        if (loaded.points[i].size() != 3) return false;
        for (size_t j = 0; j < 3; ++j) {
            if (fabs(loaded.points[i][j] - data.points[i][j]) > 1e-6) return false;
        }
    }

    // Verify data
    auto& loaded_vals = loaded.get<vector<cfloat>>();
    if (loaded_vals.size() != values.size()) return false;
    for (size_t i = 0; i < loaded_vals.size(); ++i) {
        if (abs(loaded_vals[i] - values[i]) > 1e-6f) return false;
    }

    filesystem::remove(fname);
    return true;
}

bool point_storage_with_frequency() {
    // Test point storage with frequency dependence
    BaseData data;
    data.is_complex = true;
    data.is_vector = false;
    data.with_k = true;
    data.with_w = true;
    data.as_mesh = false;  // Use points instead of mesh
    data.dimension = 2;

    // Define arbitrary k-points
    data.points = {
        {0.0, 0.0},
        {0.5, 0.5},
        {1.0, 0.0}
    };

    // Define frequency points
    data.w_points = {1.0, 2.0, 3.0};

    data.domain = {{1.0, 0.0}, {0.0, 1.0}};

    // Fill with complex values - nk * nw = 3 * 3 = 9 values (w-k ordering)
    vector<cfloat> values;
    for (int iw = 0; iw < 3; ++iw) {
        for (int ik = 0; ik < 3; ++ik) {
            float val = (ik + 1) * (iw + 1);
            values.push_back(cfloat(val, val / 10));
        }
    }
    data.data = values;

    // Verify nk() and nw() work correctly
    if (data.nk() != 3) return false;
    if (data.nw() != 3) return false;

    // Save to file
    string fname = "testpoints_freq.h5";
    save_data_to_hdf5(data, fname);

    // Reload
    BaseData loaded = load_data_from_hdf5(fname);

    // Verify metadata
    if (loaded.as_mesh != false) return false;
    if (loaded.with_k != true) return false;
    if (loaded.with_w != true) return false;
    if (loaded.nk() != 3) return false;
    if (loaded.nw() != 3) return false;

    // Verify points
    if (loaded.points.size() != 3) return false;
    for (size_t i = 0; i < loaded.points.size(); ++i) {
        if (loaded.points[i].size() != 2) return false;
        for (size_t j = 0; j < 2; ++j) {
            if (fabs(loaded.points[i][j] - data.points[i][j]) > 1e-6) return false;
        }
    }

    // Verify w_points
    if (loaded.w_points.size() != 3) return false;
    for (size_t i = 0; i < loaded.w_points.size(); ++i) {
        if (fabs(loaded.w_points[i] - data.w_points[i]) > 1e-6) return false;
    }

    // Verify data
    auto& loaded_vals = loaded.get<vector<cfloat>>();
    if (loaded_vals.size() != values.size()) return false;
    for (size_t i = 0; i < loaded_vals.size(); ++i) {
        if (abs(loaded_vals[i] - values[i]) > 1e-6f) return false;
    }

    filesystem::remove(fname);
    return true;
}

bool evaluate_tensor_complex_1d_w() {
    BaseData data;
    data.is_complex = true;
    data.is_vector = false;
    data.with_k = false;
    data.with_w = true;
    data.as_mesh = true;
    data.mesh = {mpts};
    data.dimension = 1;
    data.inds = {2, 2};  // 2x2 matrix

    data.domain = {
        {1.0}
    };
    data.w_points = {1.0, 2.0, 3.0};
    data.data = create_data(data.dimension, mpts, data.w_points);

    DataEvaluator field(data);
    Vec p(0.25, 0.25, 0.25);
    auto result = field(p, 1.5);
    if (auto* s = std::get_if<cfloat>(&result)) {
        return fabs(*s - cfloat(2.25, 0.225)) < 1e-6;
    }
    return false;
}

bool base_data_tests() {
    int num_tests = 15;
    bool all_tests[num_tests] = {
        create_destroy(),
        evaluate_real_scalar_1d_w(),
        evaluate_complex_scalar_1d_w(),
        evaluate_real_scalar_1d_k(),
        evaluate_complex_scalar_1d_k(),
        evaluate_real_scalar_2d_k(),
        evaluate_complex_scalar_2d_k(),
        evaluate_real_scalar_2d_w(),
        evaluate_complex_scalar_2d_w(),
        evaluate_real_scalar_3d_k(),
        evaluate_complex_scalar_3d_k(),
        evaluate_real_scalar_3d_w(),
        evaluate_complex_scalar_3d_w(),
        point_storage_scalar(),
        point_storage_with_frequency(),
    };
    filesystem::remove("testfield.h5");
    return print_test_results(all_tests, num_tests, "BaseData tests");
}

