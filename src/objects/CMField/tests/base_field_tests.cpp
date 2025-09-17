#include "../base_field.hpp"
#include "../field_evaluator.hpp"
#include <filesystem>
#include <iostream>

#include "../../../config/load/c_config.h"

using namespace std;

bool evaluate() {
    // 1. Create a simple BaseField
    BaseField data;
    data.is_complex = false;
    data.is_vector = false;
    data.with_k = false;
    data.with_w = true;
    data.as_mesh = false;

    data.dimension = 1; // 3-component vector

    vector<float> w_points = {1, 2, 3};

    // Fill with 2x3 complex values
    vector<cfloat> vecs = {
        cfloat(1.0, 0.1), cfloat(2.0, 0.2), cfloat(3.0, 0.3)
    };

    // Store as matrix of scalars (variant type)
    data.w_points = w_points;
    data.data = vecs;
    printf("evaluating\n");
    FieldEvaluator field(data);
    printf("getting esult\n");
    auto result = field(1.5);
    printf("gor esult\n");
    if (auto* s = std::get_if<cfloat>(&result)) {
            std::cout << "Got scalar: " << *s << "\n";
            return true;
    }
    return false;
    //return fabs(answer.real() - 1.5) < 1e-6;
}

bool create_destroy() {
    // 1. Create a simple BaseField
    BaseField field;
    field.is_complex = true;
    field.is_vector = true;
    field.with_k = true;
    field.with_w = false;
    field.as_mesh = true;

    field.mesh = {2};   // 2 k-points
    field.dimension = 3; // 3-component vector

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
    save_field_to_hdf5(field, fname);

    // 3. Reload
    BaseField loaded = load_field_from_hdf5(fname);

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

bool base_field_tests() {
    int num_tests = 2;
    bool all_tests[num_tests] = {
        create_destroy(),
        evaluate(),
    };
    filesystem::remove("testfield.h5");
    return print_test_results(all_tests, num_tests, "BaseField tests");
}

