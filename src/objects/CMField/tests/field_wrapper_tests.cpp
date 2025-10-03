#include "../field.hpp"
#include "../../../config/load/c_config.h"
#include <iostream>
#include <cmath>

using namespace std;

bool field_wrapper_basic_test() {
    // Create simple 1D k-space data
    // Data represents points at x=0, 0.5, 1.0 in DataEvaluator
    // Field centers this to x=-0.5, 0, 0.5
    vector<cfloat> data = {
        cfloat(0.0, 0.0),
        cfloat(0.5, 0.05),
        cfloat(1.0, 0.1)
    };

    vector<int> mesh = {3};
    vector<vector<float>> domain = {{1.0}};

    Field field(data, false, false, mesh, domain);

    // Query at x=0 in the centered coordinate system
    // This corresponds to x=0.5 in the original [0,1] system
    Vec v(0.0);
    auto result = field(v);

    if (auto* s = std::get_if<float>(&result)) {
        return fabs(*s - 0.5) < 1e-6;
    }
    return false;
}

bool field_wrapper_with_w_test() {
    // Create 2D data with frequency
    vector<cfloat> data;
    int mpts = 3;
    vector<float> w_points = {1.0, 2.0, 3.0};

    // Generate data: for each spatial point, then for each w
    // Data at points (0,0), (0.5,0), (1,0), (0,0.5), etc. in [0,1]x[0,1]
    for (int i = 0; i < mpts; i++) {
        for (int j = 0; j < mpts; j++) {
            for (int w = 0; w < 3; w++) {
                float x = 1.0 * i / (mpts - 1);
                float y = 1.0 * j / (mpts - 1);
                float w_val = w_points[w];
                float val = x + y + w_val;
                data.push_back(cfloat(val, val / 10));
            }
        }
    }

    vector<int> mesh = {mpts, mpts};
    vector<vector<float>> domain = {{1.0, 0.0}, {0.0, 1.0}};

    Field field(data, true, false, mesh, domain, w_points);

    // Query at p=(0, 0) in centered coordinates
    // This maps to (0.5, 0.5) in the original [0,1]x[0,1] system
    Vec p(0.0, 0.0);
    auto result = field(p, 1.5);

    if (auto* s = std::get_if<cfloat>(&result)) {
        // At original p=(0.5, 0.5), spatial contribution = 0.5 + 0.5 = 1.0
        // At w=1.5, w contribution = 1.5
        // Total = 2.5
        return fabs(*s - cfloat(2.5, 0.25)) < 1e-6;
    }
    return false;
}

bool field_wrapper_tests() {
    int num_tests = 2;
    bool all_tests[num_tests] = {
        field_wrapper_basic_test(),
        field_wrapper_with_w_test(),
    };
    return print_test_results(all_tests, num_tests, "Field Wrapper tests");
}
