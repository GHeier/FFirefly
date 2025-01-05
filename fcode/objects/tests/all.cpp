#include <vector>

#include "all.hpp"
#include "../../config/load/c_config.h"
#include "../field.hpp"

using namespace std;

extern "C" bool object_tests() {
    printf("Running Object tests\n");
    int num_tests = 1;
    bool all_tests[num_tests] = {
        python_field_test()
    };
    return print_test_results(all_tests, num_tests, "Object tests");
}


bool python_field_test() {
    vector<Vec> points = {Vec(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0), Vec(0.0, 1.0, 0.0), Vec(1.0, 1.0, 0.0), Vec(0.0, 0.0, 1.0), Vec(1.0, 0.0, 1.0), Vec(0.0, 1.0, 1.0), Vec(1.0, 1.0, 1.0)};
    vector<float> values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    ScalarField test_var(points, values, 3);
    vector<double> coords = {0.5, 0.5, 0.5};
    double result = test_var(coords);
    return fabs(result - 4.5) < 0.0001;
}
