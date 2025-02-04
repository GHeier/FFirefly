#include <string>
#include <vector>

#include "all.hpp"
#include "../../config/load/c_config.h"
#include "../../config/load/jl_interface.h"

using namespace std;

extern "C" bool object_tests() {
    printf("Running Object tests\n");
    int num_tests = 1;
    bool all_tests[num_tests] = {
        field_tests()
    };
    return print_test_results(all_tests, num_tests, "Object tests");
}

bool field_tests() {
    string folder = "objects/";
    string file = "field";
    string module = "CondensedMatterField";
    string function = "test_2dim_real_scalar_field";
    call_julia_func(folder.c_str(), file.c_str(), module.c_str(), function.c_str());
    function = "test_2dim_complex_scalar_field";
    call_julia_func(folder.c_str(), file.c_str(), module.c_str(), function.c_str());
    function = "test_2dim_complex_scalar_field_with_w";
    call_julia_func(folder.c_str(), file.c_str(), module.c_str(), function.c_str());
    return false;
}

