#include "all.hpp"
#include "py_interface_bool_tests.hpp"
#include "py_interface_typed_tests.hpp"
#include "jl_interface_typed_tests.hpp"
#include <iostream>

using namespace std;

extern "C" bool config_load_tests() {
    cout << "\nRunning Config Load tests" << endl;

    int num_tests = 3;
    bool all_tests[num_tests] = {
        py_interface_bool_tests(),
        py_interface_typed_tests(),
        jl_interface_typed_tests(),
    };

    int passed = 0;
    for (int i = 0; i < num_tests; i++) {
        if (all_tests[i]) passed++;
    }

    return (passed == num_tests);
}
