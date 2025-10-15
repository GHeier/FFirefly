#include "../c_config.h"
#include "py_interface_bool_tests.hpp"
#include "../py_interface.h"
#include <iostream>
#include <cassert>

using namespace std;

static bool test_return_true() {
    bool result = call_python_func_bool("config/tests", "test_bool_interface", "return_true");
    if (result != true) {
        return false;
    }
    return true;
}

static bool test_return_false() {
    bool result = call_python_func_bool("config/tests", "test_bool_interface", "return_false");
    if (result != false) {
        return false;
    }
    return true;
}

static bool test_return_truthy_int() {
    bool result = call_python_func_bool("config/tests", "test_bool_interface", "return_truthy_int");
    if (result != true) {
        return false;
    }
    return true;
}

static bool test_return_falsy_int() {
    bool result = call_python_func_bool("config/tests", "test_bool_interface", "return_falsy_int");
    if (result != false) {
        return false;
    }
    return true;
}

static bool test_complex_check() {
    bool result = call_python_func_bool("config/tests", "test_bool_interface", "complex_check");
    if (result != true) {
        return false;
    }
    return true;
}

static bool test_another_complex_check() {
    bool result = call_python_func_bool("config/tests", "test_bool_interface", "another_complex_check");
    if (result != false) {
        return false;
    }
    return true;
}

bool py_interface_bool_tests() {

    int num_tests = 6;
    bool all_tests[] = {
        test_return_true(),
        test_return_false(),
        test_return_truthy_int(),
        test_return_falsy_int(),
        test_complex_check(),
        test_another_complex_check(),
    };

    return print_test_results(all_tests, num_tests, "Python Interface tests");
}
