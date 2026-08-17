#include "src/config/load/c_config.h"
#include "py_interface_typed_tests.hpp"
#include "src/config/load/py_interface.h"
#include <iostream>
#include <cstring>
#include <cmath>

using namespace std;

// Integer tests
static bool test_int_42() {
    int result = call_python_func_int("config/tests", "test_typed_interface", "return_int_42");
    if (result != 42) {
        return false;
    }
    return true;
}

static bool test_int_negative() {
    int result = call_python_func_int("config/tests", "test_typed_interface", "return_int_negative");
    if (result != -100) {
        return false;
    }
    return true;
}

static bool test_int_zero() {
    int result = call_python_func_int("config/tests", "test_typed_interface", "return_int_zero");
    if (result != 0) {
        return false;
    }
    return true;
}

static bool test_compute_sum() {
    int result = call_python_func_int("config/tests", "test_typed_interface", "compute_sum");
    if (result != 60) {
        return false;
    }
    return true;
}

// Float tests
static bool test_float_pi() {
    float result = call_python_func_float("config/tests", "test_typed_interface", "return_float_pi");
    if (abs(result - 3.14159f) > 0.0001f) {
        return false;
    }
    return true;
}

static bool test_float_negative() {
    float result = call_python_func_float("config/tests", "test_typed_interface", "return_float_negative");
    if (abs(result - (-2.71828f)) > 0.0001f) {
        return false;
    }
    return true;
}

static bool test_float_zero() {
    float result = call_python_func_float("config/tests", "test_typed_interface", "return_float_zero");
    if (abs(result - 0.0f) > 0.0001f) {
        return false;
    }
    return true;
}

static bool test_compute_product() {
    float result = call_python_func_float("config/tests", "test_typed_interface", "compute_product");
    if (abs(result - 10.0f) > 0.0001f) {
        return false;
    }
    return true;
}

// Double tests
static bool test_double_large() {
    double result = call_python_func_double("config/tests", "test_typed_interface", "return_double_large");
    if (abs(result - 1.23456789012345) > 0.000000001) {
        return false;
    }
    return true;
}

static bool test_double_small() {
    double result = call_python_func_double("config/tests", "test_typed_interface", "return_double_small");
    if (abs(result - 0.00000123456) > 0.0000000001) {
        return false;
    }
    return true;
}

// String tests
static bool test_string_hello() {
    const char* result = call_python_func_string("config/tests", "test_typed_interface", "return_string_hello");
    bool success = (strcmp(result, "Hello, World!") == 0);
    free((void*)result);
    return success;
}

static bool test_string_empty() {
    const char* result = call_python_func_string("config/tests", "test_typed_interface", "return_string_empty");
    bool success = (strcmp(result, "") == 0);
    free((void*)result);
    return success;
}

static bool test_string_special() {
    const char* result = call_python_func_string("config/tests", "test_typed_interface", "return_string_special");
    bool success = (strcmp(result, "Test!@#$%^&*()") == 0);
    free((void*)result);
    return success;
}

static bool test_compute_message() {
    const char* result = call_python_func_string("config/tests", "test_typed_interface", "compute_message");
    bool success = (strcmp(result, "FFirefly v1.0") == 0);
    free((void*)result);
    return success;
}

bool py_interface_typed_tests() {

    int num_tests = 14;
    bool all_tests[] = {
        test_int_42(),
        test_int_negative(),
        test_int_zero(),
        test_compute_sum(),
        test_float_pi(),
        test_float_negative(),
        test_float_zero(),
        test_compute_product(),
        test_double_large(),
        test_double_small(),
        test_string_hello(),
        test_string_empty(),
        test_string_special(),
        test_compute_message(),
    };

    return print_test_results(all_tests, num_tests, "Python Typed Interface tests");
}
