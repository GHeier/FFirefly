#include "../c_config.h"
#include "jl_interface_typed_tests.hpp"
#include "../jl_interface.h"
#include <iostream>
#include <cstring>
#include <cmath>

using namespace std;

// Integer tests
static bool test_jl_int_42() {
    int result = call_julia_func_int("config/", "test_typed_interface", "TestTypedInterface", "return_int_42");
    if (result != 42) {
        return false;
    }
    return true;
}

static bool test_jl_int_negative() {
    int result = call_julia_func_int("config/", "test_typed_interface", "TestTypedInterface", "return_int_negative");
    if (result != -100) {
        return false;
    }
    return true;
}

static bool test_jl_int_zero() {
    int result = call_julia_func_int("config/", "test_typed_interface", "TestTypedInterface", "return_int_zero");
    if (result != 0) {
        return false;
    }
    return true;
}

static bool test_jl_compute_sum() {
    int result = call_julia_func_int("config/", "test_typed_interface", "TestTypedInterface", "compute_sum");
    if (result != 60) {
        return false;
    }
    return true;
}

// Float tests
static bool test_jl_float_pi() {
    float result = call_julia_func_float("config/", "test_typed_interface", "TestTypedInterface", "return_float_pi");
    if (fabs(result - 3.14159f) > 0.0001f) {
        return false;
    }
    return true;
}

static bool test_jl_float_negative() {
    float result = call_julia_func_float("config/", "test_typed_interface", "TestTypedInterface", "return_float_negative");
    if (fabs(result - (-2.71828f)) > 0.0001f) {
        return false;
    }
    return true;
}

static bool test_jl_float_zero() {
    float result = call_julia_func_float("config/", "test_typed_interface", "TestTypedInterface", "return_float_zero");
    if (fabs(result - 0.0f) > 0.0001f) {
        return false;
    }
    return true;
}

static bool test_jl_compute_product() {
    float result = call_julia_func_float("config/", "test_typed_interface", "TestTypedInterface", "compute_product");
    if (fabs(result - 10.0f) > 0.0001f) {
        return false;
    }
    return true;
}

// Double tests
static bool test_jl_double_large() {
    double result = call_julia_func_double("config/", "test_typed_interface", "TestTypedInterface", "return_double_large");
    if (fabs(result - 1.23456789012345) > 0.000000001) {
        return false;
    }
    return true;
}

static bool test_jl_double_small() {
    double result = call_julia_func_double("config/", "test_typed_interface", "TestTypedInterface", "return_double_small");
    if (fabs(result - 0.00000123456) > 0.0000000001) {
        return false;
    }
    return true;
}

// String tests
static bool test_jl_string_hello() {
    const char* result = call_julia_func_string("config/", "test_typed_interface", "TestTypedInterface", "return_string_hello");
    bool success = (strcmp(result, "Hello, World!") == 0);
    free((void*)result);
    return success;
}

static bool test_jl_string_empty() {
    const char* result = call_julia_func_string("config/", "test_typed_interface", "TestTypedInterface", "return_string_empty");
    bool success = (strcmp(result, "") == 0);
    free((void*)result);
    return success;
}

static bool test_jl_string_special() {
    const char* result = call_julia_func_string("config/", "test_typed_interface", "TestTypedInterface", "return_string_special");
    bool success = (strcmp(result, "Test!@#$%^&*()") == 0);
    free((void*)result);
    return success;
}

static bool test_jl_compute_message() {
    const char* result = call_julia_func_string("config/", "test_typed_interface", "TestTypedInterface", "compute_message");
    bool success = (strcmp(result, "FFirefly v1.0") == 0);
    free((void*)result);
    return success;
}

// Boolean tests
static bool test_jl_bool_true() {
    bool result = call_julia_func_bool("config/", "test_typed_interface", "TestTypedInterface", "return_bool_true");
    if (result != true) {
        return false;
    }
    return true;
}

static bool test_jl_bool_false() {
    bool result = call_julia_func_bool("config/", "test_typed_interface", "TestTypedInterface", "return_bool_false");
    if (result != false) {
        return false;
    }
    return true;
}

bool jl_interface_typed_tests() {

    int num_tests = 16;
    bool all_tests[] = {
        test_jl_int_42(),
        test_jl_int_negative(),
        test_jl_int_zero(),
        test_jl_compute_sum(),
        test_jl_float_pi(),
        test_jl_float_negative(),
        test_jl_float_zero(),
        test_jl_compute_product(),
        test_jl_double_large(),
        test_jl_double_small(),
        test_jl_string_hello(),
        test_jl_string_empty(),
        test_jl_string_special(),
        test_jl_compute_message(),
        test_jl_bool_true(),
        test_jl_bool_false(),
    };

    return print_test_results(all_tests, num_tests, "Julia Typed Interface tests");
}
