#include "../../config/load/c_config.h"
#include "../../config/load/py_interface.h"
#include "../../config/load/jl_interface.h"
#include "module_interface_tests.hpp"
#include <iostream>
#include <cstring>

using namespace std;

// Python tests
static bool test_py_vec_constructor_empty() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_vec_constructor_empty");
}

static bool test_py_vec_constructor_args() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_vec_constructor_args");
}

static bool test_py_vec_getters() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_vec_getters");
}

static bool test_py_load_config() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_load_config");
}

static bool test_py_field_r_constructor() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_r_constructor");
}

static bool test_py_field_r_call_w() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_r_call_w");
}

static bool test_py_field_r_call_k() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_r_call_k");
}

static bool test_py_field_r_call_list() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_r_call_list");
}

static bool test_py_field_c_constructor() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_c_constructor");
}

static bool test_py_field_c_call_w() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_c_call_w");
}

static bool test_py_field_c_call_k() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_c_call_k");
}

static bool test_py_field_c_call_list() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_c_call_list");
}

static bool test_py_field_rm_constructor() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_rm_constructor");
}

static bool test_py_field_rm_call() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_rm_call");
}

static bool test_py_field_cm_constructor() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_cm_constructor");
}

static bool test_py_field_cm_call() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_field_cm_call");
}

static bool test_py_hamiltonian_constructor() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_hamiltonian_constructor");
}

static bool test_py_hamiltonian_call() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_hamiltonian_call");
}

static bool test_py_save_data_scalar() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_save_data_scalar");
}

static bool test_py_save_data_scalar_complex() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_save_data_scalar_complex");
}

static bool test_py_save_data_vector() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_save_data_vector");
}

static bool test_py_save_data_matrix() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_save_data_matrix");
}

static bool test_py_save_read_scalar_real() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_save_read_scalar_real");
}

static bool test_py_save_read_scalar_complex() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_save_read_scalar_complex");
}

static bool test_py_save_read_with_frequency() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_save_read_with_frequency");
}

static bool test_py_save_data_dispatcher_real() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_save_data_dispatcher_real");
}

static bool test_py_save_data_dispatcher_complex() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_save_data_dispatcher_complex");
}

static bool test_py_unified_field_scalar_real() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_unified_field_scalar_real");
}

static bool test_py_unified_field_scalar_complex() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_unified_field_scalar_complex");
}

static bool test_py_unified_field_matrix_real() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_unified_field_matrix_real");
}

static bool test_py_unified_field_matrix_complex() {
    return call_python_func_bool("module/tests", "test_module_interface", "test_unified_field_matrix_complex");
}

// Julia tests
static bool test_jl_vec_constructor_empty() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_vec_constructor_empty");
}

static bool test_jl_vec_constructor_args() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_vec_constructor_args");
}

static bool test_jl_vec_getters() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_vec_getters");
}

static bool test_jl_load_config() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_load_config");
}

static bool test_jl_field_r_constructor() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_r_constructor");
}

static bool test_jl_field_r_call_w() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_r_call_w");
}

static bool test_jl_field_r_call_k() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_r_call_k");
}

static bool test_jl_field_r_call_list() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_r_call_list");
}

static bool test_jl_field_c_constructor() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_c_constructor");
}

static bool test_jl_field_c_call_w() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_c_call_w");
}

static bool test_jl_field_c_call_k() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_c_call_k");
}

static bool test_jl_field_c_call_list() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_c_call_list");
}

static bool test_jl_field_rm_constructor() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_rm_constructor");
}

static bool test_jl_field_rm_call() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_rm_call");
}

static bool test_jl_field_cm_constructor() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_cm_constructor");
}

static bool test_jl_field_cm_call() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_field_cm_call");
}

static bool test_jl_hamiltonian_constructor() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_hamiltonian_constructor");
}

static bool test_jl_hamiltonian_call() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_hamiltonian_call");
}

static bool test_jl_save_data_scalar() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_save_data_scalar");
}

static bool test_jl_save_data_scalar_complex() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_save_data_scalar_complex");
}

static bool test_jl_save_data_vector() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_save_data_vector");
}

static bool test_jl_save_data_matrix() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_save_data_matrix");
}

static bool test_jl_save_read_scalar_real() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_save_read_scalar_real");
}

static bool test_jl_save_read_scalar_complex() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_save_read_scalar_complex");
}

static bool test_jl_save_read_with_frequency() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_save_read_with_frequency");
}

static bool test_jl_save_data_dispatcher_real() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_save_data_dispatcher_real");
}

static bool test_jl_save_data_dispatcher_complex() {
    return call_julia_func_bool("module/tests/", "test_module_interface", "TestModuleInterface", "test_save_data_dispatcher_complex");
}

bool py_module_interface_tests() {
    int num_tests = 17;
    bool all_tests[] = {
        test_py_vec_constructor_empty(),
        test_py_vec_constructor_args(),
        test_py_vec_getters(),
        test_py_load_config(),
        test_py_save_data_scalar(),
        test_py_save_data_scalar_complex(),
        test_py_save_data_vector(),
        test_py_save_data_matrix(),
        test_py_save_read_scalar_real(),
        test_py_save_read_scalar_complex(),
        test_py_save_read_with_frequency(),
        test_py_save_data_dispatcher_real(),
        test_py_save_data_dispatcher_complex(),
        test_py_unified_field_scalar_real(),
        test_py_unified_field_scalar_complex(),
        test_py_unified_field_matrix_real(),
        test_py_unified_field_matrix_complex(),
    };

    return print_test_results(all_tests, num_tests, "Python Module Interface tests");
}

bool jl_module_interface_tests() {
    int num_tests = 13;
    bool all_tests[] = {
        test_jl_vec_constructor_empty(),
        test_jl_vec_constructor_args(),
        test_jl_vec_getters(),
        test_jl_load_config(),
        test_jl_save_data_scalar(),
        test_jl_save_data_scalar_complex(),
        test_jl_save_data_vector(),
        test_jl_save_data_matrix(),
        test_jl_save_read_scalar_real(),
        test_jl_save_read_scalar_complex(),
        test_jl_save_read_with_frequency(),
        test_jl_save_data_dispatcher_real(),
        test_jl_save_data_dispatcher_complex(),
    };

    return print_test_results(all_tests, num_tests, "Julia Module Interface tests");
}
