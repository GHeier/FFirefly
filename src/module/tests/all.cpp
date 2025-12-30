#include "all.hpp"
// #include "module_interface_tests.hpp" // Deprecated - uses embedded Python/Julia
#include <iostream>

extern "C" bool module_tests() {
    printf("\nRunning Module tests\n");
    printf("Note: Python/Julia module interface tests disabled (embedded interpreters removed)\n");
    printf("      Use Python/Julia scripts directly via run_python_method2/run_julia_method2\n");

    // Embedded interpreter tests are deprecated
    // bool py_passed = py_module_interface_tests();
    // bool jl_passed = jl_module_interface_tests();

    // Return true since there are no active tests
    return true;
}
