#include "all.hpp"
#include "module_interface_tests.hpp"
#include <iostream>

extern "C" bool module_tests() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Running Module Interface Tests" << std::endl;
    std::cout << "========================================\n" << std::endl;

    bool py_passed = py_module_interface_tests();
    bool jl_passed = jl_module_interface_tests();

    return py_passed && jl_passed;
}
