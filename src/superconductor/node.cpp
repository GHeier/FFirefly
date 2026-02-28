#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void superconductor_wrapper() {
    printv("Running superconductor_wrapper\n");
    if (calculation == "bcs" && method == "convolution") {
        if (debug) run_python_test("convolution");
        else run_python_method("convolution");
    }
    else if (calculation == "bcs" && method == "matrix") {
        if (debug) run_cpp_test("matrix");
        else run_cpp_method("matrix");
    }
    else if (calculation == "bcs_w" && method == "hmatrix") {
        if (debug) run_julia_test("hmatrix");
        else run_julia_method("hmatrix");
    }
    else if (calculation == "eliashberg" && method == "convolution") {
        if (debug) run_python_test("convolution");
        else run_python_method("convolution");
    }
    else if (calculation == "eliashberg" && method == "hmatrix") {
        if (debug) run_julia_test("hmatrix");
        else run_julia_method("hmatrix");
    }
    else if (calculation == "eliashberg" && method == "sparse_ir") {
        if (debug) run_julia_test("sparse_ir");
        else run_julia_method("sparse_ir");
    }
    else {
        printf("In superconductor category, calculation `%s` with method `%s` not recognized\n", calculation.c_str(), method.c_str());
    }
}