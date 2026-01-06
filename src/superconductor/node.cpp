#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void superconductor_wrapper() {
    printv("Running superconductor_wrapper\n");
    if (calculation == "bcs" && method == "convolution") run_python_method("convolution");
    else if (calculation == "bcs" && method == "matrix") run_cpp_method("matrix");
    else if (calculation == "eliashberg" && method == "convolution") run_python_method("convolution");
    else if (calculation == "eliashberg" && method == "hmatrix") run_julia_method("hmatrix");
    else {
        printf("In superconductor category, calculation `%s` with method `%s` not recognized\n", calculation.c_str(), method.c_str());
    }
}