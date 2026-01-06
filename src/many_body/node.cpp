#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void many_body_wrapper() {
    printv("Running many_body_wrapper\n");
    if (calculation == "many_body" && method == "triqs") run_python_method("triqs");
    else if (calculation == "many_body" && method == "sparse_ir") run_julia_method("sparse_ir");
    else if (calculation == "self_energy" && method == "sparse_ir") run_julia_method("sparse_ir");
    else if (calculation == "self_energy" && method == "triqs") run_python_method("triqs");
    else if (calculation == "vertex" && method == "from_susceptibility") run_cpp_method("from_susceptibility");
    else if (calculation == "response" && method == "bz_integral") run_julia_method("bz_integral");
    else if (calculation == "response" && method == "sparse_ir") run_julia_method("sparse_ir");
    else if (calculation == "renormalization" && method == "analytic") run_cpp_method("analytic");
    else if (calculation == "renormalization" && method == "from_sigma") run_cpp_method("from_sigma");
    else {
        printf("In many_body category, calculation `%s` with method `%s` not recognized\n", calculation.c_str(), method.c_str());
    }
}