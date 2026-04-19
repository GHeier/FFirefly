#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void many_body_wrapper() {
    printv("Running many_body_wrapper\n");
    if (calculation == "many_body" && method == "triqs") {
        if (debug) run_python_test("triqs");
        else run_python_method("triqs");
    }
    else if (calculation == "many_body" && method == "sparse_ir") {
        if (debug) run_julia_test("sparse_ir");
        else run_julia_method("sparse_ir");
    }
    else if (calculation == "self_energy" && method == "sparse_ir") {
        if (debug) run_julia_test("sparse_ir");
        else run_julia_method("sparse_ir");
    }
    else if (calculation == "self_energy" && method == "triqs") {
        if (debug) run_python_test("triqs");
        else run_python_method("triqs");
    }
    else if (calculation == "vertex" && method == "from_susceptibility") {
        if (debug) run_cpp_test("from_susceptibility");
        else run_cpp_method("from_susceptibility");
    }
    else if (calculation == "response" && method == "tetrahedra") {
        if (debug) run_julia_test("tetrahedra");
        else run_julia_method("tetrahedra");
    }
    else if (calculation == "response" && method == "sparse_ir") {
        if (debug) run_julia_test("sparse_ir");
        else run_julia_method("sparse_ir");
    }
    else if (calculation == "renormalization" && method == "analytic") {
        if (debug) run_cpp_test("analytic");
        else run_cpp_method("analytic");
    }
    else if (calculation == "renormalization" && method == "from_sigma") {
        if (debug) run_cpp_test("from_sigma");
        else run_cpp_method("from_sigma");
    }
    else if (calculation == "renormalization" && method == "FS_approx") {
        if (debug) run_python_test("FS_approx");
        else run_python_method("FS_approx");
    }
    else {
        printf("In many_body category, calculation `%s` with method `%s` not recognized\n", calculation.c_str(), method.c_str());
    }
}