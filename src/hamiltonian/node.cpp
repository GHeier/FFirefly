#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void hamiltonian_wrapper() {
    printv("Running hamiltonian_wrapper\n");
    if (calculation == "DOS" && method == "gaussian") {
        if (debug) run_python_test("gaussian");
        else run_python_method("gaussian");
    }
    else if (calculation == "DOS" && method == "tetrahedra") {
        if (debug) run_cpp_test("tetrahedra");
        else run_cpp_method("tetrahedra");
    }
    else if (calculation == "FS" && method == "tetrahedra") {
        if (debug) run_cpp_test("tetrahedra");
        else run_cpp_method("tetrahedra");
    }
    else if (calculation == "generate" && method == "hk_from_hr") {
        if (debug) run_python_test("hk_from_hr");
        else run_python_method("hk_from_hr");
    }
    else if (calculation == "generate" && method == "band_structure") {
        if (debug) run_python_test("band_structure");
        else run_python_method("band_structure");
    }
    else {
        printf("In hamiltonian category, calculation `%s` with method `%s` not recognized\n", calculation.c_str(), method.c_str());
    }
}