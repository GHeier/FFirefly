#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void hamiltonian_wrapper() {
    printv("Running hamiltonian_wrapper\n");
    if (calculation == "DOS" && method == "gaussian") run_python_method("gaussian");
    else if (calculation == "DOS" && method == "tetrahedra") run_cpp_method("tetrahedra");
    else if (calculation == "FS" && method == "tetrahedra") run_cpp_method("tetrahedra");
    else if (calculation == "generate" && method == "hk_from_hr") run_python_method("hk_from_hr");
    else {
        printf("In hamiltonian category, calculation `%s` with method `%s` not recognized\n", calculation.c_str(), method.c_str());
    }
}