#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void superconductor_wrapper() {
    printv("Running superconductor_wrapper\n");
    if (calculation == "bcs" && method == "lanczos") run_cpp_method("lanczos");
    else if (calculation == "bcs" && method == "power_iteration") run_cpp_method("power_iteration");
    else if (calculation == "eliashberg" && method == "lanczos") run_python_method("lanczos");
    else if (calculation == "eliashberg" && method == "power_iteration") run_julia_method("power_iteration");
    else {
        printf("In superconductor category, calculation `%s` with method `%s` not recognized\n", calculation.c_str(), method.c_str());
    }
}