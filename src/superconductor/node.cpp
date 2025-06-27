#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include "superconductor.hpp"
#include <iostream>

using namespace std;

/**
 * Example wrapper function below
 * This is how to connect your code to main.c
 * main.c will import {categoryname}_wrapper and run it
 * Make sure to handle incorrect input values
 */
extern "C" void superconductor_wrapper() {
    printv("Running superconductor_wrapper\n");
    if (calculation == "bcs")
        bcs();
    else if (calculation == "eliashberg")
        eliashberg();
    else if (calculation == "linearized_eliashberg")
        linearized_eliashberg();
    else if (calculation == "debug")
        debug();
    else
        cout << "calculation " << calculation << " not recognized" << endl;
}
