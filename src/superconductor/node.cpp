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
    if (method == "bcs")
        bcs();
    else if (method == "eliashberg")
        eliashberg();
    else if (method == "linearized_eliashberg")
        linearized_eliashberg();
    else if (method == "debug")
        debug();
    else
        cout << "Method " << method << " not recognized" << endl;
}
