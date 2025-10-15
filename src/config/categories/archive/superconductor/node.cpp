#include "node.hpp"
#include "../config/load/cpp_config.hpp"
#include "superconductor.hpp"
#include <iostream>

using namespace std;


/**
 * Wrapper function for superconductor category
 * Dispatches to appropriate calculation/method based on config variables
 */
extern "C" void superconductor_wrapper() {
    printv("Running superconductor_wrapper\n");
    if (calculation == "bcs") {
        bcs();
    }
    else if (calculation == "debug") {
        debug();
    }
    else if (calculation == "eliashberg") {
        eliashberg();
    }
    else if (calculation == "linearized_eliashberg") {
        linearized_eliashberg();
    }
    else
        cout << "calculation \"" << calculation << "\" not recognized for category superconductor" << endl;
}
