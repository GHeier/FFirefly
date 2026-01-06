#include "run.hpp"
#include "../../../config/load/cpp_config.hpp"
#include "../../../config/load/c_config.h"
#include "../../../objects/CMField/fields.hpp"
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "superconductor.hpp"

float run() {
    // Main function call goes here
    float max_eig = bcs();

    return max_eig; // Return something of any type that can be tested in the test suite.
}

int main() {
    // Load configuration from the build directory
    const char* config_path = "/home/g/Research/FFirefly/build/bin/input.cfg";

    // Load configuration using existing infrastructure
    read_c_config(config_path);
    load_cpp_config();

    // Run the method
    float result = run();

    // Success
    return 0;
}



