#include "run.hpp"
#include "../../../config/load/cpp_config.hpp"
#include "../../../config/load/c_config.h"
#include "../../../objects/CMField/fields.hpp"
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "vertex.hpp"


float run() {
    // Main function call goes here
    call_flex();

    return 3.14; // Return something of any type that can be tested in the test suite.
}

int main() {
    // Load configuration from the build directory
    std::string config_path = get_loc() + "input.cfg";

    // Load configuration using existing infrastructure
    read_c_config(config_path.c_str());
    load_cpp_config();

    // Run the method
    float result = run();

    // Success
    return 0;
}



