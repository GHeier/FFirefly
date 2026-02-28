#include "run.hpp"
#include "../../../config/load/cpp_config.hpp"
#include "../../../config/load/c_config.h"
#include "../../../objects/CMField/fields.hpp"
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "renormalization.hpp"



float run() {
    // Main function call goes here
    float m_star = 0.0;
    if (interaction == "FLEX") {
        m_star = FLEX_renormalization();
    }
    else {
        std::cerr << "Error: Unsupported interaction type: " << interaction << "\n";
        exit(1);
    }

    return m_star;
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



