#include "run.hpp"
#include "../../../config/load/cpp_config.hpp"
#include "../../../config/load/c_config.h"
#include "../../../objects/CMField/fields.hpp"
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "fs.hpp"



float run() {
    // Main function call goes here
    printf("Calling Fermi Surface creation\n");
    save_FS();

    return 3.14; // Return something of any type that can be tested in the test suite.
}

int main() {
    // Read configuration from stdin (supports: executable < input.cfg)
    // First, save stdin to a temporary file since read_c_config expects a file path
    const char* tmp_config = "/tmp/firefly_method_config.cfg";

    FILE* tmp_file = fopen(tmp_config, "w");
    if (!tmp_file) {
        std::cerr << "Error: Could not create temporary config file\n";
        return 1;
    }

    // Copy stdin to temporary file
    char buffer[4096];
    while (fgets(buffer, sizeof(buffer), stdin)) {
        fputs(buffer, tmp_file);
    }
    fclose(tmp_file);

    // Load configuration using existing infrastructure
    read_c_config(tmp_config);
    load_cpp_config();

    // Remove temporary file
    remove(tmp_config);

    // Run the method
    float result = run();

    // Success
    return 0;
}



