#include "src/config/load/cpp_config.hpp"
#include <iostream>

int main() {
    // Example 1: Run fly.x with a config file
    std::cout << "Example usage of run_with_config:\n";
    std::cout << "===================================\n\n";

    // This would run: build/bin/fly.x < input.cfg
    // int result = run_with_config("build/bin/fly.x", "input.cfg");

    std::cout << "Usage:\n";
    std::cout << "  int exitcode = run_with_config(\"./executable.exe\", \"input.cfg\");\n\n";

    std::cout << "This function:\n";
    std::cout << "  - Verifies executable and config file exist\n";
    std::cout << "  - Runs: executable < input.cfg\n";
    std::cout << "  - Returns the exit code of the executable\n";
    std::cout << "  - Returns -1 on error\n\n";

    std::cout << "The function is declared in: src/config/load/cpp_config.hpp\n";
    std::cout << "Implementation in: src/config/load/cpp_config.cpp\n";

    return 0;
}
