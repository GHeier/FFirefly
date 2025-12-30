#include "run.hpp"
#include "../../../config/load/cpp_config.hpp"
#include "../../../objects/CMField/fields.hpp"



float run() {
    # Main function call goes here
    printf("Hello, World! This is a Firefly run with k-mesh: [%d %d %d]\n", k_mesh[0], k_mesh[1], k_mesh[2]);

    return 3.14; // Return something of any type that can be tested in the test suite.
}

int main() { # Runs on file execution
    run();
    return 0;
}     



