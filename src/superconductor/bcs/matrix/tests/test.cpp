#include "../run.hpp"
#include <cstdio>
#include <iostream> 
#include <cmath>
#include <vector>

using namespace std;

// Set configuration variables for test run
vector<int> k_mesh = {4, 4, 4}; // Example k-mesh values

int test() { // Runs on file execution
    printf("Welcome to Testing! This is a Firefly run with k-mesh: [%d %d %d]\n", k_mesh[0], k_mesh[1], k_mesh[2]);
    float result = run();
    return fabs(result - 3.14) < 1e-6; // Example test condition
}     

int main() {
    // Load configuration from the build directory
    const char* config_path = "/home/g/Research/FFirefly/build/bin/input.cfg";

    // Load configuration using existing infrastructure
    //read_c_config(config_path);
    //load_cpp_config();

    // Run the method
    float result = test();

    // Success
    return 0;
}



