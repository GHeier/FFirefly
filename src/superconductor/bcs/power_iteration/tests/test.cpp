#include <cmath>
#include <vector>

using namespace std;

#include "src/superconductor/bcs/power_iteration/run.hpp"

// Set configuration variables for test run
vector<int> k_mesh = {4, 4, 4}; // Example k-mesh values

int test() { // Runs on file execution
    printf("Welcome to Testing! This is a Firefly run with k-mesh: [%d %d %d]\n", k_mesh[0], k_mesh[1], k_mesh[2]);
    float result = run();
    return fabs(result - 3.14) < 1e-6; // Example test condition
}     




