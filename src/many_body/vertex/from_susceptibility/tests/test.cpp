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
    return abs(result - 3.14) < 1e-6; // Example test condition
}     




