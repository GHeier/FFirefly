#include "src/objects/CMField/fields.hpp"
#include <filesystem>
#include <iostream>

#include "src/hamiltonian/models/band_structure.hpp"
#include "src/objects/CMField/hamiltonian.hpp"
#include "src/config/load/c_config.h"

using namespace std;

bool epsilon_compare() {
    Hamiltonian H_obj;

    // Test k-point
    Vec k({M_PI / 2, M_PI / 2, M_PI / 2});

    // Compute epsilon from Hamiltonian object
    auto H_matrix = H_obj(k);
    float epsilon_from_H = real(H_matrix[0][0]);

    // Compute epsilon from direct function
    float epsilon_direct = epsilon(1, k);

    // Compare results
    return abs(epsilon_from_H - epsilon_direct) < 1e-6;
}

bool hamiltonian_tests() {
    int num_tests = 1;
    bool all_tests[num_tests] = {
        epsilon_compare(),
    };
    return print_test_results(all_tests, num_tests, "Hamiltonian Object tests");
}

