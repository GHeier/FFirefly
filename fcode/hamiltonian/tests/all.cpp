#include "all.hpp"
#include "../../superconductor/cfg.hpp"
#include "../../config/load/c_config.h"
#include "../../config/load/cpp_config.hpp"
#include "../../hamiltonian/band_structure.hpp"
#include "../../superconductor/solver.hpp"
#include <cmath>

using namespace std;

extern "C" void load_cpp_config_wrapper() {
    load_cpp_config();
}

extern "C" bool hamiltonian_tests() {
    printf("Running Hamiltonian tests\n");
    load_cpp_cfg();
    int num_tests = 1;
    bool all_tests[num_tests] = {
        DOS_test()
    };
    return print_test_results(all_tests, num_tests, "Hamiltonian tests");
}

bool DOS_test() {
    // Set global variables for testing
    set_global(nbnd, 1);
    band.push_back("fermi_gas");
    eff_mass.push_back(0.5);
    set_global(fermi_energy, 1.0);
    set_global(k_mesh, {80, 80, 80});


    float answer = 0.02533;
    printv("Starting Density of States Test\n");

    vector<Vec> FS = get_FS(fermi_energy);
    float DOS = get_DOS(FS);

    if (fabs(answer - DOS) < 0.001) {
        printv("Density of States Test Passed\n");
        return true;
    }
    else {
        printf("Density of States Test Failed\nCalculated: %.3f\nExpected: %.3f\n", DOS, answer);
        return false;
    }
}


