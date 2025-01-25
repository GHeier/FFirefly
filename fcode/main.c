#include <stdio.h>
#include <omp.h>

#include "config/load/c_config.h"
#include "config/load/py_interface.h"

// Category nodes below
#include "response/DOS.h"
#include "response/response.h"
#include "superconductor/superconductor.hpp"

// Test nodes below
#include "hamiltonian/tests/all.hpp"
#include "objects/tests/all.hpp"

// Global category calls
void DOS() {
    printf("Starting DOS Calculation\n\n");
    DOS_spectrum();
}

void response() {
    printf("Starting Response Calculation\n\n");
    response_wrapper();
}

void superconductor() {
    printf("Starting Superconductor Calculation\n\n");
    superconductor_wrapper();
}

// Global test calls
void test() {
    printf("Starting Test Calculations\n\n");

    int num_tests = 2;

    bool all_tests[num_tests];
    all_tests[0] = hamiltonian_tests();
    all_tests[1] = object_tests();

    print_test_results(all_tests, num_tests, "tests");
}

int main() {
    printf("Starting Program\n");
    // Set the number of threads used in parallelization to one less than the maximum
    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs - 1);
    printf("Number of threads used in CPU parallelization: %d\n", num_procs - 1);

    load_c_config(); // Read input to load global c variables
    load_cpp_config_wrapper(); // Read input to load global cpp variables
    start_python();


    if (!strcmp(c_category, "DOS")) DOS();
    else if (!strcmp(c_category, "response")) response();
    else if (!strcmp(c_category, "superconductor")) superconductor();
    else if (!strcmp(c_category, "test")) test();
    else printf("Unknown Calculation Type\n\n");

    //unload_c_config(); // Free memory allocated for global c variables
    end_python();
    printf("Program Complete\n");
    return 0;
}
