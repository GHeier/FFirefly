#include <stdio.h>
#include <omp.h>

#include "config/load/c_config.h"
#include "superconductor/superconductor.hpp"

extern void polarization_wrapper();
void response() {
    polarization_wrapper();
}

void superconductor() {
    printf("Starting Superconductor Calculation\n");
    superconductor_wrapper();
}

int main(int argc, char *argv[]) {
    printf("Starting Program\n");
    // Set the number of threads used in parallelization to one less than the maximum
    // This allows for the last thread to be used for other tasks
    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs - 1);
    printf("Number of threads used in CPU parallelization: %d\n", num_procs - 1);

    load_c_config(); // Read input to load global c variables

    if (!strcmp(c_category, "response")) response();
    else if (!strcmp(c_category, "superconductor")) superconductor();
    else printf("Unknown Calculation Type\n");

    //unload_c_config(); // Free memory allocated for global c variables
    printf("Program Complete\n");
    return 0;
}
