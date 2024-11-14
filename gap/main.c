#include <stdio.h>
#include <omp.h>

#include "config.h"

extern void polarization_wrapper();
void response() {
    polarization_wrapper();
}

int main(int argc, char *argv[]) {
    printf("Starting Program\n");
    // Set the number of threads used in parallelization to one less than the maximum
    // This allows for the last thread to be used for other tasks
    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs - 1);
    printf("Number of threads used in CPU parallelization: %d\n", num_procs - 1);

    load_c_config();
    printf("Calculation Type: %s\n", c_calculation_type);

    if (strcmp(c_calculation_type, "response") == 0) response();
    else printf("Unknown Calculation Type\n");

    return 0;
}
