#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "config.h"

int main(int argc, char *argv[]) {
    printf("Starting Program\n");

    if (argc < 2) {
        fprintf(stderr, "Error: No file provided\n");
        exit(1);
    }

    FILE *file = fopen(argv[1], "r");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file\n");
        exit(1);
    }

    load_c_config(file);

    fclose(file);

    // Sets the number of threads used in parallelization to one less than the maximum
    // This allows for the main thread to be used for other tasks
    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs - 1);

    return 0;
}
