#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "config/load/c_config.h"
#include "config/load/py_interface.h"

// Category nodes below
#include "algorithms/electron_number.hpp"
#include "hamiltonian/fs.hpp"
#include "many_body/vertex.hpp"
#include "many_body/self_energy.hpp"
#include "many_body/renormalization.hpp"
#include "many_body/node.hpp"
#include "response/DOS.h"
#include "response/response.h"
#include "superconductor/node.hpp"

// Test nodes below
#include "hamiltonian/tests/all.hpp"
#include "objects/tests/all.hpp"
#include "algorithms/tests/all.hpp"

#define MAX_TOKENS 10 // Max number of substrings
#define MAX_LENGTH 50 // Max length of each substring

// Global category calls
void DOS() {
    printf("Starting DOS Calculation\n\n");
    DOS_spectrum();
}

void electron_number() {
    printf("Starting Electron Number Calculation\n\n");
    electron_number_wrapper();
}

void fermi_surface() {
    printf("Starting Fermi Surface Calculation\n\n");
    save_FS();
}

void response() {
    printf("Starting Response Calculation\n\n");
    response_wrapper();
}

void vertex() {
    printf("Starting Vertex Calculation\n\n");
    vertex_wrapper();
}

void self_energy() {
    printf("Starting Self Energy Calculation\n\n");
    self_energy_wrapper();
}

void renormalization() {
    printf("Starting Self Energy Calculation\n\n");
    renormalization_wrapper();
}

void many_body() {
    printf("Starting Bethe-Salpeter Calculation\n\n");
    many_body_wrapper();
}

void superconductor() {
    printf("Starting Superconductor Calculation\n\n");
    superconductor_wrapper();
}

void print_banner_top() {
    printf("╔══════════════════════════════════════╗\n");
    printf("║       🚀 FFirefly Launching...       ║\n");
}

void print_banner_bottom(time_t now) {
    struct tm *t = localtime(&now);
    char time_str[16];
    strftime(time_str, sizeof(time_str), "%H:%M:%S", t);
    printf("║        Launched at %s          ║\n", time_str);
    printf("╚══════════════════════════════════════╝\n");
}

// Global test calls
void test() {
    printf("Starting Test Calculations\n");

    int num_tests = 3;

    bool all_tests[num_tests];
    all_tests[0] = hamiltonian_tests();
    all_tests[1] = object_tests();
    all_tests[2] = algorithm_tests();

    print_test_results(all_tests, num_tests, "tests");
}

int main() {
    time_t start = time(NULL);
    print_banner_top();

    load_c_config();           // Read input to load global c variables
    start_python();

    // Handling '+' separated categories for sequential runs
    char str_copy[MAX_LENGTH];
    char meth_copy[MAX_LENGTH];
    strncpy(str_copy, c_category, MAX_LENGTH); // Copy string safely
    strncpy(meth_copy, c_method, MAX_LENGTH); // Copy string safely
    char *tokens[MAX_TOKENS];
    char *meth_tokens[MAX_TOKENS];
    int count = 0;
    int mcount = 0;

    char *saveptr1, *saveptr2;
    char *token = strtok_r(str_copy, "+", &saveptr1);
    char *meth_token = strtok_r(meth_copy, "+", &saveptr2);

    while (token != NULL && count < MAX_TOKENS) {
        tokens[count] = strdup(token);
        count++;
        token = strtok_r(NULL, "+", &saveptr1);
    }
    while (meth_token != NULL && mcount < MAX_TOKENS) {
        meth_tokens[mcount] = strdup(meth_token);
        mcount++;
        meth_token = strtok_r(NULL, "+", &saveptr2);
    }

    print_banner_bottom(start);

    // Set the number of threads used in parallelization to one less than the
    // maximum
    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs - 1);
    if (!strcmp(c_verbosity, "high"))
        printf("Number of threads used in CPU parallelization: %d\n",
               num_procs - 1);

    for (int i = 0; i < count; i++) {
        char *category = tokens[i];
        int ind = i;
        if (ind >= mcount)
            ind = mcount - 1;
        if (c_method[0] != '\0')
            c_method = meth_tokens[ind];
        load_cpp_config_wrapper(); 
    /*
        * ADDING A CATEGORY OCCURS BELOW
        * FOLLOW THE PATTERN
    */
        if (!strcmp(category, "DOS"))
            DOS();
        else if (!strcmp(category, "fermi_surface"))
            fermi_surface();
        else if (!strcmp(category, "electron_num"))
            electron_number();
        else if (!strcmp(category, "response"))
            response();
        else if (!strcmp(category, "superconductor"))
            superconductor();
        else if (!strcmp(category, "vertex"))
            vertex();
        else if (!strcmp(category, "self_energy"))
            self_energy();
        else if (!strcmp(category, "renormalization"))
            renormalization();
        else if (!strcmp(category, "many_body"))
            many_body();
        else if (!strcmp(category, "test"))
            test();
        else
            printf("Unknown Category\n");
    }

    // unload_c_config(); // Free memory allocated for global c variables
    for (int i = 0; i < count; i++) {
        free(tokens[i]);
    }
    end_python();

    time_t end_time = time(NULL);
    int runtime_secs = (int)difftime(end_time, start);
    int hours = runtime_secs / 3600;
    int minutes = (runtime_secs % 3600) / 60;
    int seconds = runtime_secs % 60;

    printf("╔══════════════════════════════════════╗\n");
    printf("║         🏁 FFirefly Landed.          ║\n");
    printf("║          Runtime: %02d:%02d:%02d           ║\n", hours, minutes, seconds);
    printf("╚══════════════════════════════════════╝\n");
    return 0;
}
