#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "config.h"


// Global Variables are listed below, with their default values

// [CONTROL]
char c_calculation_type[50] = "superconductor"; 
char* get_calculation_type() {return c_calculation_type;}
char c_outdir[50] = "./output"; char* get_outdir() {return c_outdir;}
char c_prefix[50] = "sample"; char* get_prefix() {return c_prefix;}
char c_verbosity[50] = "high"; char* get_verbosity() {return c_verbosity;}

// [SYSTEM]
char c_interaction[50] = "scalapino"; char* get_interaction() {return c_interaction;}
int c_ibrav = 0;
float c_mu = 0.0;
float c_U = 4.0;

// [SUPERCONDUCTOR]
bool c_FS_only = true;
float c_wc = 0.05;
int c_num_eigenvalues_to_save = 1;
int c_frequency_pts = 5; // Number of BCS frequency points

// [MESH]
int c_k_mesh[3] = {20, 20, 20};
int c_q_mesh[3] = {20, 20, 20};
int c_w_pts = 11; // Number of frequency points

// [CELL]
float c_cell[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
float c_brillouin_zone[3][3] = {{2 * M_PI, 0.0, 0.0}, {0.0, 2 * M_PI, 0.0}, {0.0, 0.0, 2 * M_PI}};
               
// [BANDS]

// INPUTS DERIVED FROM THE CONFIG ABOVE
float c_max_freq = 10.0; 
int c_dim = 3; // Number of dimensions)
char c_band_name[50] = "sphere"; char* get_band_name() {return c_band_name;}


void get_dimensions() {
    // Check if bottom row and right column are zero
    if (c_cell[2][0] == 0.0 && c_cell[2][1] == 0.0 && c_cell[2][2] == 0.0 
            && c_cell[0][2] == 0.0 && c_cell[1][2] == 0.0 && c_cell[2][2] == 0.0) {
        c_dim = 2;
    } else {
        c_dim = 3;
    }
}

void cell_to_BZ(const float cell[3][3], float brillouin_zone[3][3]);
void set_string(char *dest, const char *src) {
    int size = 50;
    strncpy(dest, src, size - 1);
    dest[size - 1] = '\0';  // Ensure null-termination
}

void load_c_config(FILE *file) {
    char line[256];
    char key[50];
    char value[50];
    char section[50];
    int row = 0;
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "=") != NULL) {
            sscanf(line, "%[^=]=%s", key, value);
            if (strstr(line, "[CELL]") != NULL) {
                set_string(section, "CELL");
                continue;  // Skip to the next file line
            }

            // If we're in the [CELL] section, read matrix values
            if (strstr(section, "CELL") != NULL && strlen(line) > 1) {
                sscanf(line, "%f %f %f", &c_cell[row][0], &c_cell[row][1], &c_cell[row][2]);
                row++;
                // Stop reading after filling 3 rows
                if (row == 3) break;
            }

            if (strstr(key, "calculation") != NULL) {
                set_string(c_calculation_type, value);
            } else if (strstr(key, "outdir") != NULL) {
                set_string(c_outdir, value);
            } else if (strstr(key, "interaction") != NULL) {
                set_string(c_interaction, value);
            } else if (strstr(key, "FS_only") != NULL) {
                c_FS_only = atoi(value);
            } else if (strstr(key, "surface_pnts") != NULL) {
                sscanf(line, "%d %d %d", &c_k_mesh[0], &c_k_mesh[1], &c_k_mesh[2]);
            } else if (strstr(key, "potential_pnts") != NULL) {
                sscanf(line, "%d %d %d", &c_q_mesh[0], &c_q_mesh[1], &c_q_mesh[2]);
            } else if (strstr(key, "frequency_pts") != NULL) {
                c_w_pts = atoi(value);
            } else {
                printf("Invalid key: %s\n", key);
                exit(1);
            }
        }
    }
    cell_to_BZ(c_cell, c_brillouin_zone);
    get_dimensions();
}
