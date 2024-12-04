#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "c_config.h"


// Global Variables are listed below, with their default values

// [CONTROL] 
char c_category[50] = "test"; char* get_category() {return c_category;}
char c_outdir[50] = "./output"; char* get_outdir() {return c_outdir;}
char c_prefix[50] = "sample"; char* get_prefix() {return c_prefix;}
char c_verbosity[50] = "high"; char* get_verbosity() {return c_verbosity;}
char c_datfile_in[50] = "input.dat"; char* get_datfile_in() {return c_datfile_in;}
char c_datfile_out[50] = "output.dat"; char* get_datfile_out() {return c_datfile_out;}

// [SYSTEM]
char c_interaction[50] = "scalapino"; char* get_interaction() {return c_interaction;}
int c_ibrav = 0;
float c_fermi_energy = 0.0;
float c_onsite_U = 4.0;

// [SUPERCONDUCTOR]
bool c_FS_only = true;
float c_bcs_cutoff_frequency = 0.05;
int c_num_eigenvalues_to_save = 1;
int c_frequency_pts = 5; // Number of BCS frequency points

// [RESPONSE]
char c_calculation[50] = "polarization"; char* get_calculation() {return c_calculation;}
bool c_dynamic = false;

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


void foo() {
    printf("foo\n");
}

void get_dimensions() {
    // Check if bottom row and right column are zero
    if (c_cell[2][0] == 0.0 && c_cell[2][1] == 0.0 && c_cell[2][2] == 0.0 
            && c_cell[0][2] == 0.0 && c_cell[1][2] == 0.0 && c_cell[2][2] == 0.0) {
        c_dim = 2;
    } else {
        c_dim = 3;
    }
}

//void cell_to_BZ(const float cell[3][3], float brillouin_zone[3][3]);
void strip_single_quotes(char *str) {
    int i, j = 0;
    for (i = 0; str[i] != '\0'; i++) {
        if (str[i] != '\'') {
            str[j++] = str[i]; // Copy non-single-quote characters
        }
    }
    str[j] = '\0'; // Null-terminate the modified string
}
void set_string(char *dest, const char *src) {
    int size = 50;
    strncpy(dest, src, size - 1);
    strip_single_quotes(dest);
    dest[size - 1] = '\0';  // Ensure null-termination
}

void load_c_config() {
    char line[256];
    char key[50];
    char value[50];
    char section[50];
    int row = 0;
    while (fgets(line, sizeof(line), stdin) != NULL) {
        if (strstr(line, "[CELL]") != NULL) {
            set_string(section, "CELL");
            continue;  // Skip to the next file line
        }
        // If we're in the [CELL] section, read matrix values
        if (strstr(section, "CELL") != NULL && strlen(line) > 1) {
            sscanf(line, "%f %f %f", &c_cell[row][0], &c_cell[row][1], &c_cell[row][2]);
            row++;
            // Stop reading after filling 3 rows
            if (row == 3) {
                section[0] = '\0';
                continue;
            }
        }
        if (strstr(line, "=") != NULL) {
            sscanf(line, "%[^=]=%s", key, value);

            // [CONTROL]
            if (strstr(key, "category") != NULL) {
                set_string(c_category, value);
            } 
            else if (strstr(key, "outdir") != NULL) {
                set_string(c_outdir, value);
            } 
            else if (strstr(key, "prefix") != NULL) {
                set_string(c_prefix, value);
            } 
            else if (strstr(key, "verbosity") != NULL) {
                set_string(c_verbosity, value);
            } 
            else if (strstr(key, "datfile_in") != NULL) {
                set_string(c_datfile_in, value);
            } 
            else if (strstr(key, "datfile_out") != NULL) {
                set_string(c_datfile_out, value);
            }
            // [SYSTEM]
            else if (strstr(key, "interaction") != NULL) {
                set_string(c_interaction, value);
            } 
            else if (strstr(key, "ibrav") != NULL) {
                c_ibrav = atoi(value);
            } 
            else if (strstr(key, "fermi_energy") != NULL) {
                c_fermi_energy = atof(value);
            } 
            else if (strstr(key, "onsite_U") != NULL) {
                c_onsite_U = atof(value);
            } 
            // [MESH]
            else if (strstr(key, "k_mesh") != NULL) {
                sscanf(line, " k_mesh = %d %d %d", &c_k_mesh[0], &c_k_mesh[1], &c_k_mesh[2]);
            } 
            else if (strstr(key, "q_mesh") != NULL) {
                sscanf(line, " q_mesh = %d %d %d", &c_q_mesh[0], &c_q_mesh[1], &c_q_mesh[2]);
            } 
            else if (strstr(key, "w_pts") != NULL) {
                c_w_pts = atoi(value);
            } 
            // [SUPERCONDUCTOR]
            else if (strstr(key, "bcs_cutoff_frequency") != NULL) {
                c_bcs_cutoff_frequency = atof(value);
            } 
            else if (strstr(key, "num_eigenvalues_to_save") != NULL) {
                c_num_eigenvalues_to_save = atoi(value);
            } 
            else if (strstr(key, "FS_only") != NULL) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_FS_only = true;
                } else {
                    c_FS_only = false;
                }
            } 
            // [RESPONSE]
            else if (strstr(key, "calculation") != NULL) {
                set_string(c_calculation, value);
            } 
            else if (strstr(key, "dynamic") != NULL) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_dynamic = true;
                } else {
                    c_dynamic = false;
                }
            } 
            else {
                printf("Invalid key: %s\n", key);
                exit(1);
            }
        }
    }
    //cell_to_BZ(c_cell, c_brillouin_zone);
    get_dimensions();
    printf("Loaded Config\n");
}
