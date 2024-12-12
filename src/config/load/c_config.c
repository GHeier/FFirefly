#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "c_config.h"


// Global Variables are listed below, with their default values

//[CONTROL]
char* c_category = "test";
char* get_category() {return c_category;}
char* c_calculation = "test";
char* get_calculation() {return c_calculation;}
char* c_outdir = "./output";
char* get_outdir() {return c_outdir;}
char* c_prefix = "sample";
char* get_prefix() {return c_prefix;}
char* c_verbosity = "high";
char* get_verbosity() {return c_verbosity;}
char* c_datfile_in = "input.dat";
char* get_datfile_in() {return c_datfile_in;}
char* c_datfile_out = "output.dat";
char* get_datfile_out() {return c_datfile_out;}

//[SYSTEM]
char* c_interaction = "none";
char* get_interaction() {return c_interaction;}
int c_dimension = 3;
int c_ibrav = 0;
int c_nbnd = 1;
float c_fermi_energy = 0.0;
float c_Temperature = 0.0;
float c_onsite_U = 0.0;

//[MESH]
int c_k_mesh[3] = {10, 10, 10};
int c_q_mesh[3] = {10, 10, 10};
int c_w_pts = 100;

//[CELL]
float c_cell[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
float c_brillouin_zone[3][3] = {{6.283185307179586, 0.0, 0.0}, {0.0, 6.283185307179586, 0.0}, {0.0, 0.0, 6.283185307179586}};

//[BANDS]
char** c_band = (char*[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
char** get_band() {return (char**)c_band;}
float c_eff_mass[50];
float c_t0[50];
float c_t1[50];
float c_t2[50];
float c_t3[50];
float c_t4[50];
float c_t5[50];
float c_t6[50];
float c_t7[50];
float c_t8[50];
float c_t9[50];
float c_t10[50];

//[SUPERCONDUCTOR]
bool c_FS_only = true;
float c_bcs_cutoff_frequency = 0.05;
int c_num_eigenvalues_to_save = 1;
int c_frequency_pts = 5;

//[RESPONSE]
bool c_dynamic = false;
// End of Global Variables


void get_dimensions() {
    // Check if bottom row and right column are zero
    if (c_cell[2][0] == 0.0 && c_cell[2][1] == 0.0 && c_cell[2][2] == 0.0 
            && c_cell[0][2] == 0.0 && c_cell[1][2] == 0.0 && c_cell[2][2] == 0.0) {
        c_dimension = 2;
    } else {
        c_dimension = 3;
    }
}

#define PI 3.141592653589793

void cross_product(const float a[3], const float b[3], float result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

float dot_product(const float a[3], const float b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void cell_to_BZ(float ucell[3][3], float (*bz_matrix)[3]) {
    float volume;
    float a[3], b[3], c[3], cross_bc[3], cross_ca[3], cross_ab[3];
    
    // Extract lattice vectors
    for (int i = 0; i < 3; i++) {
        a[i] = ucell[i][0];
        b[i] = ucell[i][1];
        c[i] = ucell[i][2];
    }
    
    // Calculate the volume of the unit cell using a dot product and cross product
    cross_product(b, c, cross_bc);
    volume = dot_product(a, cross_bc);
    
    // Calculate the reciprocal lattice vectors
    cross_product(c, a, cross_ca);
    cross_product(a, b, cross_ab);
    
    // Calculate the BZ matrix
    for (int i = 0; i < 3; i++) {
        bz_matrix[i][0] = 2.0f * PI * cross_bc[i] / volume;
        bz_matrix[i][1] = 2.0f * PI * cross_ca[i] / volume;
        bz_matrix[i][2] = 2.0f * PI * cross_ab[i] / volume;
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

void load_default_band_values() {
    strcpy(c_band[0], "fermi_gas");
    c_eff_mass[0] = 1.0;
}


void load_c_config() {
    char line[256];
    char key[50];
    char value[50];
    char section[50];
    int row = 0;
    int n = 0;
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
            sscanf(line, " %s = %s", key, value);

            // Read in variable values from the config file

//[CONTROL]
            if (strstr(key, "category") != NULL) {
                set_string(c_category, value);
            }
            else if (strstr(key, "calculation") != NULL) {
                set_string(c_calculation, value);
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

//[SYSTEM]
            else if (strstr(key, "interaction") != NULL) {
                set_string(c_interaction, value);
            }
            else if (strstr(key, "dimension") != NULL) {
                c_dimension = atoi(value);
            }
            else if (strstr(key, "ibrav") != NULL) {
                c_ibrav = atoi(value);
            }
            else if (strstr(key, "nbnd") != NULL) {
                c_nbnd = atoi(value);
            }
            else if (strstr(key, "fermi_energy") != NULL) {
                c_fermi_energy = atof(value);
            }
            else if (strstr(key, "Temperature") != NULL) {
                c_Temperature = atof(value);
            }
            else if (strstr(key, "onsite_U") != NULL) {
                c_onsite_U = atof(value);
            }

//[MESH]
            else if (strstr(key, "k_mesh") != NULL) {
                sscanf(line, " k_mesh = %d %d %d", &c_k_mesh[0], &c_k_mesh[1], &c_k_mesh[2]);
            }
            else if (strstr(key, "q_mesh") != NULL) {
                sscanf(line, " q_mesh = %d %d %d", &c_q_mesh[0], &c_q_mesh[1], &c_q_mesh[2]);
            }
            else if (strstr(key, "w_pts") != NULL) {
                c_w_pts = atoi(value);
            }

//[CELL]



//[BANDS]
            else if (strstr(key, "band") != NULL) {
                n = atoi(key + 4);
                set_string(c_band[n], value);
            }
            else if (strstr(key, "eff_mass") != NULL) {
                c_eff_mass[n] = atof(value);
            }
            else if (strstr(key, "t0") != NULL) {
                c_t0[n] = atof(value);
            }
            else if (strstr(key, "t1") != NULL) {
                c_t1[n] = atof(value);
            }
            else if (strstr(key, "t2") != NULL) {
                c_t2[n] = atof(value);
            }
            else if (strstr(key, "t3") != NULL) {
                c_t3[n] = atof(value);
            }
            else if (strstr(key, "t4") != NULL) {
                c_t4[n] = atof(value);
            }
            else if (strstr(key, "t5") != NULL) {
                c_t5[n] = atof(value);
            }
            else if (strstr(key, "t6") != NULL) {
                c_t6[n] = atof(value);
            }
            else if (strstr(key, "t7") != NULL) {
                c_t7[n] = atof(value);
            }
            else if (strstr(key, "t8") != NULL) {
                c_t8[n] = atof(value);
            }
            else if (strstr(key, "t9") != NULL) {
                c_t9[n] = atof(value);
            }
            else if (strstr(key, "t10") != NULL) {
                c_t10[n] = atof(value);
            }

//[SUPERCONDUCTOR]
            else if (strstr(key, "FS_only") != NULL) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_FS_only = true;
                } else {
                    c_FS_only = false;
                }
            }
            else if (strstr(key, "bcs_cutoff_frequency") != NULL) {
                c_bcs_cutoff_frequency = atof(value);
            }
            else if (strstr(key, "num_eigenvalues_to_save") != NULL) {
                c_num_eigenvalues_to_save = atoi(value);
            }
            else if (strstr(key, "frequency_pts") != NULL) {
                c_frequency_pts = atoi(value);
            }

//[RESPONSE]
            else if (strstr(key, "dynamic") != NULL) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_dynamic = true;
                } else {
                    c_dynamic = false;
                }
            }
            // End of variable reading
            else {
                printf("%s = %s\n", key, value);
                printf("Invalid key: %s\n", key);
                exit(1);
            }
        }
    }
    c_nbnd = n;
    cell_to_BZ(c_cell, c_brillouin_zone);
    get_dimensions();
    printf("Loaded Config\n");
}
