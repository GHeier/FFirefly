#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "c_config.h"

// Global Variables are listed below, with their default values

//[CONTROL]
char* c_category = "test";
char* get_category() {return c_category;}
char* c_calculation = "test";
char* get_calculation() {return c_calculation;}
char* c_method = "none";
char* get_method() {return c_method;}
char* c_outdir = "./";
char* get_outdir() {return c_outdir;}
char* c_indir = "./";
char* get_indir() {return c_indir;}
char* c_prefix = "sample";
char* get_prefix() {return c_prefix;}
char* c_verbosity = "low";
char* get_verbosity() {return c_verbosity;}
bool c_automatic_file_read = true;
bool c_write_result = true;
char* c_filetype = "h5";
char* get_filetype() {return c_filetype;}

//[SYSTEM]
char* c_interaction = "none";
char* get_interaction() {return c_interaction;}
int c_dimension = 3;
char* c_celltype = "";
char* get_celltype() {return c_celltype;}
int c_nbnd = 0;
int c_natoms = 0;
float c_fermi_energy = 0.0;
float c_Temperature = 0.0;
float c_onsite_U = 0.0;
float c_cutoff_energy = 0.05;

//[MESH]
int c_k_mesh[3] = {10, 10, 10};
int c_q_mesh[3] = {10, 10, 10};
int c_w_pts = 100;

//[CELL]
float c_cell[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

//[BRILLOUIN_ZONE]
float c_brillouin_zone[3][3] = {{6.283185307179586, 0.0, 0.0}, {0.0, 6.283185307179586, 0.0}, {0.0, 0.0, 6.283185307179586}};

//[ATOMS]
char** c_atom = (char*[]){0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
char** get_atom() {return (char**)c_atom;}
float c_position[50][3] = {0};

//[HAMILTONIAN]
char* c_type = "tight_binding";
char* get_type() {return c_type;}
int c_basis_size = 1;
float c_H0[1][1] = {{1.0}};

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
int c_num_eigenvalues_to_save = 1;
int c_frequency_pts = 5;
char* c_projections = "";
char* get_projections() {return c_projections;}

//[RESPONSE]
bool c_dynamic = false;

//[MANY_BODY]
bool c_self_consistent = false;
// End of Global Variables

void get_dimensions() {

    // Check if bottom row and right column are zero
    if (c_cell[2][0] == 0.0 && c_cell[2][1] == 0.0 && c_cell[2][2] == 0.0 &&
        c_cell[0][2] == 0.0 && c_cell[1][2] == 0.0 && c_cell[2][2] == 0.0) {
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

    // Calculate the volume of the unit cell using a dot product and cross
    // product
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

// void cell_to_BZ(const float cell[3][3], float brillouin_zone[3][3]);
void strip_single_quotes(char *str) {
    int i, j = 0;
    for (i = 0; str[i] != '\0'; i++) {
        if (str[i] != '\'') {
            str[j++] = str[i]; // Copy non-single-quote characters
        }
    }
    str[j] = '\0'; // Null-terminate the modified string
}
void set_string(char **dest, const char *src) {
    int size = 80;
    *dest = (char *)malloc(size * sizeof(char)); // Allocate memory for dest
    strncpy(*dest, src, size - 1);               // Copy up to size-1 characters
    (*dest)[size - 1] = '\0';                    // Ensure null-termination
    strip_single_quotes(*dest); // Strip single quotes if present
}

void set_section(char *dest, const char *src) {
    int size = 50;
    strncpy(dest, src, size - 1); // Copy up to size-1 characters
    dest[size - 1] = '\0';        // Ensure null-termination
    strip_single_quotes(dest);    // Strip single quotes if present
}

void load_default_band_values() {
    c_band[0] =
        (char *)malloc(50 * sizeof(char)); // Allocating space for 50 characters
    strcpy(c_band[0], "fermi_gas");
    c_eff_mass[0] = 1.0;
}

void make_save_file() {
    char path[PATH_MAX]; // Buffer to hold the executable path

    // Read the symbolic link for the executable
    ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
    path[len - 5] = '\0';      // Remove the executable name
    strcat(path, "input.cfg"); // Append the input file name
    FILE *file = fopen(path, "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }
    char line[256];
    // Check if stdin is a terminal
    if (!isatty(fileno(stdin))) {
        while (fgets(line, sizeof(line), stdin) != NULL) {
            fputs(line, file);
        }
    }
    fclose(file);
}

void read_c_config(const char *path) {
    char line[256];
    char key[50];
    char value[50];
    char section[50];
    int row = 0;
    int n = 0;
    int atom_ind = 0;
    bool got_dimension = false;
    bool got_bz = false;
    bool got_nbnd = false;
    FILE *file = fopen(path, "r");
    if (file == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }
    while (fgets(line, sizeof(line), file) != NULL) {
        if (strstr(line, "[CELL]") != NULL) {
            set_section(section, "CELL");
            continue; // Skip to the next file line
        }
        if (strstr(line, "[BRILLOUIN_ZONE]") != NULL) {
            set_section(section, "BRILLOUIN_ZONE");
            got_bz = true;
            continue; // Skip to the next file line
        }
        // If we're in the [CELL] section, read matrix values
        if (strstr(section, "CELL") != NULL && strlen(line) > 1) {
            sscanf(line, "%f %f %f", &c_cell[row][0], &c_cell[row][1],
                   &c_cell[row][2]);
            row++;
            // Stop reading after filling 3 rows
            if (row == 3) {
                section[0] = '\0';
                continue;
            }
        }
        if (strstr(section, "BRILLOUIN_ZONE") != NULL && strlen(line) > 1) {
            sscanf(line, "%f %f %f", &c_brillouin_zone[row][0],
                   &c_brillouin_zone[row][1], &c_brillouin_zone[row][2]);
            atom_ind++;
            // Stop reading after filling 3 rows
            if (line[0] == '\0') {
                section[0] = '\0';
                continue;
            }
        }
        if (strstr(section, "ATOMS") != NULL && strlen(line) > 1) {
            char temp;
            sscanf(line, "%c %f %f %f", &temp, &c_position[atom_ind][0],
                   &c_position[atom_ind][1], &c_position[atom_ind][2]);
            set_string(&c_atom[atom_ind], &temp);
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
                set_string(&c_category, value);
            }
            else if (strstr(key, "calculation") != NULL) {
                set_string(&c_calculation, value);
            }
            else if (strstr(key, "method") != NULL) {
                set_string(&c_method, value);
            }
            else if (strstr(key, "outdir") != NULL) {
                set_string(&c_outdir, value);
            }
            else if (strstr(key, "indir") != NULL) {
                set_string(&c_indir, value);
            }
            else if (strstr(key, "prefix") != NULL) {
                set_string(&c_prefix, value);
            }
            else if (strstr(key, "verbosity") != NULL) {
                set_string(&c_verbosity, value);
            }
            else if (strstr(key, "automatic_file_read") != NULL) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_automatic_file_read = true;
                } else {
                    c_automatic_file_read = false;
                }
            }
            else if (strstr(key, "write_result") != NULL) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_write_result = true;
                } else {
                    c_write_result = false;
                }
            }
            else if (strstr(key, "filetype") != NULL) {
                set_string(&c_filetype, value);
            }

//[SYSTEM]
            else if (strstr(key, "interaction") != NULL) {
                set_string(&c_interaction, value);
            }
            else if (strstr(key, "dimension") != NULL) {
                c_dimension = atoi(value);
                 got_dimension = true;
            }
            else if (strstr(key, "celltype") != NULL) {
                set_string(&c_celltype, value);
            }
            else if (strstr(key, "nbnd") != NULL) {
                c_nbnd = atoi(value);
                 got_nbnd = true;
            }
            else if (strstr(key, "natoms") != NULL) {
                c_natoms = atoi(value);
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
            else if (strstr(key, "cutoff_energy") != NULL) {
                c_cutoff_energy = atof(value);
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

//[BRILLOUIN_ZONE]

//[ATOMS]

//[BANDS]
            else if (strstr(key, "band") != NULL) {
                n = atoi(key + 4)-1;
                set_string(&c_band[n], value);
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
            else if (strstr(key, "num_eigenvalues_to_save") != NULL) {
                c_num_eigenvalues_to_save = atoi(value);
            }
            else if (strstr(key, "frequency_pts") != NULL) {
                c_frequency_pts = atoi(value);
            }
            else if (strstr(key, "projections") != NULL) {
                set_string(&c_projections, value);
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

//[MANY_BODY]
            else if (strstr(key, "self_consistent") != NULL) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_self_consistent = true;
                } else {
                    c_self_consistent = false;
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
    if (!got_nbnd && c_band[0] != 0)
        c_nbnd = n + 1;
    if (!got_bz)
        cell_to_BZ(c_cell, c_brillouin_zone);
    if (c_dimension == 2) {
        // c_cell[2][0] = 0.0;
        // c_cell[2][1] = 0.0;
        // c_cell[2][2] = 0.0;
        // c_cell[0][2] = 0.0;
        // c_cell[1][2] = 0.0;
        // c_cell[2][2] = 0.0;
        // c_brillouin_zone[2][0] = 0.0;
        // c_brillouin_zone[2][1] = 0.0;
        // c_brillouin_zone[2][2] = 0.0;
        // c_brillouin_zone[0][2] = 0.0;
        // c_brillouin_zone[1][2] = 0.0;
        // c_brillouin_zone[2][2] = 0.0;
        c_k_mesh[2] = 1;
        c_q_mesh[2] = 1;
    }
    if (!got_dimension)
        get_dimensions();
    int outdir_length = strlen(c_outdir);
    if (c_outdir[outdir_length - 1] != '/') {
        char *temp_dir = malloc(outdir_length + 2);
        strcpy(temp_dir, c_outdir);
        temp_dir[outdir_length] = '/';
        temp_dir[outdir_length + 1] = '\0';
        c_outdir = temp_dir;
    }
    int result = mkdir_p(c_outdir, 0755);
    if (result != 0)
        printf("Could not make outdir directory\n");
    result = mkdir_p(c_indir, 0755);
    if (result != 0)
        printf("Could not make indir directory\n");
}

void load_c_config() {
    char path[PATH_MAX]; // Buffer to hold the executable path
    ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
    path[len - 5] = '\0';      // Remove the executable name
    strcat(path, "input.cfg"); // Append the input file name
    make_save_file();
    read_c_config(path);
}

void unload_c_config() {
    // Unload the config file
    //
    ////[CONTROL]
    //    free(c_category);
    //    free(c_calculation);
    //    free(c_outdir);
    //    free(c_indir);
    //    free(c_prefix);
    //    free(c_verbosity);
    //    free(c_input_data_file);
    //    free(c_output_data_file);
    //    for (int i = 0; i < 50; i++) {
    //        free(c_automatic_file_read[i]);
    //    }
    //
    ////[SYSTEM]
    //    free(c_interaction);
    //    for (int i = 0; i < 50; i++) {
    //        free(c_dimension[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_ibrav[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_nbnd[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_fermi_energy[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_Temperature[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_onsite_U[i]);
    //    }
    //
    ////[MESH]
    //    for (int i = 0; i < 50; i++) {
    //        free(c_k_mesh[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_q_mesh[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_w_pts[i]);
    //    }
    //
    ////[CELL]
    //    for (int i = 0; i < 50; i++) {
    //        free(c_cell[i]);
    //    }
    //
    ////[BRILLOUIN_ZONE]
    //    for (int i = 0; i < 50; i++) {
    //        free(c_brillouin_zone[i]);
    //    }
    //
    ////[ATOMS]
    //    free(c_atom);
    //    for (int i = 0; i < 50; i++) {
    //        free(c_position[i]);
    //    }
    //
    ////[BANDS]
    //    free(c_band);
    //    for (int i = 0; i < 50; i++) {
    //        free(c_eff_mass[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t0[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t1[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t2[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t3[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t4[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t5[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t6[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t7[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t8[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t9[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_t10[i]);
    //    }
    //
    ////[SUPERCONDUCTOR]
    //    free(c_method);
    //    for (int i = 0; i < 50; i++) {
    //        free(c_FS_only[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_cutoff_energy[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_num_eigenvalues_to_save[i]);
    //    }
    //    for (int i = 0; i < 50; i++) {
    //        free(c_frequency_pts[i]);
    //    }
    //    free(c_projections);
    //
    ////[RESPONSE]
    //    for (int i = 0; i < 50; i++) {
    //        free(c_dynamic[i]);
    //    }
    // End of unloading the config file
}

#include <stdarg.h>
#include <stdio.h>

// Print function with color and printf-style formatting
void printcolor(Color color, const char *format, ...) {
    // Print color escape code
    printf("\033[1;%dm", color);

    // Handle variadic arguments
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);

    // Reset color
    printf("\033[0m");
}

bool print_test_results(bool all_tests[], int num_tests,
                        const char *test_name) {
    int tests_passed = 0;
    for (int i = 0; i < num_tests; i++)
        if (all_tests[i])
            tests_passed++;
    if (tests_passed == num_tests) {
        printcolor(GREEN, "All %s passed!\n", test_name);
        return true;
    } else {
        printcolor(RED, " - %d/%d %s passed\n", tests_passed, num_tests,
                   test_name);
        return false;
    }
}
int mkdir_p(const char *path, mode_t mode) {
    char *tmp = strdup(path);
    char *p = NULL;
    size_t len;
    int status = 0;

    if (tmp == NULL)
        return -1;

    len = strlen(tmp);
    if (len == 0) {
        free(tmp);
        return -1;
    }

    if (tmp[len - 1] == '/')
        tmp[len - 1] = '\0';

    for (p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(tmp, mode) != 0) {
                if (errno != EEXIST) {
                    status = -1;
                    break;
                }
            }
            *p = '/';
        }
    }

    if (status == 0) {
        if (mkdir(tmp, mode) != 0) {
            if (errno != EEXIST)
                status = -1;
        }
    }

    free(tmp);
    return status;
}
