#ifndef _GNU_SOURCE
#define _GNU_SOURCE  // dladdr
#endif
#include <dlfcn.h>
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
bool c_debug = false;
char* c_prefix = "sample";
char* get_prefix() {return c_prefix;}
char* c_verbosity = "low";
char* get_verbosity() {return c_verbosity;}
bool c_automatic_file_read = true;
bool c_write_result = true;
char* c_filetype = "h5";
char* get_filetype() {return c_filetype;}
int c_mpi_procs = 0;

//[SYSTEM]
char* c_interaction = "none";
char* get_interaction() {return c_interaction;}
int c_dimension = 3;
char* c_celltype = "";
char* get_celltype() {return c_celltype;}
int c_nbnd = 0;
float c_fermi_energy = 0.0;
float c_num_electrons = 0.0;
bool c_mu_from_n = false;
float c_Temperature = 0.0;
float c_cutoff_energy = 0.05;
float c_smearing = 0.02;
float c_mixing = 0.02;
int c_max_iters = 100;
float c_qp_weight = 1.0;
int c_recurse_level = 0;
float c_tol = 0.0001;

//[HAMILTONIAN]
char* c_hamiltonian = "tight_binding";
char* get_hamiltonian() {return c_hamiltonian;}
float c_eps_dx2y2 = 0.0;
float c_eps_dz2 = 0.0;
float c_eps_px = 0.0;
float c_eps_py = 0.0;
float c_eps_pz = 0.0;
float c_delta_dp = 0.0;

//[HUBBARD]
float c_U0 = 0.0;
float c_U1 = 0.0;
float c_J0 = 0.0;
float c_J1 = 0.0;

//[MESH]
int c_k_mesh[3] = {10, 10, 10};
int c_q_mesh[3] = {10, 10, 10};
int c_w_pts = 100;

//[CELL]
float c_cell[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

//[BRILLOUIN_ZONE]
float c_brillouin_zone[3][3] = {{6.283185307179586, 0.0, 0.0}, {0.0, 6.283185307179586, 0.0}, {0.0, 0.0, 6.283185307179586}};

//[BASIS]
char** c_states = (char*[]){"H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"};
char** get_states() {return (char**)c_states;}
float c_positions[50][3] = {{0.0, 0.0, 0.0}};

//[BANDS]
char* c_band = "fermi_gas";
char* get_band() {return c_band;}
float c_eff_mass = 1.0;
float c_t0 = 1.0;
float c_t1 = 0.0;
float c_t2 = 0.0;
float c_t3 = 0.0;
float c_t4 = 0.0;
float c_t5 = 0.0;
float c_t6 = 0.0;
float c_t7 = 0.0;
float c_t8 = 0.0;
float c_t9 = 0.0;
float c_t10 = 0.0;
float c_tz0 = 0.0;
float c_tz1 = 0.0;
float c_tz2 = 0.0;
float c_tz3 = 0.0;
float c_tz4 = 0.0;

//[SUPERCONDUCTOR]
bool c_FS_only = true;
int c_num_eigenvalues_to_save = 5;
int c_frequency_pts = 0;
char* c_projections = "";
char* get_projections() {return c_projections;}

//[RESPONSE]
bool c_dynamic = false;

//[MANY_BODY]
bool c_self_consistent = false;
char* c_impurity_solver = "IPT";
char* get_impurity_solver() {return c_impurity_solver;}
int c_n_cycles = 20000;
int c_len_cycles = 200;
int c_warmup_cycles = 8000;
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

void make_lowercase(char *str) {
    if (str == NULL) {
        return;
    }
    for (int i = 0; str[i] != '\0'; i++) {
        if (str[i] >= 'A' && str[i] <= 'Z') {
            str[i] = str[i] + 32; // Convert uppercase to lowercase
        }
    }
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
    // c_band is already initialized to "fermi_gas" in the defaults
    // c_eff_mass is already initialized to 1.0 in the defaults
    // This function is kept for compatibility but does nothing
}

// Directory of the running executable, with trailing '/'
static void exe_dir(char *path, size_t size) {
    ssize_t len = readlink("/proc/self/exe", path, size - 1);
    if (len <= 0) {
        strcpy(path, "./");
        return;
    }
    path[len] = '\0';
    char *slash = strrchr(path, '/');
    if (slash) slash[1] = '\0';
}

static void input_cfg_path(char *path, size_t size) {
    exe_dir(path, size);
    strncat(path, "input.cfg", size - strlen(path) - 1);
}

// Path of the config last read by read_c_config ("" = none yet)
static char loaded_config_path[PATH_MAX] = "";

bool c_config_was_read(void) {
    return loaded_config_path[0] != '\0';
}

static bool readable(const char *path) {
    return access(path, R_OK) == 0;
}

// Directory of the binary or library this code was linked into, with trailing '/'
static bool library_dir(char *path, size_t size) {
    Dl_info info;
    if (!dladdr((void *)&firefly_config_path, &info) || info.dli_fname == NULL) return false;
    char resolved[PATH_MAX];
    if (realpath(info.dli_fname, resolved) == NULL) return false;
    char *slash = strrchr(resolved, '/');
    if (slash == NULL) return false;
    slash[1] = '\0';
    snprintf(path, size, "%s", resolved);
    return true;
}

const char *firefly_config_path(void) {
    static char path[PATH_MAX];
    if (c_config_was_read()) return loaded_config_path;
    const char *env = getenv("FIREFLY_CONFIG");
    if (env && env[0] != '\0') return env;
    input_cfg_path(path, sizeof(path));
    if (readable(path)) return path;
    if (library_dir(path, sizeof(path))) {
        strncat(path, "../bin/input.cfg", sizeof(path) - strlen(path) - 1);
        if (readable(path)) return path;
    }
    return NULL;
}

// Name of the config as the user gave it, for error messages (NULL = path read)
static const char *config_display_name = NULL;

// Copy the run configuration into input.cfg next to the executable, which is
// where method executables and the Python/Julia loaders read it from.
//   source == NULL : read from stdin (if it is not a terminal)
//   source == ""   : write an empty config (all defaults)
//   otherwise      : copy the named file
static void write_input_cfg(const char *source) {
    char path[PATH_MAX];
    input_cfg_path(path, sizeof(path));
    FILE *src = NULL;
    if (source && source[0] != '\0') {
        src = fopen(source, "r");
        if (src == NULL) {
            fprintf(stderr, "Error: cannot open config file '%s': %s\n", source, strerror(errno));
            exit(1);
        }
    } else if (source == NULL && !isatty(fileno(stdin))) {
        src = stdin;
    }
    config_display_name = (source && source[0] != '\0') ? source : (src == stdin ? "<stdin>" : NULL);
    FILE *file = fopen(path, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: cannot write '%s': %s\n", path, strerror(errno));
        exit(1);
    }
    char line[256];
    if (src) {
        while (fgets(line, sizeof(line), src) != NULL) {
            fputs(line, file);
        }
        if (src != stdin) fclose(src);
    }
    fclose(file);
}

void make_save_file() {
    write_input_cfg(NULL);
}

// True for blank lines and comment lines ('#', ';' or '!' first)
static bool is_comment_or_blank(const char *line) {
    while (*line == ' ' || *line == '\t') line++;
    return *line == '#' || *line == ';' || *line == '!' || *line == '\n' || *line == '\r' || *line == '\0';
}

// Split "key = value ..." into key and the first token of value. Spaces
// around '=' are optional. Returns false if the line has no key.
static bool split_key_value(const char *line, char *key, size_t key_size,
                            char *value, size_t value_size) {
    const char *eq = strchr(line, '=');
    if (!eq) return false;
    const char *k0 = line;
    while (k0 < eq && (*k0 == ' ' || *k0 == '\t')) k0++;
    const char *k1 = eq;
    while (k1 > k0 && (k1[-1] == ' ' || k1[-1] == '\t')) k1--;
    size_t klen = k1 - k0;
    if (klen == 0 || klen >= key_size) return false;
    memcpy(key, k0, klen);
    key[klen] = '\0';
    value[0] = '\0';
    // First whitespace-delimited token, bounded by value_size
    const char *v = eq + 1;
    while (*v == ' ' || *v == '\t') v++;
    size_t vlen = strcspn(v, " \t\r\n");
    if (vlen >= value_size) vlen = value_size - 1;
    memcpy(value, v, vlen);
    value[vlen] = '\0';
    return true;
}

void read_c_config(const char *path) {
    char line[256];
    char key[256];
    char value[256];
    char section[50];
    int row = 0;
    int n = 0;
    int state_ind = 0;
    bool got_dimension = false;
    bool got_bz = false;
    bool got_nbnd = false;
    int nbnd = 0;
    FILE *file = fopen(path, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: cannot open config file '%s': %s\n", path, strerror(errno));
        exit(1);
    }
    if (realpath(path, loaded_config_path) == NULL)
        snprintf(loaded_config_path, sizeof(loaded_config_path), "%s", path);
    int line_no = 0;
    while (fgets(line, sizeof(line), file) != NULL) {
        line_no++;
        if (is_comment_or_blank(line)) continue;
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
            state_ind++;
            // Stop reading after filling 3 rows
            if (line[0] == '\0') {
                section[0] = '\0';
                continue;
            }
        }
        if (strstr(section, "BASIS") != NULL && strlen(line) > 1) {
            char state_name[256];
            if (sscanf(line, "%s %f %f %f", state_name, &c_positions[nbnd][0],
                   &c_positions[nbnd][1], &c_positions[nbnd][2]) == 4) {
                set_string(&c_states[nbnd], state_name);
                nbnd++;
                // Stop reading after filling 50 states
                if (nbnd >= 50) {
                    section[0] = '\0';
                    continue;
                }
            }
        }
        if (strstr(line, "=") != NULL) {
            if (!split_key_value(line, key, sizeof(key), value, sizeof(value))) {
                fprintf(stderr, "Error: %s:%d: cannot parse line: %s",
                        config_display_name ? config_display_name : path, line_no, line);
                exit(1);
            }

            // Read in variable values from the config file

//[CONTROL]
            if (strcmp(key, "category") == 0) {
                set_string(&c_category, value);
            }
            else if (strcmp(key, "calculation") == 0) {
                set_string(&c_calculation, value);
            }
            else if (strcmp(key, "method") == 0) {
                set_string(&c_method, value);
            }
            else if (strcmp(key, "outdir") == 0) {
                set_string(&c_outdir, value);
            }
            else if (strcmp(key, "debug") == 0) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_debug = true;
                } else {
                    c_debug = false;
                }
            }
            else if (strcmp(key, "prefix") == 0) {
                set_string(&c_prefix, value);
            }
            else if (strcmp(key, "verbosity") == 0) {
                set_string(&c_verbosity, value);
            }
            else if (strcmp(key, "automatic_file_read") == 0) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_automatic_file_read = true;
                } else {
                    c_automatic_file_read = false;
                }
            }
            else if (strcmp(key, "write_result") == 0) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_write_result = true;
                } else {
                    c_write_result = false;
                }
            }
            else if (strcmp(key, "filetype") == 0) {
                set_string(&c_filetype, value);
            }
            else if (strcmp(key, "mpi_procs") == 0) {
                c_mpi_procs = atoi(value);
            }

//[SYSTEM]
            else if (strcmp(key, "interaction") == 0) {
                set_string(&c_interaction, value);
            }
            else if (strcmp(key, "dimension") == 0) {
                c_dimension = atoi(value);
                 got_dimension = true;
            }
            else if (strcmp(key, "celltype") == 0) {
                set_string(&c_celltype, value);
            }
            else if (strcmp(key, "nbnd") == 0) {
                c_nbnd = atoi(value);
                 got_nbnd = true;
            }
            else if (strcmp(key, "fermi_energy") == 0) {
                c_fermi_energy = atof(value);
            }
            else if (strcmp(key, "num_electrons") == 0) {
                c_num_electrons = atof(value);
            }
            else if (strcmp(key, "mu_from_n") == 0) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_mu_from_n = true;
                } else {
                    c_mu_from_n = false;
                }
            }
            else if (strcmp(key, "Temperature") == 0) {
                c_Temperature = atof(value);
            }
            else if (strcmp(key, "cutoff_energy") == 0) {
                c_cutoff_energy = atof(value);
            }
            else if (strcmp(key, "smearing") == 0) {
                c_smearing = atof(value);
            }
            else if (strcmp(key, "mixing") == 0) {
                c_mixing = atof(value);
            }
            else if (strcmp(key, "max_iters") == 0) {
                c_max_iters = atoi(value);
            }
            else if (strcmp(key, "qp_weight") == 0) {
                c_qp_weight = atof(value);
            }
            else if (strcmp(key, "recurse_level") == 0) {
                c_recurse_level = atoi(value);
            }
            else if (strcmp(key, "tol") == 0) {
                c_tol = atof(value);
            }

//[HAMILTONIAN]
            else if (strcmp(key, "hamiltonian") == 0) {
                set_string(&c_hamiltonian, value);
            }
            else if (strcmp(key, "eps_dx2y2") == 0) {
                c_eps_dx2y2 = atof(value);
            }
            else if (strcmp(key, "eps_dz2") == 0) {
                c_eps_dz2 = atof(value);
            }
            else if (strcmp(key, "eps_px") == 0) {
                c_eps_px = atof(value);
            }
            else if (strcmp(key, "eps_py") == 0) {
                c_eps_py = atof(value);
            }
            else if (strcmp(key, "eps_pz") == 0) {
                c_eps_pz = atof(value);
            }
            else if (strcmp(key, "delta_dp") == 0) {
                c_delta_dp = atof(value);
            }

//[HUBBARD]
            else if (strcmp(key, "U0") == 0) {
                c_U0 = atof(value);
            }
            else if (strcmp(key, "U1") == 0) {
                c_U1 = atof(value);
            }
            else if (strcmp(key, "J0") == 0) {
                c_J0 = atof(value);
            }
            else if (strcmp(key, "J1") == 0) {
                c_J1 = atof(value);
            }

//[MESH]
            else if (strcmp(key, "k_mesh") == 0) {
                sscanf(line, " k_mesh = %d %d %d", &c_k_mesh[0], &c_k_mesh[1], &c_k_mesh[2]);
            }
            else if (strcmp(key, "q_mesh") == 0) {
                sscanf(line, " q_mesh = %d %d %d", &c_q_mesh[0], &c_q_mesh[1], &c_q_mesh[2]);
            }
            else if (strcmp(key, "w_pts") == 0) {
                c_w_pts = atoi(value);
            }

//[CELL]

//[BRILLOUIN_ZONE]

//[BASIS]

//[BANDS]
            else if (strcmp(key, "band") == 0) {
                set_string(&c_band, value);
            }
            else if (strcmp(key, "eff_mass") == 0) {
                c_eff_mass = atof(value);
            }
            else if (strcmp(key, "t0") == 0) {
                c_t0 = atof(value);
            }
            else if (strcmp(key, "t1") == 0) {
                c_t1 = atof(value);
            }
            else if (strcmp(key, "t2") == 0) {
                c_t2 = atof(value);
            }
            else if (strcmp(key, "t3") == 0) {
                c_t3 = atof(value);
            }
            else if (strcmp(key, "t4") == 0) {
                c_t4 = atof(value);
            }
            else if (strcmp(key, "t5") == 0) {
                c_t5 = atof(value);
            }
            else if (strcmp(key, "t6") == 0) {
                c_t6 = atof(value);
            }
            else if (strcmp(key, "t7") == 0) {
                c_t7 = atof(value);
            }
            else if (strcmp(key, "t8") == 0) {
                c_t8 = atof(value);
            }
            else if (strcmp(key, "t9") == 0) {
                c_t9 = atof(value);
            }
            else if (strcmp(key, "t10") == 0) {
                c_t10 = atof(value);
            }
            else if (strcmp(key, "tz0") == 0) {
                c_tz0 = atof(value);
            }
            else if (strcmp(key, "tz1") == 0) {
                c_tz1 = atof(value);
            }
            else if (strcmp(key, "tz2") == 0) {
                c_tz2 = atof(value);
            }
            else if (strcmp(key, "tz3") == 0) {
                c_tz3 = atof(value);
            }
            else if (strcmp(key, "tz4") == 0) {
                c_tz4 = atof(value);
            }

//[SUPERCONDUCTOR]
            else if (strcmp(key, "FS_only") == 0) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_FS_only = true;
                } else {
                    c_FS_only = false;
                }
            }
            else if (strcmp(key, "num_eigenvalues_to_save") == 0) {
                c_num_eigenvalues_to_save = atoi(value);
            }
            else if (strcmp(key, "frequency_pts") == 0) {
                c_frequency_pts = atoi(value);
            }
            else if (strcmp(key, "projections") == 0) {
                set_string(&c_projections, value);
            }

//[RESPONSE]
            else if (strcmp(key, "dynamic") == 0) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_dynamic = true;
                } else {
                    c_dynamic = false;
                }
            }

//[MANY_BODY]
            else if (strcmp(key, "self_consistent") == 0) {
                strip_single_quotes(value);
                if (strcmp(value, "true") == 0) {
                    c_self_consistent = true;
                } else {
                    c_self_consistent = false;
                }
            }
            else if (strcmp(key, "impurity_solver") == 0) {
                set_string(&c_impurity_solver, value);
            }
            else if (strcmp(key, "n_cycles") == 0) {
                c_n_cycles = atoi(value);
            }
            else if (strcmp(key, "len_cycles") == 0) {
                c_len_cycles = atoi(value);
            }
            else if (strcmp(key, "warmup_cycles") == 0) {
                c_warmup_cycles = atoi(value);
            }
            // End of variable reading
            else {
                fprintf(stderr, "Error: %s:%d: unknown config variable '%s'\n",
                        config_display_name ? config_display_name : path, line_no, key);
                fprintf(stderr, "  Valid variables are listed in src/config/input_variables.py\n");
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
}

void load_c_config_from(const char *source) {
    char path[PATH_MAX];
    input_cfg_path(path, sizeof(path));
    write_input_cfg(source);
    read_c_config(path);
    // Child processes (method executables, Python and Julia methods) read the same file
    setenv("FIREFLY_CONFIG", path, 1);
}

void load_c_config() {
    load_c_config_from(NULL);
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
    //        free(c_U0[i]);
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
    ////[BASIS]
    //    free(c_state);
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
    int failed_tests[num_tests];
    int fail_count = 0;
    for (int i = 0; i < num_tests; i++) {
        if (all_tests[i])
            tests_passed++;
        else
            failed_tests[fail_count++] = (i + 1);
    }
    if (tests_passed == num_tests) {
        printcolor(GREEN, "All %d %s passed!\n", num_tests, test_name);
        return true;
    } else {
        printcolor(RED, " - %d/%d %s passed\n", tests_passed, num_tests,
                   test_name);
        for (int i = 0; i < fail_count; i++) {
            printcolor(RED, "   - Test %d failed\n", failed_tests[i]);
        }
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
