#pragma once

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Global Variables are listed below, with their default values

//[CONTROL]
extern char* c_category; char* get_category();
extern char* c_calculation; char* get_calculation();
extern char* c_method; char* get_method();
extern char* c_outdir; char* get_outdir();
extern bool c_debug;
extern char* c_prefix; char* get_prefix();
extern char* c_verbosity; char* get_verbosity();
extern bool c_automatic_file_read;
extern bool c_write_result;
extern char* c_filetype; char* get_filetype();

//[SYSTEM]
extern char* c_interaction; char* get_interaction();
extern int c_dimension;
extern char* c_celltype; char* get_celltype();
extern int c_nbnd;
extern float c_fermi_energy;
extern float c_num_electrons;
extern bool c_mu_from_n;
extern float c_Temperature;
extern float c_cutoff_energy;
extern float c_smearing;
extern float c_mixing;
extern int c_max_iters;
extern float c_qp_weight;
extern int c_recurse_level;

//[HAMILTONIAN]
extern char* c_hamiltonian; char* get_hamiltonian();
extern float c_eps_dx2y2;
extern float c_eps_dz2;
extern float c_eps_px;
extern float c_eps_py;
extern float c_eps_pz;
extern float c_delta_dp;

//[HUBBARD]
extern float c_U0;
extern float c_U1;
extern float c_J0;
extern float c_J1;

//[MESH]
extern int c_k_mesh[3];
extern int c_q_mesh[3];
extern int c_w_pts;

//[CELL]
extern float c_cell[3][3];

//[BRILLOUIN_ZONE]
extern float c_brillouin_zone[3][3];

//[BASIS]
extern char** c_states; char** get_states();
extern float c_positions[50][3];

//[BANDS]
extern char* c_band; char* get_band();
extern float c_eff_mass;
extern float c_t0;
extern float c_t1;
extern float c_t2;
extern float c_t3;
extern float c_t4;
extern float c_t5;
extern float c_t6;
extern float c_t7;
extern float c_t8;
extern float c_t9;
extern float c_t10;
extern float c_tz0;
extern float c_tz1;
extern float c_tz2;
extern float c_tz3;
extern float c_tz4;

//[SUPERCONDUCTOR]
extern bool c_FS_only;
extern int c_num_eigenvalues_to_save;
extern int c_frequency_pts;
extern char* c_projections; char* get_projections();

//[RESPONSE]
extern bool c_dynamic;

//[MANY_BODY]
extern bool c_self_consistent;
extern char* c_impurity_solver; char* get_impurity_solver();
// End of Global Variables

void get_dimensions();

void cell_to_BZ(float ucell[3][3], float (*bz_matrix)[3]);
void make_lowercase(char *str);
void set_string(char **dest, const char *src);
void set_section(char *dest, const char *src);
void load_default_band_values();

void make_save_file();
void load_c_config();
extern void unload_c_config();

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    RESET = 0,
    RED = 31,
    GREEN = 32,
    YELLOW = 33,
    BLUE = 34,
    MAGENTA = 35,
    CYAN = 36,
    WHITE = 37
} Color;

void read_c_config(const char *path);
void printcolor(Color color, const char *format, ...);
bool print_test_results(bool all_tests[], int num_tests, const char *test_name);
void load_cpp_config();

#ifdef __cplusplus
}
#endif
int mkdir_p(const char *path, mode_t mode);
