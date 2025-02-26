#pragma once

#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

// Global Variables are listed below, with their default values

//[CONTROL]
extern char* c_category; char* get_category();
extern char* c_calculation; char* get_calculation();
extern char* c_outdir; char* get_outdir();
extern char* c_prefix; char* get_prefix();
extern char* c_verbosity; char* get_verbosity();
extern char* c_input_data_file; char* get_input_data_file();
extern char* c_output_data_file; char* get_output_data_file();

//[SYSTEM]
extern char* c_interaction; char* get_interaction();
extern int c_dimension;
extern int c_ibrav;
extern int c_nbnd;
extern float c_fermi_energy;
extern float c_Temperature;
extern float c_onsite_U;

//[MESH]
extern int c_k_mesh[3];
extern int c_q_mesh[3];
extern int c_w_pts;

//[CELL]
extern float c_cell[3][3];

//[BRILLOUIN_ZONE]
extern float c_brillouin_zone[3][3];

//[BANDS]
extern char** c_band; char** get_band();
extern float c_eff_mass[50];
extern float c_t0[50];
extern float c_t1[50];
extern float c_t2[50];
extern float c_t3[50];
extern float c_t4[50];
extern float c_t5[50];
extern float c_t6[50];
extern float c_t7[50];
extern float c_t8[50];
extern float c_t9[50];
extern float c_t10[50];

//[SUPERCONDUCTOR]
extern char* c_method; char* get_method();
extern bool c_FS_only;
extern float c_bcs_cutoff_frequency;
extern int c_num_eigenvalues_to_save;
extern int c_frequency_pts;
extern char* c_projections; char* get_projections();

//[RESPONSE]
extern bool c_dynamic;
// End of Global Variables

void get_dimensions();

void cell_to_BZ(float ucell[3][3], float (*bz_matrix)[3]);
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

void read_c_config(const char* path);
void printcolor(Color color, const char* format, ...);
bool print_test_results(bool all_tests[], int num_tests, const char* test_name);

#ifdef __cplusplus
}
#endif

