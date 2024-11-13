#pragma once
#ifndef CONFIG_H
#define CONFIG_H

#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>


// Global Variables are listed below, with their default values

// [CONTROL]
extern char c_calculation_type[50];
extern char* get_calculation_type();
extern char c_outdir[50]; char* get_outdir();
extern char c_prefix[50]; char* get_prefix();
extern char c_verbosity[50]; char* get_verbosity();

// [SYSTEM]
extern char c_interaction[50]; char* get_interaction();
extern int c_ibrav;
extern float c_mu;
extern float c_U;

// [SUPERCONDUCTOR]
extern bool c_FS_only;
extern float c_wc;
extern int c_num_eigenvalues_to_save;
extern int c_frequency_pts;

// [MESH]
extern int c_k_mesh[3];
extern int c_q_mesh[3];
extern int c_w_pts;

// [CELL]
extern float c_cell[3][3];
extern float c_brillouin_zone[3][3];
               
// [BANDS]

// INPUTS DERIVED FROM THE CONFIG ABOVE
extern float c_max_freq;
extern int c_dim;
extern char c_band_name[50]; char* get_band_name();


void get_dimensions();

void convert_to_BZ(const float cell[3][3], float brillouin_zone[3][3]);
void set_string(char *dest, const char *src);

void load_c_config(FILE *file);

#endif // CONFIG_H
