#pragma once
#ifndef C_CONFIG_H
#define C_CONFIG_H

#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>


// Global Variables are listed below, with their default values

// [CONTROL]
extern char c_category[50];
extern char* get_category();
extern char c_outdir[50]; char* get_outdir();
extern char c_prefix[50]; char* get_prefix();
extern char c_verbosity[50]; char* get_verbosity();
extern char c_datfile_in[50]; char* get_datfile_in();
extern char c_datfile_out[50]; char* get_datfile_out();

// [SYSTEM]
extern char c_interaction[50]; char* get_interaction();
extern int c_ibrav;
extern float c_fermi_energy;
extern float c_onsite_U;

// [SUPERCONDUCTOR]
extern bool c_FS_only;
extern float c_bcs_cutoff_frequency;
extern int c_num_eigenvalues_to_save;
extern int c_frequency_pts;

// [RESPONSE]
extern char c_calculation[50]; char* get_calculation();
extern bool dynamic;

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

//extern void cell_to_BZ(const float cell[3][3], float brillouin_zone[3][3]);
void set_string(char *dest, const char *src);

extern void load_c_config();
void foo();

#endif // CONFIG_H
