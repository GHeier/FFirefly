#pragma once
#ifndef CPP_CONFIG_H
#define CPP_CONFIG_H

#include <string>
#include <unordered_map>
#include "../../objects/vec.h"

using namespace std;
extern int n;
extern int m;
extern int l;
extern float mu;
extern float k_max;

extern float t;
extern float tn;
extern float tnn;
extern float U;
extern float wc;
// Global Variables are listed below'

//[CONTROL]
extern string category;
extern string calculation;
extern string outdir;
extern string prefix;
extern string verbosity;
extern string datfile_in;
extern string datfile_out;

//[SYSTEM]
extern string interaction;
extern int dimension;
extern int ibrav;
extern int nbnd;
extern float fermi_energy;
extern float onsite_U;

//[MESH]
extern int k_mesh[3];
extern int q_mesh[3];
extern int w_pts;

//[CELL]
extern float cell[3][3];
extern float brillouin_zone[3][3];

//[BANDS]
extern string band[50];

extern float eff_mass[50];

//[SUPERCONDUCTOR]
extern bool FS_only;
extern float bcs_cutoff_frequency;
extern int num_eigenvalues_to_save;
extern int frequency_pts;

//[RESPONSE]
extern bool dynamic;
// End of Global Variables

void load_cpp_config();
void init_config(float &mu, float &U, float &t, float &tn, float &w_D, float new_mu, float new_U, float new_t, float new_tn, float new_w_D);
void change_global_constant(float &a, float b);

#endif
