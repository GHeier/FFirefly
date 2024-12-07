/**
 * @file cpp_config.c
 *
 * @brief Configuration file, containing constants, relevant number of divisions, and functions
 * that can be called, pointing to the specific model required
 *
 * @author Griffin Heier
 */

#include <math.h>
#include <string>
#include "cpp_config.h"
#include "c_config.h"

using namespace std;

// Global Variables are listed below

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
extern string bands;
extern float eff_mass;

//[SUPERCONDUCTOR]
extern bool FS_only;
extern float bcs_cutoff_frequency;
extern int num_eigenvalues_to_save;
extern int frequency_pts;

//[RESPONSE]
extern bool dynamic;
// End of Global Variables 
int n = k_mesh[0];
int m = q_mesh[0];
int l = 5;
float wc = bcs_cutoff_frequency;
float mu = fermi_energy;
float U = onsite_U;
// Constants
float t = 1.0;
float tn = 0.0;
float tnn = 0.0;
float k_max = M_PI;

void load_cpp_config() {
    // Load the C++ configuration file

    //[CONTROL]
    category = c_category;
    calculation = c_calculation;
    outdir = c_outdir;
    prefix = c_prefix;
    verbosity = c_verbosity;
    datfile_in = c_datfile_in;
    datfile_out = c_datfile_out;

    //[SYSTEM]
    interaction = c_interaction;
    dimension = c_dimension;
    ibrav = c_ibrav;
    fermi_energy = c_fermi_energy;
    onsite_U = c_onsite_U;

    //[MESH]
    for (int i = 0; i < 3; i++) k_mesh[i] = c_k_mesh[i];
    for (int i = 0; i < 3; i++) q_mesh[i] = c_q_mesh[i];
    w_pts = c_w_pts;

    //[CELL]

    //[BANDS]
    bands = c_bands;
    eff_mass = c_eff_mass;

    //[SUPERCONDUCTOR]
    FS_only = c_FS_only;
    bcs_cutoff_frequency = c_bcs_cutoff_frequency;
    num_eigenvalues_to_save = c_num_eigenvalues_to_save;
    frequency_pts = c_frequency_pts;

    //[RESPONSE]
    dynamic = c_dynamic;
    // End of Global Functions 
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cell[i][j] = c_cell[i][j];
            brillouin_zone[i][j] = c_brillouin_zone[i][j];
        }
    }
}

void init_config(float &mu, float &U, float &t, float &tn, float &wc, float new_mu, float new_U, float new_t, float new_tn, float new_wc) {
    mu = new_mu;
    U = new_U;
    t = new_t;
    tn = new_tn;
    wc = new_wc;
}

void change_global_constant(float &a, float b) {
    a = b;
}

