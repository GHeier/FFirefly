/**
 * @file cfg.c
 *
 * @brief Configuration file, containing constants, relevant number of divisions, and functions
 * that can be called, pointing to the specific model required
 *
 * @author Griffin Heier
 */

#include <math.h>
#include <string>
#include "cfg.h"
#include "config.h"

using namespace std;

// Global Variables are listed below, with their default values
// Config Parameters
string calculation(c_calculation);
string category(c_category);
string outdir(c_outdir);
string prefix(c_prefix);
string verbosity(c_verbosity);

string potential_name(c_interaction);
int ibrav = c_ibrav;
int n = c_k_mesh[0];
int m = c_q_mesh[0];
int l = 5;
int w_pts = c_w_pts;
float max_freq = c_max_freq;

float wc = c_bcs_cutoff_frequency;
float mu = c_fermi_energy;
float U = c_onsite_U;

float cell[3][3]; 
float brillouin_zone[3][3]; 


int dim = 3; // Number of dimensions)
string band_name = "sphere";
int num_eigenvalues_to_save = 1;
bool FS_only = true;
               
// Constants
float t = 1.0;
float tn = 0.0;
float tnn = 0.0;
float k_max = M_PI;

void load_config() {
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

