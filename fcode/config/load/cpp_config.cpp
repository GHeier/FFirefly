#include <type_traits>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include "cpp_config.hpp"
#include "c_config.h"

using namespace std;

// Global Variables are listed below

//[CONTROL]
string category;
string calculation;
string outdir;
string prefix;
string verbosity;
string datfile_in;
string datfile_out;

//[SYSTEM]
string interaction;
int dimension;
int ibrav;
int nbnd;
float fermi_energy;
float Temperature;
float onsite_U;

//[MESH]
vector<int> k_mesh(3);
vector<int> q_mesh(3);
int w_pts;

//[CELL]
vector<vector<float>> cell(3, vector<float>(3));

//[BRILLOUIN_ZONE]
vector<vector<float>> brillouin_zone(3, vector<float>(3));

//[BANDS]
vector<string> band;

vector<float> eff_mass;
vector<float> t0;
vector<float> t1;
vector<float> t2;
vector<float> t3;
vector<float> t4;
vector<float> t5;
vector<float> t6;
vector<float> t7;
vector<float> t8;
vector<float> t9;
vector<float> t10;

//[SUPERCONDUCTOR]
string method;
bool FS_only;
float bcs_cutoff_frequency;
int num_eigenvalues_to_save;
int frequency_pts;

//[RESPONSE]
bool dynamic;
// End of Global Variables 


extern "C" void load_cpp_config() {
    printf("Loading C++ Configuration File\n");
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
    nbnd = c_nbnd;
    fermi_energy = c_fermi_energy;
    Temperature = c_Temperature;
    onsite_U = c_onsite_U;

//[MESH]
    for (int i = 0; i < 3; i++) k_mesh[i] = c_k_mesh[i];
    for (int i = 0; i < 3; i++) q_mesh[i] = c_q_mesh[i];
    w_pts = c_w_pts;

//[CELL]
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) cell[i][j] = c_cell[i][j];

//[BRILLOUIN_ZONE]
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) brillouin_zone[i][j] = c_brillouin_zone[i][j];

    printf("band[0] = %s\n", c_band[0]);
//[BANDS]
    for (int i = 0; i < nbnd; i++) band.push_back(c_band[i]);
    for (int i = 0; i < nbnd; i++) eff_mass.push_back(c_eff_mass[i]);
    for (int i = 0; i < nbnd; i++) t0.push_back(c_t0[i]);
    for (int i = 0; i < nbnd; i++) t1.push_back(c_t1[i]);
    for (int i = 0; i < nbnd; i++) t2.push_back(c_t2[i]);
    for (int i = 0; i < nbnd; i++) t3.push_back(c_t3[i]);
    for (int i = 0; i < nbnd; i++) t4.push_back(c_t4[i]);
    for (int i = 0; i < nbnd; i++) t5.push_back(c_t5[i]);
    for (int i = 0; i < nbnd; i++) t6.push_back(c_t6[i]);
    for (int i = 0; i < nbnd; i++) t7.push_back(c_t7[i]);
    for (int i = 0; i < nbnd; i++) t8.push_back(c_t8[i]);
    for (int i = 0; i < nbnd; i++) t9.push_back(c_t9[i]);
    for (int i = 0; i < nbnd; i++) t10.push_back(c_t10[i]);

//[SUPERCONDUCTOR]
    method = c_method;
    FS_only = c_FS_only;
    bcs_cutoff_frequency = c_bcs_cutoff_frequency;
    num_eigenvalues_to_save = c_num_eigenvalues_to_save;
    frequency_pts = c_frequency_pts;

//[RESPONSE]
    dynamic = c_dynamic;
    // End of Global Functions 
    printf("C++ Configuration File Loaded\n");
}


//void set_global(string &a, string b) {
//    a = b;
//}

void set_global(string &a, const char* b) {
    //printf("Setting %s to %s\n", a.c_str(), b);
    a = b;
}

void set_nbnd(int nbnd_) {
    nbnd = nbnd;
}

void set_band(int n, const char* band_) {
    band[n] = band_;
}

void set_eff_mass(int n, float eff_mass_) {
    eff_mass[n] = eff_mass_;
}

void set_t0(int n, float t0_) {
    t0[n] = t0_;
}

void set_t1(int n, float t1_) {
    t1[n] = t1_;
}

