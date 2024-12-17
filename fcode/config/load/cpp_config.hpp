#pragma once

#include <string>
#include <vector>

using namespace std;
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
extern float Temperature;
extern float onsite_U;

//[MESH]
extern vector<int> k_mesh;
extern vector<int> q_mesh;
extern int w_pts;

//[CELL]
extern vector<vector<float>> cell;
extern vector<vector<float>> brillouin_zone;

//[BANDS]
extern vector<string> band;

extern vector<float> eff_mass;
extern vector<float> t0;
extern vector<float> t1;
extern vector<float> t2;
extern vector<float> t3;
extern vector<float> t4;
extern vector<float> t5;
extern vector<float> t6;
extern vector<float> t7;
extern vector<float> t8;
extern vector<float> t9;
extern vector<float> t10;

//[SUPERCONDUCTOR]
extern string method;
extern bool FS_only;
extern float bcs_cutoff_frequency;
extern int num_eigenvalues_to_save;
extern int frequency_pts;

//[RESPONSE]
extern bool dynamic;
// End of Global Variables

void load_cpp_config();
void set_global_constant(float &a, float b);

