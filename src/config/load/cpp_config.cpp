#include <type_traits>
#include <iostream>
#include <filesystem>
#include <math.h>
#include <string>
#include <vector>
#include "cpp_config.hpp"
#include "c_config.h"

using namespace std;
namespace fs = std::filesystem;
// Global Variables are listed below

//[CONTROL]
string category;
string calculation;
string method;
string outdir;
string indir;
string prefix;
string verbosity;
bool automatic_file_read;
bool write_result;
string filetype;

//[SYSTEM]
string hamiltonian;
string interaction;
int dimension;
string celltype;
int nbnd;
int nstates;
float fermi_energy;
float num_electrons;
float Temperature;
float onsite_U;
float cutoff_energy;
float smearing;
float mixing;
int max_iters;

//[MESH]
vector<int> k_mesh(3);
vector<int> q_mesh(3);
int w_pts;

//[CELL]
vector<vector<float>> cell(3, vector<float>(3));

//[BRILLOUIN_ZONE]
vector<vector<float>> brillouin_zone(3, vector<float>(3));

//[BASIS]
vector<string> states;
vector<vector<float>> positions(50, vector<float>(3));

//[BANDS]
string band;
float eff_mass;
float t0;
float t1;
float t2;
float t3;
float t4;
float t5;
float t6;
float t7;
float t8;
float t9;
float t10;

//[SUPERCONDUCTOR]
bool FS_only;
int num_eigenvalues_to_save;
int frequency_pts;
string projections;

//[RESPONSE]
bool dynamic;

//[MANY_BODY]
bool self_consistent;
// End of Global Variables

// Track if config has been loaded
static bool cpp_config_loaded = false;

extern "C" void load_cpp_config() {
    cpp_config_loaded = true;
    // Load the C++ configuration file

//[CONTROL]
    category = c_category;
    calculation = c_calculation;
    method = c_method;
    outdir = c_outdir;
    indir = c_indir;
    prefix = c_prefix;
    verbosity = c_verbosity;
    automatic_file_read = c_automatic_file_read;
    write_result = c_write_result;
    filetype = c_filetype;

//[SYSTEM]
    hamiltonian = c_hamiltonian;
    interaction = c_interaction;
    dimension = c_dimension;
    celltype = c_celltype;
    nbnd = c_nbnd;
    nstates = c_nstates;
    fermi_energy = c_fermi_energy;
    num_electrons = c_num_electrons;
    Temperature = c_Temperature;
    onsite_U = c_onsite_U;
    cutoff_energy = c_cutoff_energy;
    smearing = c_smearing;
    mixing = c_mixing;
    max_iters = c_max_iters;

//[MESH]
    for (int i = 0; i < 3; i++) k_mesh[i] = c_k_mesh[i];
    for (int i = 0; i < 3; i++) q_mesh[i] = c_q_mesh[i];
    w_pts = c_w_pts;

//[CELL]
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) cell[i][j] = c_cell[i][j];

//[BRILLOUIN_ZONE]
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) brillouin_zone[i][j] = c_brillouin_zone[i][j];

//[BASIS]
    for (int i = 0; i < nstates; i++) states.push_back(c_states[i]);
    for (int i = 0; i < nstates; i++) for (int j = 0; j < 3; j++) positions[i][j] = c_positions[i][j];

//[BANDS]
    band = c_band;
    eff_mass = c_eff_mass;
    t0 = c_t0;
    t1 = c_t1;
    t2 = c_t2;
    t3 = c_t3;
    t4 = c_t4;
    t5 = c_t5;
    t6 = c_t6;
    t7 = c_t7;
    t8 = c_t8;
    t9 = c_t9;
    t10 = c_t10;

//[SUPERCONDUCTOR]
    FS_only = c_FS_only;
    num_eigenvalues_to_save = c_num_eigenvalues_to_save;
    frequency_pts = c_frequency_pts;
    projections = c_projections;

//[RESPONSE]
    dynamic = c_dynamic;

//[MANY_BODY]
    self_consistent = c_self_consistent;
    // End of Global Functions 
    if (!isDirectoryExisting(outdir)) {
        if (fs::create_directory(outdir)) {
            std::cout << "Directory " << outdir << " created successfully.\n";
        } else {
            std::cout << "Failed to create directory " << outdir << "\n";
        }
    }
}


//void set_global(string &a, string b) {
//    a = b;
//}

void set_global(string &a, const char* b) {
    //printf("Setting %s to %s\n", a.c_str(), b);
    a = b;
}

void set_global(vector<int> &a, vector<int> b) {
    printv("Setting mesh to $d $d $d\n", b[0], b[1], b[2]);
    a = b;
}

void set_nbnd(int nbnd_) {
    nbnd = nbnd;
}

void set_band(int n, const char* band_) {
    band = band_;
}

void set_eff_mass(int n, float eff_mass_) {
    eff_mass = eff_mass_;
}

void set_t0(int n, float t0_) {
    t0 = t0_;
}

void set_t1(int n, float t1_) {
    t1 = t1_;
}

void read_c_config_wrapper(string path) {
    read_c_config(path.c_str());
    //load_cpp_config();
}

bool isDirectoryExisting(const std::string& path) {
    std::filesystem::path dirPath(path);
    return std::filesystem::exists(dirPath) && std::filesystem::is_directory(dirPath);
}

void ensure_cpp_config_loaded() {
    if (!cpp_config_loaded) {
        read_c_config("/home/g/Research/FFirefly/build/bin/input.cfg");
        load_cpp_config();
    }
}
