#pragma once

#include <string>
#include <vector>

using namespace std;

// Global Variables are listed below'

//[CONTROL]
extern string category;
extern string calculation;
extern string method;
extern string outdir;
extern string indir;
extern string prefix;
extern string verbosity;
extern bool automatic_file_read;
extern bool write_result;
extern string filetype;

//[SYSTEM]
extern string interaction;
extern int dimension;
extern string celltype;
extern int nbnd;
extern int natoms;
extern float fermi_energy;
extern float Temperature;
extern float onsite_U;
extern float cutoff_energy;

//[MESH]
extern vector<int> k_mesh;
extern vector<int> q_mesh;
extern int w_pts;

//[CELL]
extern vector<vector<float>> cell;

//[BRILLOUIN_ZONE]
extern vector<vector<float>> brillouin_zone;

//[ATOMIC_POSITIONS]
extern vector<string> atom;

extern vector<vector<float>> position;


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
extern bool FS_only;
extern int num_eigenvalues_to_save;
extern int frequency_pts;
extern string projections;

//[RESPONSE]
extern bool dynamic;

//[MANY_BODY]
extern bool self_consistent;
// End of Global Variables


extern "C" void load_cpp_config();

template <typename T>
void set_global(T &a, T b) {
    a = b;
}

inline void set_global(float &a, float b) {
    a = b;
}
//
//void set_global(string &a, string b);
void set_global(string &a, const char* b);
void set_global(vector<int> &a, vector<int> b);

extern "C" {
    void set_nbnd(int nbnd_);
    void set_band(int n, const char* band_);
    void set_eff_mass(int n, float eff_mass_);
    void set_t0(int n, float t0_);
    void set_t1(int n, float t1_);
}

//inline void set_global(int &a, int b) {
//    a = b;
//}
//
//inline void set_global(bool &a, bool b) {
//    a = b;
//}

template <typename... Args>
void printv(const std::string& format, Args... args) {
    if (verbosity == "high") printf(format.c_str(), args...);
}
void read_c_config_wrapper(string path);
bool isDirectoryExisting(const std::string& path);
