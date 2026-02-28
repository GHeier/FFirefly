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
extern bool debug;
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
extern float fermi_energy;
extern float num_electrons;
extern bool mu_from_n;
extern float Temperature;
extern float cutoff_energy;
extern float smearing;
extern float mixing;
extern int max_iters;
extern float qp_weight;

//[HAMILTONIAN]
extern string hamiltonian;

//[HUBBARD]
extern float U0;
extern float U1;
extern float J0;
extern float J1;

//[MESH]
extern vector<int> k_mesh;
extern vector<int> q_mesh;
extern int w_pts;

//[CELL]
extern vector<vector<float>> cell;

//[BRILLOUIN_ZONE]
extern vector<vector<float>> brillouin_zone;

//[BASIS]
extern vector<string> states;
extern vector<vector<float>> positions;

//[BANDS]
extern string band;
extern float eff_mass;
extern float t0;
extern float t1;
extern float t2;
extern float t3;
extern float t4;
extern float t5;
extern float t6;
extern float t7;
extern float t8;
extern float t9;
extern float t10;

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

// Ensure config is loaded (call this before using any C++ objects from Python/Julia)
void ensure_cpp_config_loaded();

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

// Run executable with config file piped to stdin
int run_with_config(const std::string& executable, const std::string& config_file);

// Helper to get executable location
std::string get_loc();

// Run method executables (used by category node.cpp files)
int run_cpp_method(const std::string& method_name);
int run_python_method(const std::string& method_name);
int run_julia_method(const std::string& method_name);

int run_cpp_test(const std::string& test_name);
int run_python_test(const std::string& test_name);
int run_julia_test(const std::string& test_name);

// Alternative shell-based execution methods (method2)
int run_python_method2(const std::string& method_name);
int run_julia_method2(const std::string& method_name);
