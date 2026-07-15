#include <type_traits>
#include <iostream>
#include <filesystem>
#include <math.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <sys/wait.h>
#include <linux/limits.h>
#include <unistd.h>
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
bool debug;
string prefix;
string verbosity;
bool automatic_file_read;
bool write_result;
string filetype;

//[SYSTEM]
string interaction;
int dimension;
string celltype;
int nbnd;
float fermi_energy;
float num_electrons;
bool mu_from_n;
float Temperature;
float cutoff_energy;
float smearing;
float mixing;
int max_iters;
float qp_weight;
int recurse_level;

//[HAMILTONIAN]
string hamiltonian;
float eps_dx2y2;
float eps_dz2;
float eps_px;
float eps_py;
float eps_pz;
float delta_dp;

//[HUBBARD]
float U0;
float U1;
float J0;
float J1;

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
float tz0;
float tz1;
float tz2;
float tz3;
float tz4;

//[SUPERCONDUCTOR]
bool FS_only;
int num_eigenvalues_to_save;
int frequency_pts;
string projections;

//[RESPONSE]
bool dynamic;

//[MANY_BODY]
bool self_consistent;
string impurity_solver;
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
    debug = c_debug;
    prefix = c_prefix;
    verbosity = c_verbosity;
    automatic_file_read = c_automatic_file_read;
    write_result = c_write_result;
    filetype = c_filetype;

//[SYSTEM]
    interaction = c_interaction;
    dimension = c_dimension;
    celltype = c_celltype;
    nbnd = c_nbnd;
    fermi_energy = c_fermi_energy;
    num_electrons = c_num_electrons;
    mu_from_n = c_mu_from_n;
    Temperature = c_Temperature;
    cutoff_energy = c_cutoff_energy;
    smearing = c_smearing;
    mixing = c_mixing;
    max_iters = c_max_iters;
    qp_weight = c_qp_weight;
    recurse_level = c_recurse_level;

//[HAMILTONIAN]
    hamiltonian = c_hamiltonian;
    eps_dx2y2 = c_eps_dx2y2;
    eps_dz2 = c_eps_dz2;
    eps_px = c_eps_px;
    eps_py = c_eps_py;
    eps_pz = c_eps_pz;
    delta_dp = c_delta_dp;

//[HUBBARD]
    U0 = c_U0;
    U1 = c_U1;
    J0 = c_J0;
    J1 = c_J1;

//[MESH]
    for (int i = 0; i < 3; i++) k_mesh[i] = c_k_mesh[i];
    for (int i = 0; i < 3; i++) q_mesh[i] = c_q_mesh[i];
    w_pts = c_w_pts;

//[CELL]
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) cell[i][j] = c_cell[i][j];

//[BRILLOUIN_ZONE]
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) brillouin_zone[i][j] = c_brillouin_zone[i][j];

//[BASIS]
    for (int i = 0; i < nbnd; i++) states.push_back(c_states[i]);
    for (int i = 0; i < nbnd; i++) for (int j = 0; j < 3; j++) positions[i][j] = c_positions[i][j];

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
    tz0 = c_tz0;
    tz1 = c_tz1;
    tz2 = c_tz2;
    tz3 = c_tz3;
    tz4 = c_tz4;

//[SUPERCONDUCTOR]
    FS_only = c_FS_only;
    num_eigenvalues_to_save = c_num_eigenvalues_to_save;
    frequency_pts = c_frequency_pts;
    projections = c_projections;

//[RESPONSE]
    dynamic = c_dynamic;

//[MANY_BODY]
    self_consistent = c_self_consistent;
    impurity_solver = c_impurity_solver;
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
        read_c_config((get_loc() + "input.cfg").c_str());
        load_cpp_config();
    }
}

int run_with_config(const std::string& executable, const std::string& config_file) {
    // Verify the executable exists
    if (!fs::exists(executable)) {
        std::cerr << "Error: Executable '" << executable << "' not found.\n";
        return -1;
    }

    // Verify the config file exists
    if (!fs::exists(config_file)) {
        std::cerr << "Error: Config file '" << config_file << "' not found.\n";
        return -1;
    }

    // Build the command: executable < config_file
    std::string command = executable + " < " + config_file;

    //if (verbosity == "high") {
    //    std::cout << "Running: " << command << std::endl;
    //}

    // Execute the command and return the exit code
    int result = std::system(command.c_str());

    if (result == -1) {
        std::cerr << "Error: Failed to execute command.\n";
        return -1;
    }

    // Extract the actual exit status
    if (WIFEXITED(result)) {
        return WEXITSTATUS(result);
    } else {
        std::cerr << "Error: Process did not terminate normally.\n";
        return -1;
    }
}

string get_loc() {
    char path[PATH_MAX]; // Buffer to hold the executable path

    // Read the symbolic link for the executable
    ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
    if (len == -1) {
        return "./build/bin/";  // Fallback
    }
    path[len] = '\0';

    // Find the last '/' to get directory
    string exe_path(path);
    size_t last_slash = exe_path.find_last_of('/');
    if (last_slash != string::npos) {
        return exe_path.substr(0, last_slash + 1);  // Include the trailing slash
    }
    return "./";
}

int run_cpp_method(const string& method_name) {
    string loc = get_loc();
    string exe = loc + category + "_" + calculation + "_" + method_name + ".exe";
    string config_path = loc + "input.cfg";
    return run_with_config(exe, config_path);
}

// Helper function to run a command and pipe its output to stdout
static int run_and_pipe_output(const string& cmd) {
    string full_cmd = cmd + " 2>&1";  // Redirect stderr to stdout
    FILE* pipe = popen(full_cmd.c_str(), "r");
    if (!pipe) {
        fprintf(stderr, "Error: popen() failed for command: %s\n", cmd.c_str());
        return -1;
    }

    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        printf("%s", buffer);
        fflush(stdout);  // Ensure output is immediately visible
    }

    int status = pclose(pipe);
    return WEXITSTATUS(status);
}

int run_python_method(const string& method_name) {
    string loc = get_loc();
    string config_path = loc + "input.cfg";
    // Go up from build/bin/ to project root, then to src/
    size_t build_pos = loc.find("/build/bin/");
    if (build_pos != string::npos) {
        // Found build/bin/, go to project root
        string src_dir = loc.substr(0, build_pos) + "/src/";
        string script_path = src_dir + category + "/" + calculation + "/" + method_name + "/run.py";
        string exe = "python3 " + script_path;
        int result = std::system(exe.c_str());
        return result;
    }
    // Fallback: try simple /bin/ pattern (for non-build locations)
    size_t bin_pos = loc.find("/bin/");
    if (bin_pos != string::npos) {
        string src_dir = loc.substr(0, bin_pos) + "/src/";
        string script_path = src_dir + category + "/" + calculation + "/" + method_name + "/run.py";
        string exe = "python3 " + script_path;
        int result = std::system(exe.c_str());
        return result;
    }
    return -1;  // Error: couldn't find src directory
}

int run_julia_method(const string& method_name) {
    string loc = get_loc();
    string config_path = loc + "input.cfg";
    // Go up from build/bin/ to project root, then to src/
    size_t build_pos = loc.find("/build/bin/");
    if (build_pos != string::npos) {
        // Found build/bin/, go to project root
        string src_dir = loc.substr(0, build_pos) + "/src/";
        string script_path = src_dir + category + "/" + calculation + "/" + method_name + "/run.jl";
        string exe = "julia " + script_path;
        int result = std::system(exe.c_str());
        return result;
    }
    // Fallback: try simple /bin/ pattern (for non-build locations)
    size_t bin_pos = loc.find("/bin/");
    if (bin_pos != string::npos) {
        string src_dir = loc.substr(0, bin_pos) + "/src/";
        string script_path = src_dir + category + "/" + calculation + "/" + method_name + "/run.jl";
        string exe = "julia " + script_path;
        int result = std::system(exe.c_str());
        return result;
    }
    return -1;  // Error: couldn't find src directory
}

int run_cpp_test(const string& test_name) {
    string loc = get_loc();
    string exe = loc + category + "_" + calculation + "_" + test_name + ".exe";
    string config_path = loc + "input.cfg";
    return run_with_config(exe, config_path);
}

int run_python_test(const string& test_name) {
    string loc = get_loc();
    string config_path = loc + "input.cfg";
    // Go up from build/bin/ to project root, then to src/
    size_t build_pos = loc.find("/build/bin/");
    if (build_pos != string::npos) {
        // Found build/bin/, go to project root
        string src_dir = loc.substr(0, build_pos) + "/src/";
        string script_path = src_dir + category + "/" + calculation + "/" + test_name + "/tests/test.py";
        string exe = "python3 " + script_path;
        int result = std::system(exe.c_str());
        return result;
        //return run_with_config(exe, config_path);
    }
    // Fallback: try simple /bin/ pattern (for non-build locations)
    size_t bin_pos = loc.find("/bin/");
    if (bin_pos != string::npos) {
        string src_dir = loc.substr(0, bin_pos) + "/src/";
        string script_path = src_dir + category + "/" + calculation + "/" + test_name + "/tests/test.py";
        string exe = "python3 " + script_path;
        int result = std::system(exe.c_str());
        return result;
        //return run_with_config(exe, config_path);
    }
    return -1;  // Error: couldn't find src directory
}

int run_julia_test(const string& test_name) {
    string loc = get_loc();
    string config_path = loc + "input.cfg";
    // Go up from build/bin/ to project root, then to src/
    size_t build_pos = loc.find("/build/bin/");
    if (build_pos != string::npos) {
        // Found build/bin/, go to project root
        string src_dir = loc.substr(0, build_pos) + "/src/";
        string script_path = src_dir + category + "/" + calculation + "/" + test_name + "/tests/test.jl";
        string exe = "julia " + script_path;
        int result = std::system(exe.c_str());
        return result;
        //return run_with_config(exe, config_path);
    }
    // Fallback: try simple /bin/ pattern (for non-build locations)
    size_t bin_pos = loc.find("/bin/");
    if (bin_pos != string::npos) {
        string src_dir = loc.substr(0, bin_pos) + "/src/";
        string script_path = src_dir + category + "/" + calculation + "/" + test_name + "tests/test.jl";
        string exe = "julia " + script_path;
        int result = std::system(exe.c_str());
        return result;
        //return run_with_config(exe, config_path);
    }
    return -1;  // Error: couldn't find src directory
}

// Method2 implementations - pure shell-based execution
int run_python_method2(const string& method_name) {
    string loc = get_loc();
    string src_dir;

    // Go up from build/bin/ to project root, then to src/
    size_t build_pos = loc.find("/build/bin/");
    if (build_pos != string::npos) {
        src_dir = loc.substr(0, build_pos) + "/src/";
    } else {
        // Fallback: try simple /bin/ pattern (for non-build locations)
        size_t bin_pos = loc.find("/bin/");
        if (bin_pos != string::npos) {
            src_dir = loc.substr(0, bin_pos) + "/src/";
        } else {
            return -1;  // Error: couldn't find src directory
        }
    }

    string script_path = src_dir + category + "/" + calculation + "/" + method_name + "/run.py";
    string command = "python3 " + script_path;

    if (verbosity == "high") {
        std::cout << "Running (method2): " << command << std::endl;
    }

    int result = std::system(command.c_str());
    if (result == -1) {
        std::cerr << "Error: Failed to execute Python script.\n";
        return -1;
    }

    if (WIFEXITED(result)) {
        return WEXITSTATUS(result);
    } else {
        std::cerr << "Error: Python process did not terminate normally.\n";
        return -1;
    }
}

int run_julia_method2(const string& method_name) {
    string loc = get_loc();
    string src_dir;

    // Go up from build/bin/ to project root, then to src/
    size_t build_pos = loc.find("/build/bin/");
    if (build_pos != string::npos) {
        src_dir = loc.substr(0, build_pos) + "/src/";
    } else {
        // Fallback: try simple /bin/ pattern (for non-build locations)
        size_t bin_pos = loc.find("/bin/");
        if (bin_pos != string::npos) {
            src_dir = loc.substr(0, bin_pos) + "/src/";
        } else {
            return -1;  // Error: couldn't find src directory
        }
    }

    string script_path = src_dir + category + "/" + calculation + "/" + method_name + "/run.jl";
    string command = "julia " + script_path;

    if (verbosity == "high") {
        std::cout << "Running (method2): " << command << std::endl;
    }

    int result = std::system(command.c_str());
    if (result == -1) {
        std::cerr << "Error: Failed to execute Julia script.\n";
        return -1;
    }

    if (WIFEXITED(result)) {
        return WEXITSTATUS(result);
    } else {
        std::cerr << "Error: Julia process did not terminate normally.\n";
        return -1;
    }
}
