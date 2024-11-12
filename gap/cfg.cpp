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
#include <sstream>
#include <fstream>
#include "cfg.h"

using namespace std;

Config cfg; // Global configuration object to be used by all languages

// Global Variables are listed below, with their default values
// Config Parameters
string calculation_type = "projection";
string outdir = "./output";
string prefix = "sample";
string verbosity = "high";

string potential_name = "scalapino";
int ibrav = 0;

// Mesh parameters
int n = 20; // Number of k points
int m = 20; // Number of chi points
int l = 5; // Number of BCS frequency points
int w_pts = 11; // Number of matsubara frequency points
float max_freq = 10.0; // Maximum frequency for matsubara

float cell[3][3] = {{1.0, 0.0, 0.0},
                     {0.0, 1.0, 0.0},
                     {0.0, 0.0, 1.0}};

float brillouin_zone[3][3] = {{2.0 * M_PI, 0.0, 0.0},
                               {0.0, 2.0 * M_PI, 0.0},
                               {0.0, 0.0, 2.0 * M_PI}};

int dim = 3; // Number of dimensions)
string band_name = "sphere";
int num_eigenvalues_to_save = 1;
bool FS_only = true;
               
// Constants
float t = 1.0;
float tn = 0.0;
float tnn = 0.0;
float U = 4.0;
float k_max = M_PI;
float mu = 0.0;
float wc = 0.05;

using namespace std;

extern "C" {
#include <string.h>
    void convert_to_BZ(const float cell[3][3], float brillouin_zone[3][3]);
    void set_string(char *dest, const char *src, int size) {
        strncpy(dest, src, size - 1);
        dest[size - 1] = '\0';  // Ensure null-termination
    }
}

// Modify the potential_name field
void fill_cfg(Config &c) {
    c.n = n;
    c.m = m;
    c.l = l;
    c.w_pts = w_pts;
    c.max_freq = max_freq;
    c.dim = dim;
    set_string(c.potential_name, potential_name.c_str(), 50);
    set_string(c.band_name, band_name.c_str(), 50);
    set_string(c.calculation_type, calculation_type.c_str(), 50);
    c.num_eigenvalues_to_save = num_eigenvalues_to_save;
    c.FS_only = FS_only;

    c.mu = mu;
    c.k_max = k_max;

    c.t = t;
    c.tn = tn;
    c.tnn = tnn;
    c.U = U;
    c.wc = wc;
}

void load_config(ifstream &file) {
    string line;
    string a, b, c;
    string section;
    int row = 0;
    while (getline(file, line)) {
        cout << "Line: " << line << endl;
        if (line.find("=") != string::npos) {
            string key = line.substr(0, line.find("="));
            string value = line.substr(line.find("=") + 1);
            
            if (line == "[CELL]") {
                section = "CELL";
                continue;  // Skip to the next line
            }

            // If we're in the [CELL] section, read matrix values
            if (section == "CELL" && !line.empty()) {
                istringstream iss(line);
                if (row < 3) {
                    iss >> cell[row][0] >> cell[row][1] >> cell[row][2];
                    row++;
                }

                // Stop reading after filling 3 rows
                if (row == 3) break;
            }

            if (key == "calculation") {
              calculation_type = value;
            } else if (key == "outdir") {
                outdir = value;
            } else if (key == "potential") {
                potential_name = value;
            } else if (key == "FS_only") {
                FS_only = stoi(value);
            } else if (key == "surface_pnts") {
                n = stoi(value);
            } else if (key == "potential_pnts") {
                m = stoi(value);
            } else if (key == "frequency_pts") {
                l = stoi(value);
            } else {
                cout << "Invalid key: " << key << endl;
                exit(1);
            }
        }
    }
    convert_to_BZ(cell, brillouin_zone);
    fill_cfg(cfg);
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

