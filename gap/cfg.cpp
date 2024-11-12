/**
 * @file cfg.cpp
 *
 * @brief Configuration file, containing constants, relevant number of divisions, and functions
 * that can be called, pointing to the specific model required
 *
 * @author Griffin Heier
 */

#include <math.h>
#include <string>
#include <fstream>
#include "cfg.h"

using namespace std;

// Global Variables
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

void load_config(ifstream &file) {
    string line;
    string a, b, c;
    while (getline(file, line)) {
        cout << "Line: " << line << endl;
        if (line.find("=") != string::npos) {
            string key = line.substr(0, line.find("="));
            string value = line.substr(line.find("=") + 1);

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

