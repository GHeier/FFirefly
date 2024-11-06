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
#include <cassert>
#define assertm(exp, msg) assert(((void)msg, exp))
#include "cfg.h"

using namespace std;

// Global Variables
int n = 20; // Number of k points
int m = 20; // Number of chi points
int l = 5; // Number of BCS frequency points
int w_pts = 11; // Number of matsubara frequency points
float max_freq = 10.0; // Maximum frequency for matsubara
int dim = 3; // Number of dimensions)
string potential_name = "scalapino";
string band_name = "sphere";
int num_eigenvalues_to_save = 1;
bool FS_only = true;
               
// Constants
float t = 1.0;
float tn = 0.0;
float tnn = 0.0;
float U = 4.0;
float k_max = M_PI;
float mu = 1.0;
float wc = 0.05;


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

// Gaussian integration constants
float weights_0th[1] = {2.0}; float * w0 = weights_0th;
float weights_1st[2] = {1.0, 1.0}; float * w1 = weights_1st;
float weights_2nd[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0}; float * w2 = weights_2nd;
float weights_3rd[4] = {0.347855, 0.652145, 0.652145, 0.347855}; float * w3 = weights_3rd;
float weights_4th[5] = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927}; float * w4 = weights_4th;
float *weights[5] = {w0, w1, w2, w3, w4};

float points_0th[1] = {0}; float *p0 = points_0th;
float points_1st[2] = {-1/pow(3,0.5), 1/pow(3,0.5)}; float *p1 = points_1st;
float points_2nd[3] = {-pow(3/5,0.5), 0, pow(3/5,0.5)}; float *p2 = points_2nd;
float points_3rd[4] = {-0.861136, -0.339981, 0.339981, 0.861136}; float *p3 = points_3rd;
float points_4th[5] = {-0.90618, -0.538469, 0, 0.538469, 0.90618}; float *p4 = points_4th;

float *points[5] = {p0, p1, p2, p3, p4};
