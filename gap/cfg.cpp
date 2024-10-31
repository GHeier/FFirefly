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
#include <unordered_map>
#include "vec.h"
#include "cfg.h"
#include "band_structure.h"
#include "potential.h"

using namespace std;

// Global Variables
int n = 20; // Number of k points
int s_div = (dim == 3) ? 40 : 300; // Number of integral surface divisions
int s_pts = (dim == 3) ? 50 : 1000; // Number of integral surfaces
int m = 20; // Number of chi points
int l = 5; // Number of frequency points
int dim = 3; // Number of dimensions
string potential_name = "phonon_coulomb";
string band_name = "simple_cubic";
int num_eigenvalues_to_save = 1;
               
// Constants
float ε = 55.263494;
float e = 1;
float t = 1.0;
float tn = 0.0;
float tnn = 0.0;
float U = 4.0;
float k_max = M_PI;
float mu = 1.0;
float wc = 0.5;
float N = 1;
float M = 15.999;
float C = 1;
float Vol = 1;

// Debugging DOS REMEMBER TO FIX HEADER WHEN U COMMENT OUT
// float DOS = 1

Vec c(5000, 4000, 0); // lntdal ultsonic wave velocity ([001], [110], 0)


void init_config(float &mu, float &U, float &t, float &tn, float &w_D, float new_mu, float new_U, float new_t, float new_tn, float new_w_D) {
    mu = new_mu;
    U = new_U;
    t = new_t;
    tn = new_tn;
    w_D = new_w_D;
}

void change_global_constant(float &a, float b) {
    a = b;
}

// Energy band functions
float epsilon(const Vec k) {
    if (band_name == "simple_cubic_layered")
        return epsilon_SC_layered(k);
    if (band_name == "simple_cubic")
        return epsilon_SC(k, t, tn);
    if (band_name == "sphere")
        return epsilon_sphere(k);
    else {
        cout << "Unknown Band structure: " << band_name << endl;
        exit(1);
        return 0;
    }
}

// Difference functions are all used for surface integration schemes
float e_diff(const Vec k, const Vec q) {
    return epsilon(k+q) - epsilon(k);
}

float e_base_avg(const Vec k, const Vec q) {
    return epsilon(k);
}

float vp_diff(const Vec k, const Vec q) {
    Vec v;
    if (band_name == "simple_cubic_layered")
        v = fermi_velocity_SC_layered(k+q) - fermi_velocity_SC_layered(k);
    else if (band_name == "simple_cubic")
        v = fermi_velocity_SC(k+q) - fermi_velocity_SC(k);
    else if (band_name == "sphere")
        v = fermi_velocity_sphere(k+q) - fermi_velocity_sphere(k);
    else {
        cout << "No band structure specified\n";
        assert(1==2);
        return 0;
    }
    return v.norm();
}

// Fermi Velocity corresponds to energy band functions above
float vp(const Vec k) {
    if (band_name == "simple_cubic_layered")
        return fermi_velocity_SC_layered(k).norm();
    if (band_name == "simple_cubic") {
        return fermi_velocity_SC(k).norm();
    }
    if (band_name == "sphere")
        return fermi_velocity_sphere(k).norm();
    else {
        cout << "No band structure specified\n";
        exit(1);
        return 0;
    }
}

// Original function for V (NO CHANGE FROM ORIGINAL)
float V(const Vec k1, const Vec k2, float w, const float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube) {
    if (potential_name == "const")
        return potential_const(k1, k2);
    if (potential_name == "scalapino") {
        float t2 = potential_scalapino_cube(k1, k2, w, T, chi_cube);
        return t2;
    }
    if (potential_name == "scalapino_triplet")
        return potential_scalapino_triplet(k1, k2, T, w, chi_cube);
    if (potential_name == "phonon_coulomb") {
        Vec q = k1 - k2;
        return phonon_coulomb(q, c);
    }
    if (potential_name == "test")
        return potential_test(k1, k2);// / pow(2*M_PI, dim);
    else {
        cout << "Unknown Potential Function: " << potential_name << endl;
        exit(1);
        return 0;
    }
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
