#include <math.h>
#include <string>
#include "vec.h"
#include "cfg.h"
#include "band_structure.h"
#include "potential.h"

using namespace std;

// Global Variables
int n = 10; // Number of k points
int m = 41; // Number of chi points
int dim = 3; // Number of dimensions)
string potential_name = "scalapino";
string band_name = "simple_cubic";
               
// Constants
double t = 1.0;
double tn = 0.0;
double U = 2.0;
double k_max = M_PI;
double mu = 0.9;
double w_D = 1.0;

void init_config(double &mu, double &U, double &t, double &tn, double new_mu, double new_U, double new_t, double new_tn) {
    mu = new_mu;
    U = new_U;
    t = new_t;
    tn = new_tn;
}

// Energy band functions
double epsilon(const Vec k) {
    if (band_name == "simple_cubic_layered")
        return epsilon_SC_layered(k);
    if (band_name == "simple_cubic")
        return epsilon_SC(k, t, tn);
    if (band_name == "sphere")
        return epsilon_sphere(k);
    else {
        cout << "Unknown Band structure: " << band_name << endl;
        assert(1==2);
        return 0;
    }
}
// Fermi Velocity corresponds to energy band functions above
double vp(const Vec k) {
    if (band_name == "simple_cubic_layered")
        return fermi_velocity_SC_layered(k).vals.norm();
    if (band_name == "simple_cubic")
        return fermi_velocity_SC(k).vals.norm();
    if (band_name == "sphere")
        return fermi_velocity_sphere(k).vals.norm();
    else {
        cout << "No band structure specified\n";
        assert(1==2);
        return 0;
    }
}
// Potential functions
double V(const Vec k1, const Vec k2, const double T, const vector<vector<vector<double>>> &chi_cube) {
    if (potential_name == "const") 
        return potential_const(k1, k2);
    if (potential_name == "scalapino") 
        //return potential_scal(k1, k2, T);
        return potential_scalapino_cube(k1, k2, T, chi_cube);
    if (potential_name == "scalapino_triplet") 
        return potential_scalapino_triplet(k1, k2, T, chi_cube);
    if (potential_name == "test") 
        return potential_test(k1, k2) / pow(2*M_PI, dim);
    else {
        cout << "Unknown Potential Function: " << potential_name << endl;
        assert(1==2);
        return 0;
    }
}

