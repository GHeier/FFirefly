#include <math.h>
#include <string>
#include "vec.h"
#include "cfg.h"
#include "band_structure.h"
#include "potential.h"

using namespace std;

// Global Variables
int n = 10; // Number of k points
int m = 40; // Number of chi points
int l = 5; // Number of frequency points
int dim = 3; // Number of dimensions)
string potential_name = "const";
string band_name = "simple_cubic";
               
// Constants
double t = 1.0;
double tn = 0.0;
double tnn = 0.0;
double U = 4.0;
double k_max = M_PI;
double mu = 1.0;
double w_D = 0.5;


void init_config(double &mu, double &U, double &t, double &tn, double &w_D, double new_mu, double new_U, double new_t, double new_tn, double new_w_D) {
    mu = new_mu;
    U = new_U;
    t = new_t;
    tn = new_tn;
    w_D = new_w_D;
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
        exit(1);
        return 0;
    }
}

double e_diff(const Vec k, const Vec q) {
    return epsilon(k+q) - epsilon(k);
}

double e_base_avg(const Vec k, const Vec q) {
    return epsilon(k);
}

double e_base(const Vec k, const Vec q) {
    return epsilon(k+q) + epsilon(k-q);
}

double e_split(const Vec k, const Vec q) {
    return epsilon(k+q) - epsilon(k-q);
}

double e_surface(const Vec k, const Vec q) {
    return epsilon(k) - mu;
} 

double vp_diff(const Vec k, const Vec q) {
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
    return v.vals.norm();
}

double vp_split(const Vec k, const Vec q) {
    Vec v;
    if (band_name == "simple_cubic_layered")
        v = fermi_velocity_SC_layered(k+q/2) - fermi_velocity_SC_layered(k-q/2);
    else if (band_name == "simple_cubic")
        v = fermi_velocity_SC(k+q/2) - fermi_velocity_SC(k-q/2);
    else if (band_name == "sphere")
        v = fermi_velocity_sphere(k+q/2) - fermi_velocity_sphere(k-q/2);
    else {
        cout << "No band structure specified\n";
        assert(1==2);
        return 0;
    }
    return v.vals.norm();
}

double vp_surface(const Vec k, const Vec q) {
    if (band_name == "simple_cubic_layered")
        return fermi_velocity_SC_layered(k).vals.norm();
    else if (band_name == "simple_cubic")
        return fermi_velocity_SC(k).vals.norm();
    else if (band_name == "sphere")
        return fermi_velocity_sphere(k).vals.norm();
    else {
        cout << "No band structure specified\n";
        assert(1==2);
        return 0;
    }
}

// Fermi Velocity corresponds to energy band functions above
double vp(const Vec k) {
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
// Potential functions
double V(const Vec k1, const Vec k2, double w, const double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube) {
    if (potential_name == "const") 
        return potential_const(k1, k2);
    if (potential_name == "scalapino") 
        //return potential_scal(k1, k2, T);
        return potential_scalapino_cube(k1, k2, w, T, chi_cube);
    if (potential_name == "scalapino_triplet") 
        return potential_scalapino_triplet(k1, k2, T, w, chi_cube);
    if (potential_name == "test") 
        return potential_test(k1, k2);// / pow(2*M_PI, dim);
    else {
        cout << "Unknown Potential Function: " << potential_name << endl;
        exit(1);
        return 0;
    }
}

// NOTE: Changes delta value as well
int get_num_points_from_delta(double &delta) {
    if (delta == 0) {
        delta = 0.0001;
    }
    int pts = 10*k_max/delta + 1;
    if (delta > 0.1) {
        delta = 0;
    }
    return pts;
}

// Gaussian integration constants
double weights_0th[1] = {2.0}; double * w0 = weights_0th;
double weights_1st[2] = {1.0, 1.0}; double * w1 = weights_1st;
double weights_2nd[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0}; double * w2 = weights_2nd;
double weights_3rd[4] = {0.347855, 0.652145, 0.652145, 0.347855}; double * w3 = weights_3rd;
double weights_4th[5] = {0.236927, 0.478629, 0.568889, 0.478629, 0.236927}; double * w4 = weights_4th;
double *weights[5] = {w0, w1, w2, w3, w4};

double points_0th[1] = {0}; double *p0 = points_0th;
double points_1st[2] = {-1/pow(3,0.5), 1/pow(3,0.5)}; double *p1 = points_1st;
double points_2nd[3] = {-pow(3/5,0.5), 0, pow(3/5,0.5)}; double *p2 = points_2nd;
double points_3rd[4] = {-0.861136, -0.339981, 0.339981, 0.861136}; double *p3 = points_3rd;
double points_4th[5] = {-0.90618, -0.538469, 0, 0.538469, 0.90618}; double *p4 = points_4th;

double *points[5] = {p0, p1, p2, p3, p4};

float sin_arr[31416]; 
float cos_arr[31416]; 
struct sin_arr_init {
    sin_arr_init() {
        for (int i = 0; i < 31416; i++) {
            sin_arr[i] = sin( (double)i / 10000.0 );
        }
    }
} sin_arr_init;
struct cos_arr_init {
    cos_arr_init() {
        for (int i = 0; i < 31416; i++) {
            cos_arr[i] = cos( (double)i / 10000.0 );
        }
    }
} cos_arr_init;

double get_sin(double x) {
    return sin_arr[ (int)(x*10000.0) % 31416 ];
}
double get_cos(double x) {
    return cos_arr[ (int)(x*10000.0) % 31416 ];
}
