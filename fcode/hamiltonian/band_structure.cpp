
#include <math.h>

#include "../objects/vec.hpp"
#include "../config/load/cpp_config.hpp"
#include "../config/load/c_config.h"
#include "../objects/surfaces.hpp"
#include "band_structure.hpp"

using namespace std;

// Energy band functions
float epsilon(int n, Vec k) {
    n--;
    if (n < 0) {
        printf("\nThe band index is negative/too small. Counting starts at 1\n");
        throw("The band index is negative/too small. Counting starts at 1\n");
    }
    if (band[n] == "simple_cubic_layered")
        return epsilon_SC_layered(n, k);
    if (band[n] == "tight_binding" && ibrav == 1)
        return epsilon_SC(n, k);
    if (band[n] == "fermi_gas") 
        return epsilon_fermi_gas(n, k);
    if (band[n] == "noband") {
        throw("The 0 band index is empty. Counting starts at 1\n");
    }
    else {
        cout << "Unknown Band structure: " << band[n] << endl;
        exit(1);
    }
}

// Difference functions are all used for surface integration schemes
float e_diff(int n, Vec k, Vec q) {
    return epsilon(n, k+q) - epsilon(n, k);
}

// Fermi Velocity corresponds to energy band functions above
float vp(int n, Vec k) {
    n--;
    if (band[n] == "simple_cubic_layered")
        return fermi_velocity_SC_layered(n, k).norm();
    if (band[n] == "simple_cubic") {
        return fermi_velocity_SC(n, k).norm();
    }
    if (band[n] == "fermi_gas")
        return fermi_velocity_fermi_gas(n, k).norm();
    else {
        cout << "Fermi velocity not available for band structure: " << band[n] << endl;
        exit(1);
    }
}

float vp_diff(int n, Vec k, Vec q) {
    Vec v;
    if (band[n] == "simple_cubic_layered")
        v = fermi_velocity_SC_layered(n, k+q) - fermi_velocity_SC_layered(n, k);
    else if (band[n] == "simple_cubic")
        v = fermi_velocity_SC(n, k+q) - fermi_velocity_SC(n, k);
    else if (band[n] == "fermi_gas")
        v = fermi_velocity_fermi_gas(n, k+q) - fermi_velocity_fermi_gas(n, k);
    else {
        cout << "No band structure specified\n";
        exit(1);
    }
    return v.norm();
}

/* ======================================================================
* ======================== Energy Band Functions ========================
*/

// Fermi gas
float epsilon_fermi_gas(int n, Vec k) {
    if (dimension < 3) k.z = 0;
    if (dimension < 4) k.w = 0;
    return pow(k.norm(), 2) / (2*eff_mass[n]);
}

Vec fermi_velocity_fermi_gas(int n, Vec k) {
    if (dimension < 3) k.z = 0;
    if (dimension < 4) k.w = 0;
    return 2*k / (2*eff_mass[n]);
}

// Cubic Lattice
float epsilon_SC(int n, Vec k) {
    float val = 0.0;
    for (int i = 0; i < dimension; i++) {
        val += -2*t0[n]*cos(k(i));
        val += -2*t2[n]*cos(k(i));
    }
    val += -4*t1[n]*cos(k(0))*cos(k(1));
    return val;
}

Vec fermi_velocity_SC(int n, Vec k) {
    Vec v;
    for (int i = 0; i < dimension; i++) {
        v(i) = -sin(k(i));
    }
    v = -2*t0[n]*v;
    return v;
}

// Cubic lattice with different hopping in z-direction
float epsilon_SC_layered(int n, Vec k) {
    float val = 0.0;
    for (int i = 0; i < dimension; i++) {
        if (i < 2) 
            val += (-2*t0[n])*(cos(k(i)));
        else
            val += (-2*t1[n])*(cos(k(i)));
    }
    return val;
}

Vec fermi_velocity_SC_layered(int n, Vec k) {
    Vec v;
    for (int i = 0; i < dimension; i++) {
        if (i < 2) 
            v(i) = (-2*t0[n])*(-sin(k(i)));
        else
            v(i) = (-2*t1[n])*(-sin(k(i)));
    }
    return v;
}

/* ======================================================================
* ======================== C Version of Band Structure ========================
*/

double epsilon_c(int n, double k[3]) {
    n--;
    Vec k_vec = Vec(k[0], k[1], k[2]);
    n++;
    return epsilon(n, k_vec);
}

double vp_c(int n, double k[3]) {
    Vec k_vec = Vec(k[0], k[1], k[2]);
    return vp(n, k_vec);
}


// E indicates the chemical potential at which to find the Fermi surface
vector<Vec> get_FS(float E) {
    printv("Finding Fermi Surface for E = %.2f\n", E);
    vector<Vec> FS;
    for (int i = 1; i <= nbnd; i++) {
        auto func = [i](Vec k) { return epsilon(i, k); };
        vector<Vec> temp = tetrahedron_method(func, E);
        for (auto k : temp) {
            k.n = i;
        }
        FS.insert(FS.end(), temp.begin(), temp.end());
    }
    return FS;
}

