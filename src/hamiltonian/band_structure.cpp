
#include <fstream>
#include <math.h>

#include "../config/load/c_config.h"
#include "../config/load/cpp_config.hpp"
#include "../objects/surfaces.hpp"
#include "../objects/vec.hpp"
#include "band_structure.hpp"

using namespace std;

// Energy band functions
float epsilon(int n, Vec k) {
    n--;
    if (n < 0) {
        printf(
            "\nThe band index is negative/too small. Counting starts at 1\n");
        throw("The band index is negative/too small. Counting starts at 1\n");
    }
    if (band == "tight_binding" && celltype == "SC")
        return epsilon_SC(n, k);
    if (band == "tight_binding" && celltype == "FCC")
        return epsilon_FCC(n, k);
    if (band == "tight_binding" && celltype == "BCC")
        return epsilon_BCC(n, k);
    if (band == "fermi_gas")
        return epsilon_fermi_gas(n, k);
    if (band == "noband") {
        throw("The 0 band index is empty. Counting starts at 1\n");
    } else {
        cout << "Unknown Band structure: " << band << endl;
        exit(1);
    }
}

// Difference functions are all used for surface integration schemes
float e_diff(int n, Vec k, Vec q) { return epsilon(n, k + q) - epsilon(n, k); }

// Fermi Velocity corresponds to energy band functions above
float vp(int n, Vec k) {
    n--;
    if (band == "simple_cubic_layered")
        return fermi_velocity_SC_layered(n, k).norm();
    if (band == "tight_binding" && celltype == "SC") {
        return fermi_velocity_SC(n, k).norm();
    }
    if (band == "tight_binding" && celltype == "BCC") {
        return fermi_velocity_BCC(n, k).norm();
    }
    if (band == "tight_binding" && celltype == "FCC") {
        return fermi_velocity_FCC(n, k).norm();
    }
    if (band == "fermi_gas")
        return fermi_velocity_fermi_gas(n, k).norm();
    else {
        cout << "Fermi velocity not available for band structure: " << band
             << endl;
        exit(1);
    }
}

float vp_diff(int n, Vec k, Vec q) {
    Vec v;
    if (band == "simple_cubic_layered")
        v = fermi_velocity_SC_layered(n, k + q) -
            fermi_velocity_SC_layered(n, k);
    else if (band == "simple_cubic")
        v = fermi_velocity_SC(n, k + q) - fermi_velocity_SC(n, k);
    else if (band == "fermi_gas")
        v = fermi_velocity_fermi_gas(n, k + q) - fermi_velocity_fermi_gas(n, k);
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
    if (dimension < 3)
        k.z = 0;
    if (dimension < 4)
        k.w = 0;
    return pow(k.norm(), 2) / (2 * eff_mass);
}

Vec fermi_velocity_fermi_gas(int n, Vec k) {
    if (dimension < 3)
        k.z = 0;
    if (dimension < 4)
        k.w = 0;
    return 2 * k / (2 * eff_mass);
}

// Cubic Lattice
float epsilon_SC(int n, Vec k) {
    float val = 0.0;
    for (int i = 0; i < dimension; i++) {
        val += -2 * t0 * cos(k(i));
        val += -2 * t2 * cos(k(i));
    }
    val += -4 * t1 * cos(k(0)) * cos(k(1));
    return val;
}

Vec fermi_velocity_SC(int n, Vec k) {
    Vec v;
    for (int i = 0; i < dimension; i++) {
        v(i) = -sin(k(i));
    }
    v = -2 * t0 * v;
    return v;
}

// Cubic lattice with different hopping in z-direction
float epsilon_SC_layered(int n, Vec k) {
    float val = 0.0;
    for (int i = 0; i < dimension; i++) {
        if (i < 2)
            val += (-2 * t0) * (cos(k(i)));
        else
            val += (-2 * t1) * (cos(k(i)));
    }
    return val;
}

Vec fermi_velocity_SC_layered(int n, Vec k) {
    Vec v;
    for (int i = 0; i < dimension; i++) {
        if (i < 2)
            v(i) = (-2 * t0) * (-sin(k(i)));
        else
            v(i) = (-2 * t1) * (-sin(k(i)));
    }
    return v;
}

float epsilon_BCC(int n, Vec k) {
    float term = -8 * t0 * cos(k(0) / 2) * cos(k(1) / 2) * cos(k(2) / 2);
    for (int i = 0; i < dimension; i++) {
        term += -2 * t1 * cos(k(i));
    }
    return term;
}

Vec fermi_velocity_BCC(int n, Vec k) {
    float term = cos(k(0) / 2) * cos(k(1) / 2) * cos(k(2) / 2);
    Vec v;
    for (int i = 0; i < dimension; i++) {
        v(i) = -4 * t0 * sin(k(i) / 2) * term / cos(k(i) / 2);
        v(i) += -2 * t1 * sin(k(i));
    }
    return v;
}

float epsilon_FCC(int n, Vec k) {
    float val = 0.0;
    for (int i = 0; i < dimension; i++) {
        val += -4 * t0 * cos(k(i) / 2) * cos(k((i + 1) % dimension) / 2);
        val += -2 * t1 * cos(k(i));
    }
    return val;
}

Vec fermi_velocity_FCC(int n, Vec k) {
    Vec v;
    for (int i = 0; i < dimension; i++) {
        v(i) = 2 * t0 * (sin(k(i) / 2) * cos(k((i + 1) % dimension) / 2) +
                 sin(k((i - 1 + dimension) % dimension) / 2) *
                     cos(k(i) / 2));
        v(i) += -2 * t1 * sin(k(i));
    }
    return v;
}

/* ======================================================================
 * ======================== C Version of Band Structure ========================
 */

double epsilon_c(int n, double k[3]) {
    Vec k_vec = Vec(k[0], k[1], k[2]);
    return epsilon(n, k_vec);
}

double epsilon_c2d(int n, double k[2]) {
    Vec k_vec = Vec(k[0], k[1]);
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
    vector<Vec> temp;
    for (int i = 1; i <= nbnd; i++) {
        auto func = [i](Vec k) { return epsilon(i, k); };
        if (dimension == 3)
            temp = tetrahedron_method(func, E);
        else if (dimension == 2)
            temp = tetrahedron_method_2D(func, E);
        else {
            cout << "Dimension not supported\n";
            exit(1);
        }
        for (auto k : temp)
            k.n = i;
        FS.insert(FS.end(), temp.begin(), temp.end());
    }
    return FS;
}

float get_DOS(vector<Vec> &FS) {
    float sum = 0;
    for (auto k : FS) {
        sum += k.area / vp(k.n, k);
    }
    sum /= pow(2 * M_PI, dimension);
    printv("Density of States: %.5f\n", sum);
    return sum;
}
