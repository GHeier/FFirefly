
#include <math.h>

#include "../objects/vec.hpp"
#include "../config/load/cpp_config.hpp"
#include "../objects/surfaces.hpp"
#include "band_structure.hpp"

using namespace std;

// Energy band functions
float epsilon(int n, Vec k) {
    if (band[n] == "simple_cubic_layered")
        return epsilon_SC_layered(k);
    if (band[n] == "simple_cubic")
        return epsilon_SC(k, t, tn);
    if (band[n] == "sphere")
        return epsilon_sphere(k);
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
    if (band[n] == "simple_cubic_layered")
        return fermi_velocity_SC_layered(k).norm();
    if (band[n] == "simple_cubic") {
        return fermi_velocity_SC(k).norm();
    }
    if (band[n] == "sphere")
        return fermi_velocity_sphere(k).norm();
    else {
        cout << "Fermi velocity not available for band structure: " << band[n] << endl;
        exit(1);
    }
}

float vp_diff(int n, Vec k, Vec q) {
    Vec v;
    if (band[n] == "simple_cubic_layered")
        v = fermi_velocity_SC_layered(k+q) - fermi_velocity_SC_layered(k);
    else if (band[n] == "simple_cubic")
        v = fermi_velocity_SC(k+q) - fermi_velocity_SC(k);
    else if (band[n] == "sphere")
        v = fermi_velocity_sphere(k+q) - fermi_velocity_sphere(k);
    else {
        cout << "No band structure specified\n";
        exit(1);
    }
    return v.norm();
}

// Fermi gas
float epsilon_sphere(Vec k) {
    if (dim < 3) k.z = 0;
    if (dim < 4) k.w = 0;
    return pow(k.norm(), 2);
}

Vec fermi_velocity_sphere(Vec k) {
    if (dim < 3) k.z = 0;
    if (dim < 4) k.w = 0;
    return 2*k;
}

// Cubic Lattice
float epsilon_SC(Vec k, float t, float tn) {
    float val = 0.0;
    for (int i = 0; i < dim; i++) {
        val += -2*t*cos(k(i));
        val += -2*tnn*cos(k(i));
    }
    val += -4*tn*cos(k(0))*cos(k(1));
    return val;
}

Vec fermi_velocity_SC(Vec k) {
    Vec v;
    for (int i = 0; i < dim; i++) {
        v(i) = -sin(k(i));
    }
    v = -2*t*v;
    return v;
}

// Cubic lattice with different hopping in z-direction
float epsilon_SC_layered(Vec k) {
    float val = 0.0;
    for (int i = 0; i < dim; i++) {
        if (i < 2) 
            val += (-2*t)*(cos(k(i)));
        else
            val += (-2*tn)*(cos(k(i)));
    }
    return val;
}

Vec fermi_velocity_SC_layered(Vec k) {
    Vec v;
    for (int i = 0; i < dim; i++) {
        if (i < 2) 
            v(i) = (-2*t)*(-sin(k(i)));
        else
            v(i) = (-2*tn)*(-sin(k(i)));
    }
    return v;
}

// E indicates the chemical potential at which to find the Fermi surface
vector<Vec> get_FS(float E) {
    printf("Finding Fermi Surface for E = %.2f\n", E);
    vector<Vec> FS;
    for (int i = 0; i < nbnd; i++) {
        auto func = [i](Vec k) { return epsilon(i, k); };
        vector<Vec> temp = tetrahedron_method(func, E);
        for (auto k : temp) {
            k.n = i;
        }
        FS.insert(FS.end(), temp.begin(), temp.end());
    }
    return FS;
}

