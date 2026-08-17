
#include <fstream>
#include <math.h>

#include "src/config/load/c_config.h"
#include "src/config/load/cpp_config.hpp"
#include "src/objects/surfaces.hpp"
#include "src/objects/vec.hpp"
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
    if (band == "tight_binding" && celltype == "LSCO")
        return epsilon_LSCO(n, k);
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
    if (band == "tight_binding" && celltype == "SC") {
        return fermi_velocity_SC(n, k).norm();
    }
    if (band == "tight_binding" && celltype == "BCC") {
        return fermi_velocity_BCC(n, k).norm();
    }
    if (band == "tight_binding" && celltype == "FCC") {
        return fermi_velocity_FCC(n, k).norm();
    }
    if (band == "tight_binding" && celltype == "LSCO") {
        return fermi_velocity_LSCO(n, k).norm();
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
    if (band == "LSCO")
        v = fermi_velocity_LSCO(n, k + q) -
            fermi_velocity_LSCO(n, k);
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
        val += -2 * t0 * cos(k(i));        // NN hopping
        val += -2 * t2 * cos(2 * k(i));    // NNNN hopping (fixed: was cos(k(i)))
    }
    val += -4 * t1 * cos(k(0)) * cos(k(1));  // NNN hopping (diagonal)
    return val;
}

Vec fermi_velocity_SC(int n, Vec k) {
    Vec v;
    for (int i = 0; i < dimension; i++) {
        v(i) += 2 * t0 * sin(k(i));
        v(i) += 4 * t2 * sin(2 * k(i));    // NNNN hopping (fixed: was cos(k(i)))
    }
    v(0) += 4 * t1 * sin(k(0)) * cos(k(1));
    v(1) += 4 * t1 * cos(k(0)) * sin(k(1));
    return v;
}

// Cubic lattice with different hopping in z-direction
float epsilon_LSCO(int n, Vec k) {
    float a = cell[0][0];
    float c = cell[2][2];
    float val = 0.0;
// 2D dispersion
    val += -2 * t0 * (cos(k(0)) + cos(k(1)));
    val += -4 * t1 * cos(k(0)) * cos(k(1));
    val += -2 * t2 * (cos(2*k(0)) + cos(2*k(1)));
    val += -2 * t3 * (cos(2*k(0))*cos(k(1)) + cos(k(0))*cos(2*k(1)));
// kz dispersion 
    // Option 1 (suboptimal fit)
    //val += -2 * tz0 * pow((cos(k(0)) - cos(k(1))),2) * cos(k(0)/2) * cos(k(1)/2) * cos(k(2)/2 * c/a);
    // Option 2 (empirically better fit from https://arxiv.org/pdf/cond-mat/0503064)
    float Sxy = cos(k(0)/2) * cos(k(1)/2);
    float a0 = 0.083;
    val += -2 * (tz0 * cos(k(2)/2) + tz1 * pow(cos(k(2)/2), 2)) * ((cos(k(0)) - cos(k(1))) * 2 + a0 * Sxy*Sxy) * Sxy;
    return val;
}

Vec fermi_velocity_LSCO(int n, Vec k) {
    double kx = k(0);
    double ky = k(1);
    double kz = k(2);

    Vec v;

    // 2D part (common to both kz options)
    v(0) = 2*t0*sin(kx) + 4*t1*sin(kx)*cos(ky) + 4*t2*sin(2*kx)
         + 4*t3*sin(2*kx)*cos(ky) + 2*t3*sin(kx)*cos(2*ky);
    v(1) = 2*t0*sin(ky) + 4*t1*cos(kx)*sin(ky) + 4*t2*sin(2*ky)
         + 2*t3*cos(2*kx)*sin(ky) + 4*t3*cos(kx)*sin(2*ky);
    v(2) = 0;

    // Option 1 (suboptimal fit)
    // double a = cell[0][0];
    // double c = cell[2][2];
    // double A = cos(kx) - cos(ky);
    // double B = cos(kx/2);
    // double C = cos(ky/2);
    // double D = cos(kz*c/(2*a));
    // v(0) += 2*tz0*D*C * (2*A*sin(kx)*B + 0.5*A*A*sin(kx/2));
    // v(1) -= 2*tz0*D*B * (2*A*sin(ky)*C - 0.5*A*A*sin(ky/2));
    // v(2) += tz0*(c/a) * A*A*B*C * sin(kz*c/(2*a));

    // Option 2 (empirically better fit from https://arxiv.org/pdf/cond-mat/0503064)
    double a0    = 0.083;
    double Sxy   = cos(kx/2) * cos(ky/2);
    double A2    = cos(kx) - cos(ky);
    double g     = 2*A2 + a0*Sxy*Sxy;
    double F     = tz0 * cos(kz/2) + tz1 * cos(kz/2)*cos(kz/2);
    double coeff = A2 + 1.5*a0*Sxy*Sxy;  // = a0*Sxy^2 + g/2
    v(0) += 2*F * (2*sin(kx)*Sxy + coeff * sin(kx/2)*cos(ky/2));
    v(1) += 2*F * (-2*sin(ky)*Sxy + coeff * sin(ky/2)*cos(kx/2));
    v(2) += (tz0 + 2*tz1*cos(kz/2)) * sin(kz/2) * g * Sxy;

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
