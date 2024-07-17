/**
 * @file band_structure.cpp
 *
 * @brief Band structure functions for sipmle systems
 *
 * @author Griffin Heier
 */


#include <math.h>
#include <string>

#include "vec.h"
#include "cfg.h"

using namespace std;

// Fermi gas
float epsilon_sphere(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    if (dim == 2) q.vals[2] = 0;
<<<<<<< HEAD
    return pow(q.norm(), 2);
=======
    return pow(q.norm(),2);
>>>>>>> origin/main
}

Vec fermi_velocity_sphere(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    if (dim == 2) q.vals[2] = 0;
    return 2*q;
}

// Cubic Lattice
float epsilon_SC(const Vec k, float t, float tn) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    float val = 0.0;
    for (int i = 0; i < dim; i++) {
        val += -2*t*cos(q.vals[i]);
<<<<<<< HEAD
        val += -2*tnn*cos(q.vals[i]);
=======
>>>>>>> origin/main
    }
    val += -4*tn*cos(q.vals[0])*cos(q.vals[1]);
    return val;
}

Vec fermi_velocity_SC(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    Vec v;
    for (int i = 0; i < dim; i++) {
        v.vals[i] = -sin(q.vals[i]);
    }
    v = -2*t*v;
    return v;
}

// Cubic lattice with different hopping in z-direction
float epsilon_SC_layered(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    float val = 0.0;
    for (int i = 0; i < dim; i++) {
        if (i < 2) 
            val += (-2*t)*(cos(q.vals[i]));
        else
            val += (-2*tn)*(cos(q.vals[i]));
    }
    return val;
}

Vec fermi_velocity_SC_layered(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    Vec v;
    for (int i = 0; i < dim; i++) {
        if (i < 2) 
            v.vals[i] = (-2*t)*(-sin(q.vals[i]));
        else
            v.vals[i] = (-2*tn)*(-sin(q.vals[i]));
    }
    return v;
}


