#include <math.h>
#include <string>

#include "vec.h"
#include "cfg.h"

using namespace std;

double epsilon_sphere(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    if (dim == 2) q.vals[2] = 0;
    return pow(q.norm(),2);
}

Vec fermi_velocity_sphere(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    if (dim == 2) q.vals[2] = 0;
    return 2*q;
}

double epsilon_SC(const Vec k, double t, double tn) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    double val = 0.0;
    for (int i = 0; i < dim; i++) {
        val += -2*t*get_cos(q.vals[i]);
    }
    val += -4*tn*get_cos(q.vals[0])*get_cos(q.vals[1]);
    return val;
}

Vec fermi_velocity_SC(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    Vec v;
    for (int i = 0; i < dim; i++) {
        v.vals[i] = -get_sin(q.vals[i]);
    }
    v = -2*t*v;
    return v;
}

double epsilon_SC_layered(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    double val = 0.0;
    for (int i = 0; i < dim; i++) {
        if (i < 2) 
            val += (-2*t)*(get_cos(q.vals[i]));
        else
            val += (-2*tn)*(get_cos(q.vals[i]));
    }
    return val;
}

Vec fermi_velocity_SC_layered(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    Vec v;
    for (int i = 0; i < dim; i++) {
        if (i < 2) 
            v.vals[i] = (-2*t)*(-get_sin(q.vals[i]));
        else
            v.vals[i] = (-2*tn)*(-get_sin(q.vals[i]));
    }
    return v;
}


