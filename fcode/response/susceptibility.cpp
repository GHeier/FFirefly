/**
 * Calculates the susceptibility of the system. Slower than response and sparse_ir algorithms
 *
 * Author: Griffin Heier
 */

#include <iostream>
#include <math.h>
#include <complex>
#include <string>
#include <algorithm>
#include <functional>

#include <omp.h>
#include <boost/functional/hash.hpp>

#include "../objects/vec.hpp"
#include "../objects/fastfield.hpp"
#include "../config/load/cpp_config.hpp"
#include "susceptibility.hpp"
#include "../hamiltonian/band_structure.hpp"

using namespace std;

Susceptibility::Susceptibility() {}

Susceptibility::Susceptibility(vector<Vec> points, vector<complex<float>> values, int dimension, bool is_complex) {
    chi = FastScalarField(points, values, dimension+1, is_complex);
}

Susceptibility::Susceptibility(string filename, int dimension, bool is_complex) {
    chi = FastScalarField(filename, dimension+1, is_complex);
}

complex<float> Susceptibility::operator() (Vec point, float w) {
    point = to_IBZ(point);
    point.dimension++;
    point.w = fabs(w);
    return chi(point);
}

float Susceptibility::operator() (Vec point) {
    cout << "point: " << point << endl;
    point = to_IBZ(point);
    cout << "point: " << point << endl;
    return chi(point);
}

Susceptibility Susceptibility::operator= (Susceptibility other) {
    chi = other.chi;
    return *this;
}

//Fermi-Dirac distribution function
float fermi_dirac(float E, float T) {
    if (T == 0) {
        if (E < 0) return 1;
        return 0;
    }
    return 1 / (1 + exp(E/T));
}

float ratio(Vec k, Vec q, float w, float T) {
    float e_k = epsilon(k.n, k) - fermi_energy;
    float e_qk = epsilon(k.n, k+q) - fermi_energy;
    float dE = e_qk - e_k;
    float f_k = fermi_dirac(e_k, T);
    float f_qk = fermi_dirac(e_qk, T);

    if (fabs(dE) < 0.0001 and fabs(w) < 0.0001) {
        if (T == 0 or exp(e_k/T) > 1e6)
            return e_k < 0;
        float temp = 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        return temp;
    }
    if (fabs(w - dE) == 0) return 0;
    return (f_qk - f_k) / (w - dE);
}

Vec to_IBZ(const Vec k) {
    Vec q = k;
    float x = q(0), y = q(1), z = q(2);
    x = abs(x); y = abs(y); z = abs(z);
    if (x > M_PI) x = - (x - 2*M_PI);
    if (y > M_PI) y = - (y - 2*M_PI);
    if (z > M_PI) z = - (z - 2*M_PI);
    if (dimension == 3) {
        float arr[] = {x, y, z};
        sort(arr, arr+3, greater<float>());
        auto& [a, b, c] = arr;
        Vec result(a, b, c);
        return result;
    }
    else if (dimension == 2) {
        float arr[] = {x, y};
        sort(arr, arr+2, greater<float>());
        auto& [a, b] = arr;
        Vec result(a, b, z);
        return result;
    }
    else {
        cout << "Wrong Dimension\n";
        return q;
    }
}

