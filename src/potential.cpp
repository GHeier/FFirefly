/**
 * @file potential.cpp
 *
 * @brief List of all potentials available for use in code
 *
 * @author Griffin Heier
 */

#include <iostream>
#include <math.h>
#include <vector>

#include <omp.h>
#include <boost/functional/hash.hpp>
#include <unordered_map>

#include "potential.h"
#include "objects/vec.h"
#include "cfg.h"
#include "response/susceptibility.h"
#include "superconductor/utilities.h"

using namespace std;

// Potential functions
float V(const Vec k1, const Vec k2, float w, const float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube) {
    if (potential_name == "const") 
        return potential_const(k1, k2);
    if (potential_name == "scalapino") {
        float t2 = potential_scalapino_cube(k1, k2, w, T, chi_cube);
        return t2;
    }
    if (potential_name == "scalapino_triplet") 
        return potential_scalapino_triplet(k1, k2, T, w, chi_cube);
    if (potential_name == "test") 
        return potential_test(k1, k2);// / pow(2*M_PI, dim);
    else {
        cout << "Unknown Potential Function: " << potential_name << endl;
        exit(1);
    }
}


float potential_const(Vec k1, Vec k2) {
    return -1;
}

float potential_test(Vec k1, Vec k2) {
    Vec q1 = k1;
    Vec q2 = k2;
    return -1.0*( cos(q1(0)) - cos(q1(1)) )*( cos(q2(0)) - cos(q2(1)) ) + (-0.5)*sin(q1(0))*sin(q1(1))*sin(q2(0))*sin(q2(1));
}

float phonon_coulomb(Vec q) {
    float qx = q(0);
    float Vp = 1.0/3.0;
    if (q.norm() != 0) {
        Vp = 1/(1+2*qx*qx / pow(q.norm(), 2));
    }
    float Vc = 1 / (1 + q.norm());
    return Vp + Vc;
}

float potential_scal(Vec k1, Vec k2, float T) {
    Vec q_minus = to_IBZ_2(k1 - k2);
    Vec q_plus = to_IBZ_2(k1 + k2);
    
    float Xm = integrate_susceptibility(q_minus, T, mu, 0, 60);
    float Xp = integrate_susceptibility(q_plus, T, mu, 0, 60);

    float Vm = U*U * Xm / (1 - U*Xm) 
        + pow(U,3)*Xm*Xm / (1 - U*U * Xm*Xm);
    float Vp = U*U * Xp / (1 - U*Xp) 
        + pow(U,3)*Xp*Xp / (1 - U*U * Xp*Xp);
    return (Vm + Vp) / 2;
}

float potential_scalapino_cube(Vec k1, Vec k2, float w, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_map) {
    Vec q_minus = to_IBZ_2(k1 - k2);
    Vec q_plus = to_IBZ_2(k1 + k2);
    w = round_val(w,6);

    auto cube = chi_map.at(w);

    float Xm = calculate_chi_from_cube(cube, q_minus);
    float Xp = calculate_chi_from_cube(cube, q_plus);

    float Vm = U*U * Xm / (1 - U*Xm) 
        + pow(U,3)*Xm*Xm / (1 - U*U * Xm*Xm);
    float Vp = U*U * Xp / (1 - U*Xp) 
        + pow(U,3)*Xp*Xp / (1 - U*U * Xp*Xp);

    return 0.5 * (Vm + Vp);
}

float potential_scalapino_triplet(Vec k1, Vec k2, float w, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_map) {
    Vec q_minus = to_IBZ_2(k1 - k2);
    Vec q_plus = to_IBZ_2(k1 + k2);

    auto chi_cube = chi_map.at(w);

    float chi_minus = calculate_chi_from_cube(chi_cube, q_minus);
    float chi_plus = calculate_chi_from_cube(chi_cube, q_plus);

    float V_minus = -pow(U,2) * chi_minus / ( 1 - pow(U*chi_minus,2));
    float V_plus = -pow(U,2) * chi_plus / ( 1 - pow(U*chi_plus,2));

    return 0.5 * ( V_minus - V_plus);
}

