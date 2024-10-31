/**
 * @file potential.cpp
 *
 * @brief List of all potentials available for use in code
 *
 * @author Griffin Heier
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>

#include <omp.h>
#include <boost/functional/hash.hpp>
#include <tuple>
#include <unordered_map>
#include <algorithm>

#include "utilities.h"
#include "potential.h"
#include "vec.h"
#include "fermi_surface.h"
#include "band_structure.h"
#include "cfg.h"
#include "frequency_inclusion.hpp"
#include "susceptibility.h"

using namespace std;


float potential_const(Vec k1, Vec k2) {
    return -1;
}

float potential_test(Vec k1, Vec k2) {
    Vec q1 = k1;
    Vec q2 = k2;
    if (q1.cartesian == false) q1.to_cartesian();
    if (q2.cartesian == false) q2.to_cartesian();
    return -1*( cos(q1.vals[0]) - cos(q1.vals[1]) )*( cos(q2.vals[0]) - cos(q2.vals[1]) ) + (-0.5)*sin(q1.vals[0])*sin(q1.vals[1])*sin(q2.vals[0])*sin(q2.vals[1]);
}

float phonon_coulomb(Vec q, Vec c) {
    if (q.cartesian == false) q.to_cartesian();
    if (c.cartesian == false) c.to_cartesian();

    float cx = c.vals[0];
    float cy = c.vals[1];
    float qx = q.vals[0];
    float qy = q.vals[1];

    float α = ((cy * cy) - (cx * cx)) / (cx * cx);
    float q_s_sq = (8 * M_PI * e * e) / (Vol * ε);
    float C_sq = C * C;
    float denom_Vp = (N * M * (cx*cx) * (1 + ((α*qx*qx)/(q.norm()*q.norm()))));

    float Vp = 0;
    float Vc = 0;

    // 2D norm:
    float norm = sqrt((qx*qx)+(qy*qy));
    //if (norm != q.norm()){
	    //cout << "DEBUG: q.norm() = " << q.norm() << ", but manually computed norm is " << norm << endl;
    //}

    cout << "DEBUG: q.norm() = " << q.norm() << endl;
    cout << "DEBUG: q_s_sq = " << q_s_sq << endl;

    if (q.norm() == 0) {
        Vp = 0.0;
        cout << "DEBUG: norm of q is zero, setting Vp to " << Vp << endl;

        if (q_s_sq == 0) {
            Vc = 0;
            cout << "DEBUG: norm of q and q_s_sq are zero, setting Vc to " << Vc << endl;
        }
	    }
    cout << "DEBUG: numerator of Vp = " << (C_sq) << endl;
    cout << "DEBUG: denominator of Vp = " << (denom_Vp) << endl;
    cout << "DEBUG: N = " << N << endl;
    if (q.norm() != 0 && q_s_sq != 0) {
    	Vp = ((C_sq) / (denom_Vp));
    	Vc = (4 * M_PI * e * e) / (Vol * ε * ((q.norm() * q.norm()) + q_s_sq));
    	cout << "DEBUG: Calculated Vp = " << Vp << endl;
    	cout << "DEBUG: Calculated Vc = " << Vc << endl;
    }

    cout << "DEBUG: Returning Vp + Vc = " << (Vp + Vc) << endl;
    cout << endl;
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
    w = round(w,6);

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

