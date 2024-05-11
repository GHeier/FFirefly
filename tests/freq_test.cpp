#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "cfg.h"
#include "potential.h"
#include "frequency_inclusion.hpp"

double analytic_fermi_gas_chi(Vec q, double w, double mu) {
    double q_plus_sqr = w + q.vals.squaredNorm();
    double q_minus_sqr = w - q.vals.squaredNorm();
    vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
    double N_3 = get_DOS(FS);
    return N_3 / mu * 0.75 * ( -1.0 + (4*q.vals.squaredNorm() - pow(q_minus_sqr,2)) / (8*pow(q.vals.norm(),3)) 
            * log( (1+q_minus_sqr/(2*q.vals.norm())) / ( 1 - q_minus_sqr/(2*q.vals.norm())))
            - (4*q.vals.squaredNorm() - pow(q_plus_sqr,2)) / (8*pow(q.vals.norm(),3)) 
            * log((1+q_plus_sqr/(2*q.vals.norm())) / (1 - q_plus_sqr / (2*q.vals.norm()))));
}

double numerical_fermi_gas_chi(Vec q, double w, double mu, unordered_map<double, vector<vector<vector<double>>>> &chi_map) {
    w = round(w, 6);
    return calculate_chi_from_cube(chi_map.at(w), q);
}

void plot_chis(double T, double w) {
    vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
    auto chi_map = chi_cube_freq(T, mu, get_DOS(FS));
    ofstream file("chi_freq_v_q.txt");
    //double mu = 1.0;
    double n = 50;
    for (int i = 1; i < n; i++) {
        double mag = M_PI * i / n;
        Vec q(mag, mag, mag);
        file << q << " " 
            << analytic_fermi_gas_chi(q, w, mu)
            << " " << numerical_fermi_gas_chi(q, w, mu, chi_map)
            << endl; 
    }
}

double integrate_chi(Vec q, double T, double w) {
    return integrate_susceptibility(q, T, mu, w, 100);
}

void plot_test_chi_analytic(double T, double w, double delta) {
    ofstream file("chi_test_analytic.txt");
    double n = 50;
    int pts = get_num_points_from_delta(delta);
    for (int i = 1; i < n; i++) {
        double mag = M_PI * i / (n-1);
        Vec q(mag, mag, mag);
        double integral = modified_integral_wrapper(q,T,mu,w,delta,pts);
        //double integral = integrate_susceptibility(q, T, mu, w, 100);
        file << q << " " << integral << endl;
    }
}
