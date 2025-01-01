#include <math.h>
#include <fstream>
#include "../gap/config/load/cpp_config.h"
#include "../gap/frequency_inclusion.hpp"
#include "../gap/utilities.h"
#include "../gap/frequency_inclusion.hpp"
#include "../gap/susceptibility.h"
#include "../gap/band_structure.h"
#include "../gap/vec.h"
#include "../gap/solver.h"

using namespace std;

float analytic_fermi_gas_chi(Vec q, float w, float mu) {
    float q_plus_sqr = w + pow(q.norm(), 2);
    float q_minus_sqr = w - pow(q.norm(), 2);
    vector<Vec> FS = get_FS(mu);
    float N_3 = get_DOS(FS);
    return N_3 / mu * 0.75 * ( -1.0 + (4*pow(q.norm(),2) - pow(q_minus_sqr,2)) / (8*pow(q.norm(),3)) 
            * log( (1+q_minus_sqr/(2*q.norm())) / ( 1 - q_minus_sqr/(2*q.norm())))
            - (4*pow(q.norm(),2) - pow(q_plus_sqr,2)) / (8*pow(q.norm(),3)) 
            * log((1+q_plus_sqr/(2*q.norm())) / (1 - q_plus_sqr / (2*q.norm()))));
}

float numerical_fermi_gas_chi(Vec q, float w, float mu, unordered_map<float, vector<vector<vector<float>>>> &chi_map) {
    w = round_val(w, 6);
    return calculate_chi_from_cube_map(chi_map, q, w);
}

void plot_chis(float T, float w) {
    vector<Vec> FS = get_FS(mu);
    auto chi_map = chi_cube_freq(T, mu);
    ofstream file("chi_freq_v_q.txt");
    //float mu = 1.0;
    float n = 50;
    for (int i = 1; i < n; i++) {
        float mag = M_PI * i / n;
        Vec q(mag, mag, mag);
        file << q << " " 
            << analytic_fermi_gas_chi(q, w, mu)
            << " " << numerical_fermi_gas_chi(q, w, mu, chi_map)
            << endl; 
    }
}

float integrate_chi(Vec q, float T, float w) {
    return integrate_susceptibility(q, T, mu, w, 100);
}

void plot_test_chi_analytic(float T, float w, float delta) {
    ofstream file("chi_test_analytic.txt");
    float n = 50;
    for (int i = 1; i < n; i++) {
        float mag = M_PI * i / (n-1);
        Vec q(mag, mag, mag);
        float integral = integrate_susceptibility(q, T, mu, w, 100);
        file << q << " " << integral << endl;
    }
}
