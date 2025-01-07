#pragma once

#include <vector>
#include <string>
#include "../objects/vec.hpp"
#include "../response/susceptibility.hpp"

using namespace std;

struct indices_values {
    vector<int> indices;
    float value;
};

__attribute__((visibility("default")))
float V(const Vec k1, const Vec k2, string spin1="", string spin2="");

float potential_const();
float potential_test(Vec k1, Vec k2, string spin1, string spin2);

float phonon_coulomb(Vec q);

extern Susceptibility chi;
void load_chi(string filename, int dimension_);
void load_chi(vector<Vec> points, vector<complex<float>> values, int dimension_);
void get_chi(double q_c[3]);

float potential_FLEX(Vec k1, Vec k2, string spin1, string spin2);
float potential_FLEX_singlet(Vec k1, Vec k2);  
float potential_FLEX_triplet(Vec k1, Vec k2);

float Vs_c(double k1_c[3], double k2_c[3], const char* spin1_c, const char* spin2_c);
float V_c(double k1_c[3], double k2_c[3]);

