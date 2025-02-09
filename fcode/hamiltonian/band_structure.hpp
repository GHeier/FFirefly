#pragma once

#include <string>
#include "../objects/vec.hpp"

using namespace std;

float epsilon(int n, Vec k);
float e_diff(int n, Vec k, Vec q);
float vp(int n, Vec k);
float vp_diff(int n, Vec k, Vec q);
float epsilon_fermi_gas(int n, Vec k);
Vec fermi_velocity_fermi_gas(int n, Vec k);
float epsilon_SC(int n, Vec k);
Vec fermi_velocity_SC(int n, Vec k);
float epsilon_SC_layered(int n, Vec k);
Vec fermi_velocity_SC_layered(int n, Vec k);

extern "C" {
    double epsilon_c(int n, double k[3]);
    double epsilon_c2d(int n, double k[2]);
    double vp_c(int n, double k[3]);
}

vector<Vec> get_FS(float E);

