#pragma once

#include <string>
#include "../objects/vec.hpp"

using namespace std;

float epsilon(int n, Vec k);
float e_diff(int n, Vec k, Vec q);
float vp(int n, Vec k);
float vp_diff(int n, Vec k, Vec q);
float epsilon_sphere(Vec k);
Vec fermi_velocity_sphere(Vec k);
float epsilon_SC(Vec k, float t, float tn);
Vec fermi_velocity_SC(Vec k);
float epsilon_SC_layered(Vec k);
Vec fermi_velocity_SC_layered(Vec k);
vector<Vec> get_FS(float E);

