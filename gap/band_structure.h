#pragma once
#ifndef BAND_STRUCTURE_H
#define BAND_STRUCTURE_H

#include <string>
#include "vec.h"

using namespace std;

float epsilon(const Vec k);
float e_diff(const Vec k, const Vec q);
float vp(const Vec k);
float vp_diff(const Vec k, const Vec q);

float epsilon_sphere(const Vec k);

Vec fermi_velocity_sphere(const Vec k);

float epsilon_SC(const Vec k, float t, float tn);

Vec fermi_velocity_SC(const Vec k);

float epsilon_SC_layered(const Vec k);

Vec fermi_velocity_SC_layered(const Vec k);

vector<Vec> get_FS(float E);

#endif
