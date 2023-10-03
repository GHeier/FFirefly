#pragma once
#ifndef BAND_STRUCTURE_H
#define BAND_STRUCTURE_H

#include <string>
#include "vec.h"

using namespace std;


double epsilon_sphere(const Vec k);
Vec fermi_velocity_sphere(const Vec k);
double epsilon_SC(const Vec k, double t, double tn);
Vec fermi_velocity_SC(const Vec k);
double epsilon_SC_layered(const Vec k);
Vec fermi_velocity_SC_layered(const Vec k);

#endif
