#pragma once
#ifndef FERMI_SURFACE_H
#define FERMI_SURFACE_H
 
#include <vector>
#include "vec.h"

using namespace std;


struct VecAndEnergy {
    Vec vec;
    double energy;
};

bool operator<(const VecAndEnergy& left, const VecAndEnergy& right);
double triangle_area_from_points(Vec k1, Vec k2, Vec k3);
vector<VecAndEnergy> points_from_indices(double (*func)(Vec k, Vec q), Vec q, int i, int j, int k);
vector<Vec> points_in_tetrahedron(double (*func)(Vec k, Vec q), Vec q, double s_val, vector<VecAndEnergy> points);
bool surface_inside_cube(double s_val, vector<VecAndEnergy> p);
bool surface_inside_tetrahedron(double s_val, vector<VecAndEnergy> ep_points);
double area_in_corners(vector<Vec> cp);
vector<Vec> tetrahedron_method(double (*func)(Vec k, Vec q), Vec q, double s_val);
vector<Vec> tetrahedron_method_bounded(double (*func)(Vec k, Vec q), Vec q, double s_val);
vector<Vec> sort_by_adjacent(vector<Vec> &FS);
vector<VecAndEnergy> FS_to_q_shifted(vector<Vec> &FS, Vec q);
void filter_FS(vector<Vec> &FS);

#endif
