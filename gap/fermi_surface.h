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
vector<VecAndEnergy> points_from_indices(double (*func)(Vec k, Vec q), Vec q, int i, int j, int k, int divs);
vector<VecAndEnergy> points_from_points(double (*func)(Vec k, Vec q), Vec q, double x1, double x2, double y1, double y2, double z1, double z2);
vector<Vec> points_in_tetrahedron(double (*func)(Vec k, Vec q), Vec q, double s_val, vector<VecAndEnergy> points);
bool surface_inside_cube(double s_val, vector<VecAndEnergy> p);
bool surface_inside_tetrahedron(double s_val, vector<VecAndEnergy> ep_points);
double area_in_corners(vector<Vec> cp);
vector<Vec> tetrahedron_method(double (*func)(Vec k, Vec q), Vec q, double s_val);
double tetrahedron_sum(double (*func)(Vec k, Vec q), double (*func_diff)(Vec k, Vec q), Vec q, double s_val, double w, double T);
pair<int, int> get_index_and_length(double L, double U, vector<double> &sortedList);
double tetrahedron_sum_continuous(double (*func)(Vec k, Vec q), double (*func_diff)(Vec k, Vec q), Vec q, vector<double> &svals, double w, double T);
double analytic_tetrahedron_sum(Vec q, double w);

#endif
