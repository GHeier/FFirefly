#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <boost/math/tools/roots.hpp>
#include <algorithm>
#include "cfg.h"
#include "vec.h"
#include "frequency_inclusion.hpp"
#include "fermi_surface.h"

using namespace std;

bool operator<(const VecAndEnergy& left, const VecAndEnergy& right) {
    return left.energy < right.energy;
}

double triangle_area_from_points(Vec k1, Vec k2, Vec k3) {
    auto triangle_area = [](double d1, double d2, double d3) { 
        double s = (d1 + d2 + d3)/2; 
        return pow(s*(s-d1)*(s-d2)*(s-d3), 0.5);
    };
    // Define distances 
    Vec k12 = k1 - k2; if (k12.cartesian) k12.to_spherical();
    Vec k23 = k2 - k3; if (k23.cartesian) k23.to_spherical();
    Vec k13 = k1 - k3; if (k13.cartesian) k13.to_spherical();
    double d12 = k12.vals(0);
    double d23 = k23.vals(0);
    double d13 = k13.vals(0);

    // Calculate areas if the triangle is on the fermi surface
    double A = 0;
    A = triangle_area(d12, d23, d13);
    if (isnan(A)) {
        return 0;
    }

    return A;
}

vector<VecAndEnergy> points_from_indices(double (*func)(Vec k, Vec q), Vec q, int i, int j, int k) {
    double x1 = 2*k_max * i / n       - k_max; 
    double x2 = 2*k_max * (i+1) / n   - k_max; 
    double y1 = 2*k_max * j / n       - k_max; 
    double y2 = 2*k_max * (j+1) / n   - k_max; 
    double z1 = 2*k_max * k / n       - k_max; 
    double z2 = 2*k_max * (k+1) / n   - k_max; 

    Vec p1(x1, y1, z1);
    Vec p2(x2, y1, z1);
    Vec p3(x2, y2, z1);
    Vec p4(x1, y2, z1);
    Vec p5(x1, y1, z2);
    Vec p6(x2, y1, z2);
    Vec p7(x2, y2, z2);
    Vec p8(x1, y2, z2);

    vector<VecAndEnergy> points(8); 
    VecAndEnergy point1 = {p1, func(p1, q)}; points[0] = point1;
    VecAndEnergy point2 = {p2, func(p2, q)}; points[1] = point2;
    VecAndEnergy point3 = {p3, func(p3, q)}; points[2] = point3;
    VecAndEnergy point4 = {p4, func(p4, q)}; points[3] = point4;
    VecAndEnergy point5 = {p5, func(p5, q)}; points[4] = point5;
    VecAndEnergy point6 = {p6, func(p6, q)}; points[5] = point6;
    VecAndEnergy point7 = {p7, func(p7, q)}; points[6] = point7;
    VecAndEnergy point8 = {p8, func(p8, q)}; points[7] = point8;

    return points;
}

vector<VecAndEnergy> points_from_points(double (*func)(Vec k, Vec q), Vec q, double x1, double x2, double y1, double y2, double z1, double z2) {
    Vec p1(x1, y1, z1);
    Vec p2(x2, y1, z1);
    Vec p3(x2, y2, z1);
    Vec p4(x1, y2, z1);
    Vec p5(x1, y1, z2);
    Vec p6(x2, y1, z2);
    Vec p7(x2, y2, z2);
    Vec p8(x1, y2, z2);

    vector<VecAndEnergy> points(8); 
    VecAndEnergy point1 = {p1, func(p1, q)}; points[0] = point1;
    VecAndEnergy point2 = {p2, func(p2, q)}; points[1] = point2;
    VecAndEnergy point3 = {p3, func(p3, q)}; points[2] = point3;
    VecAndEnergy point4 = {p4, func(p4, q)}; points[3] = point4;
    VecAndEnergy point5 = {p5, func(p5, q)}; points[4] = point5;
    VecAndEnergy point6 = {p6, func(p6, q)}; points[5] = point6;
    VecAndEnergy point7 = {p7, func(p7, q)}; points[6] = point7;
    VecAndEnergy point8 = {p8, func(p8, q)}; points[7] = point8;

    return points;
}

vector<Vec> points_in_tetrahedron(double (*func)(Vec k, Vec q), Vec q, double s_val, vector<VecAndEnergy> points) {
    sort(points.begin(), points.end());
    Vec k1 = points[0].vec, k2 = points[1].vec, k3 = points[2].vec, k4 = points[3].vec;
    double ep1 = func(k1, q), ep2 = func(k2, q), ep3 = func(k3, q), ep4 = func(k4, q);

    Vec empty;

    Vec k12 = (k2-k1) * (s_val - ep1) / (ep2 - ep1) + k1;
    Vec k13 = (k3-k1) * (s_val - ep1) / (ep3 - ep1) + k1;
    Vec k14 = (k4-k1) * (s_val - ep1) / (ep4 - ep1) + k1;
    Vec k24 = (k4-k2) * (s_val - ep2) / (ep4 - ep2) + k2;
    Vec k34 = (k4-k3) * (s_val - ep3) / (ep4 - ep3) + k3;
    Vec k23 = (k3-k2) * (s_val - ep2) / (ep3 - ep2) + k2;
    
    vector<Vec> return_points(4, empty);


    if ( s_val > ep1 and s_val <= ep2) {
        return_points[0] = k12;
        return_points[1] = k13;
        return_points[2] = k14;
    }
    if ( s_val > ep3 and s_val <= ep4) {
        return_points[0] = k14;
        return_points[1] = k24;
        return_points[2] = k34;
    }
    if ( s_val > ep2 and s_val <= ep3) {
        return_points[0] = k24;
        return_points[1] = k23;
        return_points[2] = k13;
        return_points[3] = k14;
    }

    // Run condition for when plane aligns with one of the planes of the tetrahedron
    int times_not_equal = 0;
    Vec not_equal;
    for (int i = 0; i < 4; i++)
        if (s_val != func(points[i].vec, q)) {
            times_not_equal++;
            not_equal = points[i].vec;
        }
    if (times_not_equal == 1) {
        int iter = -1;
        for (int i = 0; i < 4; i++) {
            if (points[i].vec == not_equal) continue;
            return_points[iter] = points[i].vec;
            iter++;
        }
        return_points[3] = empty;
    }

    return return_points;
}

bool surface_inside_cube(double s_val, vector<VecAndEnergy> p) {
    sort(p.begin(), p.end());
    return (p[7].energy - s_val) / (p[0].energy - s_val) < 0;
}

bool surface_inside_tetrahedron(double s_val, vector<VecAndEnergy> ep_points) {
    sort(ep_points.begin(), ep_points.end());
    return ((ep_points[3].energy)-s_val) / ((ep_points[0].energy) - s_val) < 0;
}

double area_in_corners(vector<Vec> cp) {
    Vec empty;
    Vec k1 = cp[0]; Vec k2 = cp[1]; Vec k3 = cp[2]; Vec k4 = cp[3];
    if (k4 == empty) return triangle_area_from_points(k1, k2, k3);

    double A1 = 0, A2 = 0;
    A1 = triangle_area_from_points(k1, k2, k4);
    A2 = triangle_area_from_points(k3, k2, k4);



    return A1 + A2;
}

vector<Vec> tetrahedron_method(double (*func)(Vec k, Vec q), Vec q, double s_val) {
    vector<vector<double>> tetrahedrons {
        {1, 2, 3, 5}, 
        {1, 3, 4, 5},
        {2, 5, 6, 3},
        {4, 5, 8, 3},
        {5, 8, 7, 3},
        {5, 6, 7, 3}
    };

    vector<Vec> FS;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n * (dim%2) + 1 * ((dim+1)%2); k++) {
                vector<VecAndEnergy> points = points_from_indices(func, q, i, j, k);
                if (not surface_inside_cube(s_val, points)) continue;

                bool all_big = true;
                for (int c = 0; c < 6; c++) {

                    vector<VecAndEnergy> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }

                    if (not surface_inside_tetrahedron(s_val, ep_points)) continue;
                    vector<Vec> corner_points = points_in_tetrahedron(func, q, s_val, ep_points);

                    Vec average;

                    double b = 0;
                    if (corner_points[3] == average) b = 1.0;

                    for (Vec q : corner_points) {
                        average = (q + average);
                    }
                    average = average / (4-b);

                    double A = area_in_corners(corner_points);
                    if (dim == 2) A *= n / (2*k_max);
                    Vec k_point = average; k_point.area = A;
                    k_point.freq = s_val;
                    FS.push_back(k_point);
                }
            }
        }
    }
    //printf("Iters: %d\n", iters);
    return FS;
}

double tetrahedron_sum(double (*func)(Vec k, Vec q), double (*func_diff)(Vec k, Vec q), Vec q, double s_val, double w, double T) {
    vector<vector<double>> tetrahedrons {
        {1, 2, 3, 5}, 
        {1, 3, 4, 5},
        {2, 5, 6, 3},
        {4, 5, 8, 3},
        {5, 8, 7, 3},
        {5, 6, 7, 3}
    };

    double sum = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n * (dim%2) + 1 * ((dim+1)%2); k++) {
                vector<VecAndEnergy> points = points_from_indices(func, q, i, j, k);
                if (not surface_inside_cube(s_val, points)) continue;

                for (int c = 0; c < 6; c++) {

                    vector<VecAndEnergy> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }

                    if (not surface_inside_tetrahedron(s_val, ep_points)) continue;
                    vector<Vec> corner_points = points_in_tetrahedron(func, q, s_val, ep_points);

                    Vec average;

                    double b = 0;
                    if (corner_points[3] == average) b = 1.0;

                    for (Vec q : corner_points) {
                        average = (q + average);
                    }
                    average = average / (4-b);

                    double A = area_in_corners(corner_points);
                    if (dim == 2) A *= n / (2*k_max);
                    Vec k_point = average; k_point.area = A;
                    k_point.freq = s_val;
                    sum += integrand(k_point, q, w, T) * k_point.area / func_diff(k_point, q);
                }
            }
        }
    }
    return sum;
}

pair<int, int> get_index_and_length(double L, double U, vector<double> &sortedList) {
    int index = -1, length = 0;
    int lower_index = std::lower_bound(sortedList.begin(), sortedList.end(), L) - sortedList.begin();

    if (sortedList[lower_index] < L) return {-1, 0};

    for (int i = lower_index; i < sortedList.size() and sortedList[i] <= U; i++) {
        length = i - lower_index + 1;
    }
    return {lower_index, length};
}

double tetrahedron_sum_continuous(double (*func)(Vec k, Vec q), double (*func_diff)(Vec k, Vec q), Vec q, vector<double> &svals, double w, double T) {
    vector<vector<double>> tetrahedrons {
        {1, 2, 3, 5}, 
        {1, 3, 4, 5},
        {2, 5, 6, 3},
        {4, 5, 8, 3},
        {5, 8, 7, 3},
        {5, 6, 7, 3}
    };

    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n * (dim%2) + 1 * ((dim+1)%2); k++) {
                vector<VecAndEnergy> points = points_from_indices(func, q, i, j, k);
                double min = 1000, max = -1000;
                for (VecAndEnergy p : points) {
                    if (func(p.vec, q) < min) min = func(p.vec, q);
                    if (func(p.vec, q) > max) max = func(p.vec, q);
                }
                pair<int, int> index_and_length = get_index_and_length(min, max, svals);
                if (index_and_length.first == -1) continue;

                for (int c = 0; c < 6; c++) {

                    vector<VecAndEnergy> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }
                    for (int x = 0; x < index_and_length.second; x++) {
                        int ind = x + index_and_length.first;
                        double s_val = svals[ind];

                        if (not surface_inside_tetrahedron(s_val, ep_points)) continue;
                        vector<Vec> corner_points = points_in_tetrahedron(func, q, s_val, ep_points);

                        Vec average;

                        double b = 0;
                        if (corner_points[3] == average) b = 1.0;

                        for (Vec q : corner_points) {
                            average = (q + average);
                        }
                        average = average / (4-b);

                        double A = area_in_corners(corner_points);
                        if (dim == 2) A *= n / (2*k_max);
                        Vec k_point = average; k_point.area = A;
                        k_point.freq = s_val;

                        //double dS = 1;
                        double dS;
                        if (ind == 0) dS = (svals[ind+1] - svals[ind]) / 2;
                        else if (ind == svals.size()-1) dS = (svals[ind] - svals[ind-1]) / 2;
                        else dS = (svals[ind+1] - svals[ind-1]) / 2;
                        //dS = 0.1 / 2;
                        //printf("dS: %f\n", dS);

                        sum += integrand(k_point, q, w, T) 
                            * k_point.area / func_diff(k_point, q)
                            * dS;
                    }
                }
            }
        }
    }
    return sum;
}

double get_I(double D1, double D2, double D3, double V1, double V2, double V3, double V4) {
    if (V1 == V2 and V2 == V3 and V3 == V4 and V4 != 0) return 1/V1;
    if (V1 == V2 and V2 == V3 and V3 != V4 and V1 != 0) 
        return 3 * (V4*V4/pow(V1-V4,3)*log(V1/V4) + (1.5*V4*V4 + 0.5*V1*V1-2*V1*V4)/pow(V1-V4,3));
    if (V1 == V2 and V3 == V4 and V1 != V4) 
        return 3 * (2*V1*V4/pow(V1-V4,3)*log(V1/V4) + (V1+V4)/pow(V1-V4,2));
    if (V1 == V2 and V2 != V3 and V2 != V4 and V3 != V4)
        return 3 * (V3*V3/(pow(V3-V1,2)*(V3-V4)*log(V3/V1)) + V4*V4/(pow(V4-V1,2)*(V4-V3)) * log(V4/V1)
                + V1/((V3-V1)*(V4-V1)));
    return 3*(V1*V1/D1*log(V1/V4) + V2*V2/D2*log(V2/V4) + V3*V3/D3*log(V3/V4));
}

double analytic_tetrahedron_sum(Vec q, double w) {
    vector<vector<double>> tetrahedrons {
        {1, 2, 3, 5}, 
        {1, 3, 4, 5},
        {2, 5, 6, 3},
        {4, 5, 8, 3},
        {5, 8, 7, 3},
        {5, 6, 7, 3}
    };

    double sum = 0;
    double Omega = pow(2*k_max,3) / (6*n*n*n);
    if (dim == 2) Omega = pow(2*k_max,2) / (2*n*n);
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n * (dim%2) + 1 * ((dim+1)%2); k++) {
                vector<VecAndEnergy> points = points_from_indices(e_base_avg, q, i, j, k);

                for (int c = 0; c < 6; c++) {

                    vector<VecAndEnergy> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }
                    sort(ep_points.begin(), ep_points.end());

                    double V1 = e_diff(ep_points[3].vec,q) - w;
                    double V2 = e_diff(ep_points[2].vec,q) - w;
                    double V3 = e_diff(ep_points[1].vec,q) - w;
                    double V4 = e_diff(ep_points[0].vec,q) - w;

                    double D1 = (V1 - V4) * (V1 - V3) * (V1 - V2);
                    double D2 = (V2 - V4) * (V2 - V3) * (V2 - V1);
                    double D3 = (V3 - V4) * (V3 - V2) * (V3 - V1);

                    double I = get_I(D1, D2, D3, V1, V2, V3, V4);
                    sum += Omega * I;
                    if (isnan(sum)) {
                        cout << "NAN: " << V1 << " " << V2 << " " << V3 << " " << V4 << endl;
                        assert(false);
                    }
                }
            }
        }
    }
    return sum;
}

// Note: Should only be used for 2D, fails in 3D
vector<Vec> sort_by_adjacent(vector<Vec> &FS) {
    vector<Vec> sorted_FS;
    Vec empty;
    Vec k0 = FS[0];
    sorted_FS.push_back(k0);
    FS[0] = empty;
    while (sorted_FS.size() < FS.size()) {
        double min_dist = 100000;
        Vec min_vec;
        for (Vec k : FS) {
            if (k == empty) continue;
            double dist = (k - k0).vals.norm();
            if (dist < min_dist) {
                min_dist = dist;
                min_vec = k;
            }
        }
        sorted_FS.push_back(min_vec);
        cout << sorted_FS.size() << " ";
        k0 = min_vec;
        for (unsigned int i = 0; i < FS.size(); i++) {
            if (FS[i] == min_vec) {
                FS[i] = empty;
                break;
            }
        }
    }
    return sorted_FS;
}

vector<VecAndEnergy> FS_to_q_shifted(vector<Vec> &FS, Vec q) {
    vector<VecAndEnergy> q_shifted;
    for (Vec k : FS) {
        double energy = epsilon(k+q);
        q_shifted.push_back({k, energy});
    }
    return q_shifted;
}

void filter_FS(vector<Vec> &FS) {
    sort(FS.begin(), FS.end());

    double sum = 0; for (Vec k : FS) sum += k.area;

    double small_sum = 0; unsigned int i;
    for (i = 0; i < FS.size(); i++) {
        small_sum += FS[i].area;
        if (small_sum > 0.005*sum) break;
    }
    FS.erase(FS.begin(), FS.begin()+i);
}
