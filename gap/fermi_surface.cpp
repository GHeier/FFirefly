#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <boost/math/tools/roots.hpp>
#include <algorithm>
#include "cfg.h"
#include "vec.h"
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
    double d12 = k12.vals[0];
    double d23 = k23.vals[0];
    double d13 = k13.vals[0];

    // Calculate areas if the triangle is on the fermi surface
    double A = 0;
    A = triangle_area(d12, d23, d13);
    if (isnan(A)) {
        return 0;
    }

    return A;
}

vector<VecAndEnergy> points_from_indices(int i, int j, int k) {
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
    VecAndEnergy point1 = {p1, epsilon(p1)}; points[0] = point1;
    VecAndEnergy point2 = {p2, epsilon(p2)}; points[1] = point2;
    VecAndEnergy point3 = {p3, epsilon(p3)}; points[2] = point3;
    VecAndEnergy point4 = {p4, epsilon(p4)}; points[3] = point4;
    VecAndEnergy point5 = {p5, epsilon(p5)}; points[4] = point5;
    VecAndEnergy point6 = {p6, epsilon(p6)}; points[5] = point6;
    VecAndEnergy point7 = {p7, epsilon(p7)}; points[6] = point7;
    VecAndEnergy point8 = {p8, epsilon(p8)}; points[7] = point8;

    return points;
}

vector<Vec> points_in_tetrahedron(double mu, vector<VecAndEnergy> points) {
    sort(points.begin(), points.end());
    Vec k1 = points[0].vec, k2 = points[1].vec, k3 = points[2].vec, k4 = points[3].vec;
    double ep1 = epsilon(k1), ep2 = epsilon(k2), ep3 = epsilon(k3), ep4 = epsilon(k4);

    Vec empty;

    Vec k12 = (k2-k1) * (mu - ep1) / (ep2 - ep1) + k1;
    Vec k13 = (k3-k1) * (mu - ep1) / (ep3 - ep1) + k1;
    Vec k14 = (k4-k1) * (mu - ep1) / (ep4 - ep1) + k1;
    Vec k24 = (k4-k2) * (mu - ep2) / (ep4 - ep2) + k2;
    Vec k34 = (k4-k3) * (mu - ep3) / (ep4 - ep3) + k3;
    Vec k23 = (k3-k2) * (mu - ep2) / (ep3 - ep2) + k2;
    
    vector<Vec> return_points(4, empty);


    if ( mu > ep1 and mu <= ep2) {
        return_points[0] = k12;
        return_points[1] = k13;
        return_points[2] = k14;
    }
    if ( mu > ep3 and mu <= ep4) {
        return_points[0] = k14;
        return_points[1] = k24;
        return_points[2] = k34;
    }
    if ( mu > ep2 and mu <= ep3) {
        return_points[0] = k24;
        return_points[1] = k23;
        return_points[2] = k13;
        return_points[3] = k14;
    }

    // Run condition for when plane aligns with one of the planes of the tetrahedron
    int times_not_equal = 0;
    Vec not_equal;
    for (int i = 0; i < 4; i++)
        if (mu != epsilon(points[i].vec)) {
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

bool surface_inside_cube(double mu, vector<VecAndEnergy> p) {
    sort(p.begin(), p.end());
    return (p[7].energy-mu) / (p[0].energy - mu) < 0;
}

bool surface_inside_tetrahedron(double mu, vector<VecAndEnergy> ep_points) {
    sort(ep_points.begin(), ep_points.end());
    return ((ep_points[3].energy)-mu) / ((ep_points[0].energy) - mu) < 0;
}

double area_inside(vector<double> epsilons) {
    sort(epsilons.begin(), epsilons.end());
    if ( mu > epsilons[0] and mu <= epsilons[1] ) {
        double ep01 = mu - epsilons[0];
        double ep21 = epsilons[1] - epsilons[0];
        double ep31 = epsilons[2] - epsilons[0];
        double ep41 = epsilons[3] - epsilons[0];
        return 0.5 * pow(ep01, 2) * pow( pow(ep21,2) + pow(ep31,2) + pow(ep41,2), 0.5) / (ep21*ep31*ep41);
    }
    if ( mu > epsilons[1] and mu <= epsilons[2] ) {
        double ep20 = mu - epsilons[1];
        double ep21 = epsilons[1] - epsilons[0];
        double ep31 = epsilons[2] - epsilons[0];
        double ep41 = epsilons[3] - epsilons[0];
        double ep32 = epsilons[2] - epsilons[1];
        double ep42 = epsilons[3] - epsilons[1];
        return 0.5 * pow( pow(ep21,2) + pow(ep31,2) + pow(ep41,2), 0.5) / (ep31*ep41) * (ep21 - 2*ep20 - pow(ep20,2) * (ep31 + ep42) / (ep32 * ep42));
    }
    if ( mu > epsilons[2] and mu <= epsilons[3] ) {
        double ep40 = epsilons[3] - mu;
        double ep21 = epsilons[1] - epsilons[0];
        double ep31 = epsilons[2] - epsilons[0];
        double ep41 = epsilons[3] - epsilons[0];
        double ep42 = epsilons[3] - epsilons[1];
        double ep43 = epsilons[3] - epsilons[2];
        return 0.5 * pow(ep40, 2) * pow( pow(ep21,2) + pow(ep31,2) + pow(ep41,2), 0.5) / (ep41*ep42*ep43);
    }
    cout << "Not inside of tetrahedron\n";
    assert(0==1);
    return 0;
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

vector<Vec> tetrahedron_method(double mu) {
    double surface_area = 0;
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
                vector<VecAndEnergy> points = points_from_indices(i, j, k);
                if (not surface_inside_cube(mu, points)) continue;

                for (int c = 0; c < 6; c++) {

                    vector<VecAndEnergy> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }

                    if (not surface_inside_tetrahedron(mu, ep_points)) continue;
                    vector<Vec> corner_points = points_in_tetrahedron(mu, ep_points);

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
                    //if (A > 0.0001)
                        FS.push_back(k_point);
                    surface_area += A;
                    assert(not isnan(surface_area));
                }
            }
        }
    }
    return FS;
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
            double dist = (k - k0).norm();
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
