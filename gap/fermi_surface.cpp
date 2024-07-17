/**
 * @file fermi_surface.cpp
 *
 * @brief This file contains all the aspects of the tetrahedron method that go into defining a 
 * constant energy contour surface. That function is general, but usually taken to be e(k) 
 * in this codebase
 *
 * @author Griffin Heier
 */

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <boost/math/tools/roots.hpp>
#include <algorithm>
#include "cfg.h"
#include "vec.h"
//#include "frequency_inclusion.hpp"
#include "fermi_surface.h"
#include "potential.h"
#include "frequency_inclusion.hpp"

using namespace std;

bool operator<(const VecAndEnergy& left, const VecAndEnergy& right) {
    return left.energy < right.energy;
}

// Area contained within a triangle defined by these three points
float triangle_area_from_points(Vec k1, Vec k2, Vec k3) {
    auto triangle_area = [](float d1, float d2, float d3) { 
        float s = (d1 + d2 + d3)/2; 
        return pow(s*(s-d1)*(s-d2)*(s-d3), 0.5);
    };
    // Define distances 
    Vec k12 = k1 - k2; if (k12.cartesian) k12.to_spherical();
    Vec k23 = k2 - k3; if (k23.cartesian) k23.to_spherical();
    Vec k13 = k1 - k3; if (k13.cartesian) k13.to_spherical();
    // 0 index is radial distance
    float d12 = k12.vals[0];
    float d23 = k23.vals[0];
    float d13 = k13.vals[0];

    // Triangle area formula given side lengths
    float A = 0;
    A = triangle_area(d12, d23, d13);
    if (isnan(A)) {
        return 0;
    }

    return A;
}

/**
 * @brief Calculates the 8 points of a cube in k-space from a set of indices
 *
 * @param func The function to evaluate that defines the value of a surface
 * @param q The q vector
 * @param i The first index
 * @param j The second index
 * @param k The third index
 * @param divs The number of divisions in the k-space
 *
 * @return A vector of 8 VecAndEnergy structs
 */
vector<VecAndEnergy> points_from_indices(float (*func)(Vec k, Vec q), Vec q, int i, int j, int k, int divs) {
    float x1 = 2*k_max * i / divs       - k_max; 
    float x2 = 2*k_max * (i+1) / divs   - k_max; 
    float y1 = 2*k_max * j / divs       - k_max; 
    float y2 = 2*k_max * (j+1) / divs   - k_max; 
    float z1 = 2*k_max * k / divs       - k_max; 
    float z2 = 2*k_max * (k+1) / divs   - k_max; 

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

/**
 * @brief Picks out the points that create the surface inside of a tetrahedra
 *
 * The surface is made up of a plane inside of a tetrahedron. This function approximates the 
 * surface as linear, and then extrapolates to find the points that are on the surface of the 
 * tetrahedron. These points make up the plane of the surface at the value s_val.
 *
 * @param func The function to evaluate that defines the value of a surface
 * @param q The q vector
 * @param s_val The value of the surface
 * @param points The points of the tetrahedron
 *
 * @return A vector of Vec structs, sorted by energy
 */
vector<Vec> points_in_tetrahedron(float (*func)(Vec k, Vec q), Vec q, float s_val, vector<VecAndEnergy> points) {
    sort(points.begin(), points.end());
    Vec k1 = points[0].vec, k2 = points[1].vec, k3 = points[2].vec, k4 = points[3].vec;
    float ep1 = func(k1, q), ep2 = func(k2, q), ep3 = func(k3, q), ep4 = func(k4, q);

    Vec empty;

    // y = m * x + b to find the points on the surface
    Vec k12 = (k2-k1) * (s_val - ep1) / (ep2 - ep1) + k1;
    Vec k13 = (k3-k1) * (s_val - ep1) / (ep3 - ep1) + k1;
    Vec k14 = (k4-k1) * (s_val - ep1) / (ep4 - ep1) + k1;
    Vec k24 = (k4-k2) * (s_val - ep2) / (ep4 - ep2) + k2;
    Vec k34 = (k4-k3) * (s_val - ep3) / (ep4 - ep3) + k3;
    Vec k23 = (k3-k2) * (s_val - ep2) / (ep3 - ep2) + k2;
    
    vector<Vec> return_points(4, empty);


    // Assigns which type of surface is chosen based on which points the surface lies between
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

// Same sign means the surface is not inside the cube
bool surface_inside_cube(float s_val, vector<VecAndEnergy> p) {
    sort(p.begin(), p.end());
    return (p[7].energy - s_val) / (p[0].energy - s_val) < 0;
}

// Same sign means the surface is not inside the tetrahedron
bool surface_inside_tetrahedron(float s_val, vector<VecAndEnergy> ep_points) {
    sort(ep_points.begin(), ep_points.end());
    return ((ep_points[3].energy)-s_val) / ((ep_points[0].energy) - s_val) < 0;
}

float area_in_corners(vector<Vec> cp) {
    Vec empty;
    Vec k1 = cp[0]; Vec k2 = cp[1]; Vec k3 = cp[2]; Vec k4 = cp[3];
    if (k4 == empty) return triangle_area_from_points(k1, k2, k3);

    float A1 = 0, A2 = 0;
    A1 = triangle_area_from_points(k1, k2, k4);
    A2 = triangle_area_from_points(k3, k2, k4);

    return A1 + A2;
}

// This is the method that defines the surface; It is the culmination and the point of this file
vector<Vec> tetrahedron_method(float (*func)(Vec k, Vec q), Vec q, float s_val) {
    vector<vector<float>> tetrahedrons {
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
                vector<VecAndEnergy> points = points_from_indices(func, q, i, j, k, n);
                if (not surface_inside_cube(s_val, points)) continue;

                // Finds every k-space point in all possible tetrahedra, along with its associated
                // area in the surface, using the above functions
                for (int c = 0; c < 6; c++) {

                    vector<VecAndEnergy> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }

                    if (not surface_inside_tetrahedron(s_val, ep_points)) continue;
                    vector<Vec> corner_points = points_in_tetrahedron(func, q, s_val, ep_points);

                    // Averages the corner points to find the center of the triangle
                    Vec average;
                    float b = 0;
                    if (corner_points[3] == average) b = 1.0;

                    for (Vec q : corner_points) {
                        average = (q + average);
                    }
                    average = average / (4-b);

                    float A = area_in_corners(corner_points);
                    if (dim == 2) A *= n / (2*k_max);
                    Vec k_point = average; k_point.area = A;
                    k_point.freq = s_val;
                    FS.push_back(k_point);
                }
            }
        }
    }
    return FS;
}

// -1,0 is returned if there is no surface in the cube
// length is the number of surfaces in the cube
pair<int, int> get_index_and_length(float L, float U, vector<float> &sortedList) {
    int index = -1, length = 0;
    // Binary search for lower index
    int lower_index = std::lower_bound(sortedList.begin(), sortedList.end(), L) - sortedList.begin();

    // Linear search for upper index
    if (sortedList[lower_index] < L) return {-1, 0};

    for (int i = lower_index; i < sortedList.size() and sortedList[i] <= U; i++) {
        length = i - lower_index + 1;
    }
    return {lower_index, length};
}

// Same algorithm as above, but instead of defining a surface, we simply sum across all of them
// Sum is done over the function "integrand"
float tetrahedron_sum_continuous(float (*func)(Vec k, Vec q), float (*func_diff)(Vec k, Vec q), Vec q, vector<float> &svals, float w, float T) {
    vector<vector<float>> tetrahedrons {
        {1, 2, 3, 5}, 
        {1, 3, 4, 5},
        {2, 5, 6, 3},
        {4, 5, 8, 3},
        {5, 8, 7, 3},
        {5, 6, 7, 3}
    };

    float sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < s_div; i++) {
        for (int j = 0; j < s_div; j++) {
            for (int k = 0; k < s_div * (dim%2) + 1 * ((dim+1)%2); k++) {
                vector<VecAndEnergy> points = points_from_indices(func, q, i, j, k, s_div);
                float min = 1000, max = -1000;
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
                        float s_val = svals[ind];

                        if (not surface_inside_tetrahedron(s_val, ep_points)) continue;
                        vector<Vec> corner_points = points_in_tetrahedron(func, q, s_val, ep_points);

                        Vec average;

                        float b = 0;
                        if (corner_points[3] == average) b = 1.0;

                        for (Vec q : corner_points) {
                            average = (q + average);
                        }
                        average = average / (4-b);

                        float A = area_in_corners(corner_points);
                        if (dim == 2) A *= s_div / (2*k_max);
                        Vec k_point = average; k_point.area = A;
                        k_point.freq = s_val;

                        // Surfaces being summed over. Depending on whether the surface is 
                        // at the beginning of the end of the list, we weight it differently
                        // This is trapezoidal integration
                        float dS;
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

// Integral value of the tetrahedron method when interpolated linearly across each small cube
float get_I(float D1, float D2, float D3, float V1, float V2, float V3, float V4) {
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

// Computes the sum analytically, which should be quite a bit faster
// Done at zero temperature is the only caveat
float analytic_tetrahedron_sum(Vec q, float w) {
    vector<vector<float>> tetrahedrons {
        {1, 2, 3, 5}, 
        {1, 3, 4, 5},
        {2, 5, 6, 3},
        {4, 5, 8, 3},
        {5, 8, 7, 3},
        {5, 6, 7, 3}
    };

    float sum = 0;
    float Omega = pow(2*k_max,3) / (6*n*n*n);
    if (dim == 2) Omega = pow(2*k_max,2) / (2*n*n);
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n * (dim%2) + 1 * ((dim+1)%2); k++) {
                vector<VecAndEnergy> points = points_from_indices(e_base_avg, q, i, j, k, n);

                for (int c = 0; c < 6; c++) {

                    vector<VecAndEnergy> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }
                    sort(ep_points.begin(), ep_points.end());

                    float V1 = e_diff(ep_points[3].vec,q) - w;
                    float V2 = e_diff(ep_points[2].vec,q) - w;
                    float V3 = e_diff(ep_points[1].vec,q) - w;
                    float V4 = e_diff(ep_points[0].vec,q) - w;

                    float D1 = (V1 - V4) * (V1 - V3) * (V1 - V2);
                    float D2 = (V2 - V4) * (V2 - V3) * (V2 - V1);
                    float D3 = (V3 - V4) * (V3 - V2) * (V3 - V1);

                    float I = get_I(D1, D2, D3, V1, V2, V3, V4);
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

