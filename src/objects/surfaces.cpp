/**
 * @file surfaces.cpp
 *
 * @brief This file contains all the aspects of the tetrahedron method that go into defining a 
 * constant energy contour surface. That function is general, but usually taken to be e(k) 
 * in this codebase
 *
 * @author Griffin Heier
 */

#include <math.h>
#include <boost/math/tools/roots.hpp>
#include <algorithm>
#include <functional>

#include "../config/load/cpp_config.h"
#include "vec.h"

using namespace std;

// Area contained within a triangle defined by these three points
float triangle_area_from_points(Vec k1, Vec k2, Vec k3) {
    auto triangle_area = [](float d1, float d2, float d3) { 
        float s = (d1 + d2 + d3)/2; 
        return pow(s*(s-d1)*(s-d2)*(s-d3), 0.5);
    };
    // Define distances 
    Vec k12 = k1 - k2; 
    Vec k23 = k2 - k3; 
    Vec k13 = k1 - k3; 
    // 0 index is radial distance
    float d12 = k12.norm();
    float d23 = k23.norm();
    float d13 = k13.norm();

    // Triangle area formula given side lengths
    float A = 0;
    A = triangle_area(d12, d23, d13);
    if (isnan(A)) {
        return 0;
    }

    return A;
}

vector<Vec> points_from_indices(function<float(Vec)> func, int i, int j, int k, int divs) {
    float x1 = 2*k_max * i / divs       - k_max; 
    float x2 = 2*k_max * (i+1) / divs   - k_max; 
    float y1 = 2*k_max * j / divs       - k_max; 
    float y2 = 2*k_max * (j+1) / divs   - k_max; 
    float z1 = 2*k_max * k / divs       - k_max; 
    float z2 = 2*k_max * (k+1) / divs   - k_max; 

    Vec p1(x1, y1, z1); p1.w = func(p1);  
    Vec p2(x2, y1, z1); p2.w = func(p2);
    Vec p3(x2, y2, z1); p3.w = func(p3);
    Vec p4(x1, y2, z1); p4.w = func(p4);
    Vec p5(x1, y1, z2); p5.w = func(p5);
    Vec p6(x2, y1, z2); p6.w = func(p6);
    Vec p7(x2, y2, z2); p7.w = func(p7);
    Vec p8(x1, y2, z2); p8.w = func(p8);

    vector<Vec> points(8); 
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    points[3] = p4;
    points[4] = p5;
    points[5] = p6;
    points[6] = p7;
    points[7] = p8;

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
vector<Vec> points_in_tetrahedron(function<float(Vec k)> func, float s_val, vector<Vec> points) {
    sort(points.begin(), points.end());
    Vec k1 = points[0], k2 = points[1], k3 = points[2], k4 = points[3];
    float ep1 = func(k1), ep2 = func(k2), ep3 = func(k3), ep4 = func(k4);

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
        if (s_val != func(points[i])) {
            times_not_equal++;
            not_equal = points[i];
        }
    if (times_not_equal == 1) {
        int iter = -1;
        for (int i = 0; i < 4; i++) {
            if (points[i] == not_equal) continue;
            return_points[iter] = points[i];
            iter++;
        }
        return_points[3] = empty;
    }

    return return_points;
}

// Same sign means the surface is not inside the cube
bool surface_inside_cube(float s_val, vector<Vec> p) {
    sort(p.begin(), p.end());
    return (p[7].w - s_val) / (p[0].w - s_val) < 0;
}

// Same sign means the surface is not inside the tetrahedron
bool surface_inside_tetrahedron(float s_val, vector<Vec> ep_points) {
    sort(ep_points.begin(), ep_points.end());
    return ((ep_points[3].w)-s_val) / ((ep_points[0].w) - s_val) < 0;
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
vector<Vec> tetrahedron_method(function<float(Vec k)> func, float s_val) {
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
                vector<Vec> points = points_from_indices(func, i, j, k, n);
                if (not surface_inside_cube(s_val, points)) continue;

                // Finds every k-space point in all possible tetrahedra, along with its associated
                // area in the surface, using the above functions
                for (int c = 0; c < 6; c++) {

                    vector<Vec> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }

                    if (not surface_inside_tetrahedron(s_val, ep_points)) continue;
                    vector<Vec> corner_points = points_in_tetrahedron(func, s_val, ep_points);

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
                    k_point.w = s_val;
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

