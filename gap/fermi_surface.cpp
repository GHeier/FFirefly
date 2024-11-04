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
#include <set>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <boost/math/tools/roots.hpp>
#include <algorithm>
#include <memory>
#include "cfg.h"
#include "vec.h"
#include "fermi_surface.h"
#include "potential.h"
#include "frequency_inclusion.hpp"
#include "surfaces.h"
#include "integration.h"

using namespace std;

// This is the method that defines the surface; It is the culmination and the point of this file
vector<Vec> tetrahedron_method(float (*func)(Vec k, Vec q), Vec q, float s_val) {
    auto func2 = [func, q](Vec k) { return func(k, q); };
    return tetrahedron_method2(func2, s_val);
}

// Same algorithm as above, but instead of defining a surface, we simply sum across all of them
// Sum is done over the function "integrand"
float tetrahedron_sum_continuous(float (*func)(Vec k, Vec q), float (*func_diff)(Vec k, Vec q), Vec q, vector<float> &svals, float w, float T) {
    auto int_numerator = [q, w, T](Vec k) -> float {
        return ratio(k, q, w, T);
    };
    auto func2 = [func, q](Vec k) { return func(k, q); };
    auto func_diff2 = [func_diff, q](Vec k) { return func_diff(k, q); };
    return surface_transform_integral(int_numerator, func2, func_diff2, svals);
}

void sanitize_I_vals(float &V1, float &V2, float &V3, float &V4) {
    if (fabs(V1 - V2) < 1e-3) V2 = V1;
    if (fabs(V1 - V3) < 1e-3) V3 = V1;
    if (fabs(V1 - V4) < 1e-3) V4 = V1;
    if (fabs(V2 - V3) < 1e-3) V3 = V2;
    if (fabs(V2 - V4) < 1e-3) V4 = V2;
    if (fabs(V3 - V4) < 1e-3) V4 = V3;
}

vector<float> getUnique(float a, float b, float c, float d) {
    // Use a set to find unique values 
    set<float, greater<float>> uniqueValues = {a, b, c, d};
    // Copy the sorted unique values to a vector
    vector<float> result(uniqueValues.begin(), uniqueValues.end());
    return result;
}

bool check_two_equal(float V1, float V2, float V3, float V4) {
    return (V1 == V2 and V3 == V4) or (V1 == V3 and V2 == V4) or (V1 == V4 and V2 == V3);
}

// Integral value of the tetrahedron method when interpolated linearly across each small cube
float get_I(float D1, float D2, float D3, float V1, float V2, float V3, float V4) {
    sanitize_I_vals(V1, V2, V3, V4);
    vector<float> V = getUnique(V1, V2, V3, V4);
    if (find(V.begin(), V.end(), 0) != V.end()) {
        return 0;
    }
    if (V.size() == 1) {
        float r = 1 / V[0];
        return r;
    }
    if (V.size() == 2 and check_two_equal(V1, V2, V3, V4)) {
        float t1 = 2 * V[0] * V[1] / pow(V[0] - V[1], 3) * log(fabs(V[1] / V[0]));
        float t2 = (V[0] + V[1]) / (pow(V[0] - V[1], 2));
        float r = 3 * (t1 + t2);
        return 3 * (t1 + t2);
    }
    else if (V.size() == 2) {
        float t1 = V[1]*V[1] / (pow(V[0] - V[1],3) * log(fabs(V[0]/V[1])));
        float t2 = (1.5 * V[1]*V[1] + 0.5 * V[0]*V[0] - 2 * V[0]*V[1]) / pow(V[0] - V[1],3);
        float r = 3 * (t1 + t2);
        return r;
    }
    if (V.size() == 3) {
        float t1 = V[1] * V[1] / (pow(V[1] - V[0], 2) * (V[1] - V[2])) * log(fabs(V[1] / V[0]));
        float t2 = V[2] * V[2] / (pow(V[2] - V[0], 2) * (V[2] - V[1])) * log(fabs(V[2] / V[0]));
        float t3 = V[0] / ((V[1] - V[0]) * (V[2] - V[0]));
        float r = 3 * (t1 + t2 + t3);
        return r;
    }

    float t1 = (V1*V1/D1*log(fabs(V1/V4)));
    float t2 = (V2*V2/D2*log(fabs(V2/V4)));
    float t3 = (V3*V3/D3*log(fabs(V3/V4)));
    float r = 3 * (t1 + t2 + t3);
    return r;
}

// Computes the sum analytically, which should be quite a bit faster
// Done at zero temperature is the only caveat
float analytic_tetrahedron_sum(Vec q, float w, int num_pts) {
    vector<vector<float>> tetrahedrons {
        {1, 2, 3, 5}, 
        {1, 3, 4, 5},
        {2, 5, 6, 3},
        {4, 5, 8, 3},
        {5, 8, 7, 3},
        {5, 6, 7, 3}
    };

    float sum = 0;
    float Omega = pow(2*k_max,3) / (6*num_pts*num_pts*num_pts);
    if (dim == 2) Omega = pow(2*k_max,2) / (2*n*n);
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < num_pts; i++) {
        for (int j = 0; j < num_pts; j++) {
            for (int k = 0; k < num_pts * (dim%2) + 1 * ((dim+1)%2); k++) {
                vector<Vec> points = points_from_indices2(epsilon, i, j, k, num_pts);

                for (int c = 0; c < 6; c++) {

                    vector<Vec> ep_points(4);
                    for (int p = 0; p < 4; p++) {
                        ep_points[p] = points[tetrahedrons[c][p]-1];
                    }
                    // Sort by freq size
                    sort(ep_points.begin(), ep_points.end(), [](Vec a, Vec b) { 
                            return a.freq < b.freq; 
                    });

                    float V1 = e_diff(ep_points[3],q) - w;
                    float V2 = e_diff(ep_points[2],q) - w;
                    float V3 = e_diff(ep_points[1],q) - w;
                    float V4 = e_diff(ep_points[0],q) - w;

                    float D1 = (V1 - V4) * (V1 - V3) * (V1 - V2);
                    float D2 = (V2 - V4) * (V2 - V3) * (V2 - V1);
                    float D3 = (V3 - V4) * (V3 - V2) * (V3 - V1);

                    float I = get_I(D1, D2, D3, V1, V2, V3, V4);
                    sum += Omega * I;
                }
            }
        }
    }
    assert(not isnan(sum));
    return sum;
}
