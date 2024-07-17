/**
 * @file susceptibility.cpp
 *
 * @brief Calculates the susceptibility of the system.
 *
 * @author Griffin Heier
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>

#include <omp.h>
#include <boost/functional/hash.hpp>
#include <tuple>
#include <unordered_map>

#include "utilities.h"
#include "potential.h"
#include "vec.h"
#include "fermi_surface.h"
#include "band_structure.h"
#include "cfg.h"
#include "frequency_inclusion.hpp"
#include "susceptibility.h"

using namespace std;

float f(float E, float T) {
    if (T == 0) {
        if (E < 0) return 1;
        return 0;
    }
    return 1 / (1 + exp(E/T));
}

float integrand(Vec k, Vec q, float w, float T) {
    float dE = k.freq;
    float e_k = epsilon(k) - mu;
    //float e_qk = epsilon(k+q) - MU;
    float e_qk = e_k + dE;
    float f_k = f(e_k, T);
    float f_qk = f(e_qk, T);


    if (fabs(dE) < 0.0001 and fabs(w) < 0.0001) {
        if (T == 0 or exp(e_k/T) > 1e6)
            return e_k < 0;
        float temp = 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        return temp;
    }
    if (fabs(w - dE) == 0) return 0;
    return (f_qk - f_k) / (w - dE);
}

float integrate_susceptibility(Vec q, float T, float mu, float w, int num_points) {
    //return imaginary_integration(q, T, mu, w, num_points, 0.000);
    if (q.norm() < 0.0001) {
        vector<Vec> FS = tetrahedron_method(e_base_avg, q, mu);
        float sum = 0; for (auto x : FS) sum += x.area / vp(x);
        return sum / pow(2*k_max,dim);
    }
    if (q.norm() < 0.0001) q = Vec(0.01,0.01,0.01);

    float a, b; get_bounds(q, b, a, denominator);
    vector<float> spacing; get_spacing_vec(spacing, w, a, b, num_points);
    float c = tetrahedron_sum_continuous(denominator, denominator_diff, q, spacing, w, T);
    return c / pow(2*k_max,dim);
}

float trapezoidal_integration(auto &f, float x0, float x1, float y0, float y1, float z0, float z1, int num_points) {
    float sum = 0;
    float dx = (x1 - x0) / (num_points - 1);
    float dy = (y1 - y0) / (num_points - 1);
    float dz = (z1 - z0) / (num_points - 1);
    if (dim == 2) dz = 1;
    int counter = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < num_points; i++) {
        float x = x0 + i*dx;
        for (int j = 0; j < num_points; j++) {
            float y = y0 + j*dy;
            for (float k = 0; k < num_points * (dim%2) + 1 * ((dim+1)%2); k++) {
                float z = z0 + k*dz;
                float w = 1.0;
                if (i == 0 or i == num_points - 1) w /= 2.0;
                if (j == 0 or j == num_points - 1) w /= 2.0;
                if ( (k == 0 or k == num_points - 1) and dim == 3) w /= 2.0;

                //if (x*x + y*y + z*z > k_max*k_max) continue;
                //cout << x << " " << y << " " << z << " " << f(x,y,z) << endl;
                sum += w * f(x,y,z);
            }
        }
    }
    return sum * dx * dy * dz;
}

vector<vector<vector<float>>> chi_cube(float T, float mu, float w, string message) {
    int m_z = m*(dim%2) + 3*((dim+1)%2);
    vector<vector<vector<float>>> cube(m, vector<vector<float>> (m, vector<float> (m_z)));
    unordered_map<string, float> map;
    float empty_val = -98214214;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < m_z; k++) {
                Vec q((2*k_max*i)/(m-1), (2*k_max*j)/(m-1), (2*k_max*k)/(m_z-1));
                Vec q2 = to_IBZ_2(q);
                if (map.find(vec_to_string(q2)) == map.end())
                    map[vec_to_string(q2)] = empty_val;
            }
        }
    }
    
    //cout << "Taking " << map.size() << " integrals in " << dim << " dimensions.\n";
    //#pragma omp parallel for
    for(unsigned int i = 0; i < map.size(); i++) {
        auto datIt = map.begin();
        advance(datIt, i);
        string key = datIt->first;
        map[key] = integrate_susceptibility(string_to_vec(key), T, mu, w, s_pts);
        progress_bar(1.0 * i / (map.size()-1), message);
    }
    cout << endl;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < m_z; k++) {
                Vec q((2*k_max*i)/(m-1), (2*k_max*j)/(m-1), (2*k_max*k)/(m_z-1));
                Vec q2 = to_IBZ_2(q);
                cube[i][j][k] = map[vec_to_string(q2)];
            }
        }
    }

    return cube;
}

float calculate_chi_from_cube(const vector<vector<vector<float>>> &chi_cube, Vec q) {
    Vec v = to_IBZ_2(q);
    float d = 2*k_max/(m-1);

    float x = v.vals[0], y = v.vals[1], z = v.vals[2];
    if (dim == 2) z = 0;

    int i = floor(x / d);
    int j = floor(y / d);
    int k = floor(z / d);

    float x1 = i * d; 
    float y1 = j * d; 
    float z1 = k * d; 

    float x2 = x1 + d; 
    float y2 = y1 + d; 
    float z2 = z1 + d; 

    float dx = 0, dy = 0, dz = 0, wx = 0, wy = 0, wz = 0, w0 = 0;

    // Make sure there's no issue with indexing
    //cout << q << q.vals[2] << endl;
    //int s = chi_cube.size()-1; 
    //assert( i < s and j < s and k < chi_cube[0][0].size()-1);

    float f1 = chi_cube[i][j][k], f2 = chi_cube[i+1][j][k];
    float f3 = chi_cube[i+1][j+1][k], f4 = chi_cube[i][j+1][k];
    float f5 = chi_cube[i][j][k+1], f6 = chi_cube[i+1][j][k+1];
    float f7 = chi_cube[i+1][j+1][k+1], f8 = chi_cube[i][j+1][k+1];

    if (x - x1 <= z2 - z and x - x1 >= y - y1) {// blue @ 1
        w0 = f1;
        wx = (f2 - f1) / d; dx = x - x1;
        wy = (f3 - f2) / d; dy = y - y1;
        wz = (f5 - f1) / d; dz = z - z1;
    }

    else if (y + z <= z1 + y2 and x - x1 <= y - y1) {// orange @ 1
        w0 = f1;
        wx = (f3 - f4) / d; dx = x - x1;
        wy = (f4 - f1) / d; dy = y - y1;
        wz = (f5 - f1) / d; dz = z - z1;
    }

    else if (x + z >= z1 + x2 and y + z <= z1 + y2) {// red @ 2
        w0 = f2;
        wx = (f6 - f5) / d; dx = x - x2;
        wy = (f3 - f2) / d; dy = y - y1;
        wz = (f6 - f2) / d; dz = z - z1;
    }

    else if (y + z >= z1 + y2 and x + z <= z1 + x2) {// purple @ 4
        w0 = f4;
        wx = (f3 - f4) / d; dx = x - x1;
        wy = (f8 - f5) / d; dy = y - y2;
        wz = (f8 - f4) / d; dz = z - z1;
    }

    else if (x - x1 >= y - y1 and y + z >= z1 + y2) {// teal @ 7
        w0 = f7;
        wx = (f6 - f5) / d; dx = x - x2;
        wy = (f7 - f6) / d; dy = y - y2;
        wz = (f7 - f3) / d; dz = z - z2;
    }

    else if (x - x1 <= y - y1 and y + z >= z1 + y2) {// green @ 7
        w0 = f7;
        wx = (f7 - f8) / d; dx = x - x2;
        wy = (f8 - f5) / d; dy = y - y2;
        wz = (f7 - f3) / d; dz = z - z2;
    }
    else return f1 + (f2-f1)/d*x + (f4-f1)/d*y + (f5-f1)/d*z;

    return w0 + wx*dx + wy*dy + wz*dz;
}

Vec to_IBZ_2(const Vec k) {
    Vec q = k;
    if (q.cartesian == false) q.to_cartesian();
    float x = q.vals[0], y = q.vals[1], z = q.vals[2];
    x = abs(x); y = abs(y); z = abs(z);
    if (x > M_PI) x = - (x - 2*M_PI);
    if (y > M_PI) y = - (y - 2*M_PI);
    if (z > M_PI) z = - (z - 2*M_PI);
    if (dim == 3) {
        float arr[] = {x, y, z};
        sort(arr, arr+3, greater<float>());
        auto& [a, b, c] = arr;
        Vec result(a, b, c);
        return result;
    }
    else if (dim == 2) {
        float arr[] = {x, y};
        sort(arr, arr+2, greater<float>());
        auto& [a, b] = arr;
        Vec result(a, b, z);
        return result;
    }
    else {
        cout << "Wrong Dimension\n";
        return q;
    }
}

