/**
 * @file frequency_inclusion.cpp
 *
 * @brief This file is used to calculate the frequency dependent susceptibility. 
 *
 * @details It includes all the modifications taken to the standard BCS approach in order to allow 
 * for frequency-dependent interactions that occur some energy away from the fermi surface, denoted 
 * by the cutoff frequency, wc
 *
 * @author Griffin Heier
 */

#include <iostream>
#include <unistd.h>
#include <ctime>

#include <string>
#include <math.h>
#include <complex>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/gauss.hpp>

#include "../config/load/cpp_config.h"
#include "../algorithms/integration.h"
#include "../objects/vec.h"
#include "utilities.h"
#include "frequency_inclusion.hpp"
#include "../response/susceptibility.h"
#include "../hamiltonian/band_structure.h"

using std::cout;
using std::endl;
using std::sort;
using std::vector;
using std::unordered_map;
//using lambda_lanczos::LambdaLanczos;

// Gets total matrix size

// Defines all energy surfaces around_val FS
vector<vector<Vec>> freq_tetrahedron_method(float MU) {
    assert( l % 2 != 0); // N must be odd that way frequencies are evenly spaced
    vector<vector<Vec>> basis;

    float points_0th[1] = {0}; float *p0 = points_0th;
    float points_1st[2] = {-1/pow(3,0.5), 1/pow(3,0.5)}; float *p1 = points_1st;
    float points_2nd[3] = {-pow(3/5,0.5), 0, pow(3/5,0.5)}; float *p2 = points_2nd;
    float points_3rd[4] = {-0.861136, -0.339981, 0.339981, 0.861136}; float *p3 = points_3rd;
    float points_4th[5] = {-0.90618, -0.538469, 0, 0.538469, 0.90618}; float *p4 = points_4th;

    float *points[5] = {p0, p1, p2, p3, p4};

    for (int i = 0; i < l; i++) {
        float ep = wc * points[l-1][i] + MU;
        vector<Vec> layer = get_FS(ep);
        basis.push_back(layer);
    }
    return basis;
}


// Creates a map of multiple "chi cubes" that are used to calculate the frequency dependent 
// susceptibility
unordered_map <float, vector<vector<vector<float>>>> chi_cube_freq(float T, float MU) {
    vector<float> des;
    for (int i = 0; i < l; i++) {
        float p1 = wc * points[l-1][i];
        for (int j = 0; j < l; j++) {
            float p2 = wc * points[l-1][j];
            float w = p1 - p2;
            w = round_val(w, 6);
            if (count(des.begin(), des.end(), w) == 0)
                des.push_back(w);
        }
    }

    unordered_map <float, vector<vector<vector<float>>>> cube_freq_map;
    for (int i = 0; i < des.size(); i++) {
        string message = "Chi Cube " + to_string(i+1) + " / " + to_string(des.size());
        float w = des[i];
        auto cube = chi_cube(T, MU, w, message);
        cube_freq_map.insert(pair<float, vector<vector<vector<float>>>>(w, cube));
    }
//    for (auto x : cube_freq_map) {
//        cout << "Frequency: " << x.first << endl;
//    }
    return cube_freq_map;
}

MatCube create_matsubara_cube(float T, float MU, int m_pts, int w_pts, float w_min, float w_max, int num_integral_pts) {
    // Create the Matsubara cube
    // The Matsubara cube is a 4D cube that contains the susceptibility at each point in k-space and frequency space
    // The cube is created by integrating the susceptibility over the entire BZ
    // The cube is then used to interpolate the susceptibility at each point in the BZ
    // The cube is then used to calculate the V matrix
    int m_z = m * (dim%2) + 1*((dim+1)%2); // If dim == 2, m_z = 1, if dim == 3, m_z = m
    MatCube matsubara_cube(m_pts, m_pts, m_z, w_pts, 0, k_max, 0, k_max, 0, k_max, w_min, w_max, num_integral_pts);
    // First determine irreducible BZ points
    unordered_map<string, complex<float>> map;
    float empty_val = -98214214.0;
    for (int i = 0; i < m_pts; i++) {
        for (int j = 0; j < m_pts; j++) {
            for (int k = 0; k < m_z; k++) {
                Vec q((k_max*i)/(1.0*m_pts-1), (k_max*j)/(1.0*m_pts-1), (k_max*k)/(1.0*m_z-1));
                if (m_z == 1) q = Vec(2*k_max*i/(m_pts-1) , 2*k_max*j/(m_pts-1), 0);
                for (int l = 0; l < w_pts; l++) {
                    float w = w_min + l * (w_max - w_min) / (w_pts - 1);
                    q = to_IBZ_2(q);
                    string key = vec_to_string(q) + " " + to_string(w);
                    if (map.find(key) == map.end())
                        map[key] = empty_val;
                }
            }
        }
    }
    // Now take the integrals at these points
    cout << "Taking " << map.size() << " integrals in " << dim << "+1 dimensions.\n";
    //#pragma omp parallel for
    for(unsigned int i = 0; i < map.size(); i++) {
        auto datIt = map.begin();
        advance(datIt, i);
        string key = datIt->first;
        vector<float> split_key = unpack_string(key);
        Vec q(split_key[0], split_key[1], split_key[2]);
        complex<float> w = complex<float>(0, split_key[3]);
        map[key] = complex_susceptibility_integration(q, T, MU, w, num_integral_pts);
        progress_bar(1.0 * i / (map.size()-1), "MatCube Creation");
    }
    cout << endl;

    // Now fill the cube with the values. The cube is filled over the entirety of the BZ, not just the IBZ. This is out of laziness
    for (int i = 0; i < m_pts; i++) {
        for (int j = 0; j < m_pts; j++) {
            for (int k = 0; k < m_z; k++) {
                for (int l = 0; l < w_pts; l++) {
                    Vec q((k_max*i)/(m_pts-1), (k_max*j)/(m_pts-1), (k_max*k)/(m_z-1));
                    if (m_z == 1) q = Vec(k_max*i/(m-1) , k_max*j/(m-1), 0);
                    q = to_IBZ_2(q);
                    float w = w_min + l * (w_max - w_min) / (w_pts - 1);
                    string key = vec_to_string(q) + " " + to_string(w);
                    matsubara_cube.cube[i][j][k][l] = map[key];
                }
            }
        }
    }

    return matsubara_cube;
}

float calculate_chi_from_cube_map(const unordered_map<float, vector<vector<vector<float>>>> &chi_cube_map, Vec q, float w) {
    Vec v = to_IBZ_2(q);
    float d = 2*k_max/(m-1);

    float x = v(0), y = v(1), z = v(2);
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
    //cout << q << q.vals(2) << endl;
    //int s = chi_cube.size()-1; 
    //assert( i < s and j < s and k < chi_cube[0][0].size()-1);

    float f1 = chi_cube_map.at(w)[i][j][k], f2 = chi_cube_map.at(w)[i+1][j][k];
    float f3 = chi_cube_map.at(w)[i+1][j+1][k], f4 = chi_cube_map.at(w)[i][j+1][k];
    float f5 = chi_cube_map.at(w)[i][j][k+1], f6 = chi_cube_map.at(w)[i+1][j][k+1];
    float f7 = chi_cube_map.at(w)[i+1][j+1][k+1], f8 = chi_cube_map.at(w)[i][j+1][k+1];

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

