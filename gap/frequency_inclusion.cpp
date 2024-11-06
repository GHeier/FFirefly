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

#include "calculations.h"
#include "cfg.h"
#include "vec.h"
#include "utilities.h"
#include "frequency_inclusion.hpp"
#include "susceptibility.h"
#include "band_structure.h"

using std::cout;
using std::endl;
using std::sort;
using std::vector;
using std::unordered_map;
//using lambda_lanczos::LambdaLanczos;

// Gets total matrix size
int matrix_size_from_freq_FS(vector<vector<Vec>> &freq_FS) {
    int size = 0;
    for (int i = 0; i < freq_FS.size(); i++) {
        size += freq_FS[i].size();
    }
    return size;
}

// Defines all energy surfaces around FS
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

// Creates the P matrix based around the multiple energy surfaces calculated above
void create_P_freq(Matrix &P, vector<vector<Vec>> &k, float T, const unordered_map<float, vector<vector<vector<float>>>> &chi_cube2) {

    cout << "Creating P Matrix with frequency\n";
    for (int i = 0; i < k.size(); i++) {

        int ind1 = 0;
        for (int temp = 0; temp < i; temp++)
            ind1 += k[temp].size();

        for (int j = 0; j < k[i].size(); j++) {
            Vec k1 = k[i][j];
            for (int x = 0; x < k.size(); x++) {

                int ind2 = 0;
                for (int temp = 0; temp < x; temp++)
                    ind2 += k[temp].size();

                #pragma omp parallel for
                for (int y = 0; y < k[x].size(); y++) {
                    Vec k2 = k[x][y];
                    float d1 = pow(k1.area/vp(k1),0.5); 
                    float d2 = pow(k2.area/vp(k2),0.5); 
                    // f * d_epsilon
                    float fde1 = f_singlet(wc * points[l-1][i], T) * weights[l-1][i];
                    float fde2 = f_singlet(wc * points[l-1][x], T) * weights[l-1][x];
                    float w = wc * (points[l-1][x] - points[l-1][i]);

                    P(ind1 + j,ind2 + y) = (float)(- d1 * d2 * pow(fde1*fde2,0.5) * V(k1, k2, w, T, chi_cube2)); 
                }
            }
            string message = "Portion " + to_string(i) + " of " + to_string(k.size());
            progress_bar(1.0 * (ind1 + j) / P.size, message);
        }
    }
    cout << "P Matrix Created\n";

    //return P * 2 * wc / (l * k_size);
    P *= wc * (2 / pow(2*M_PI, dim)); 
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
            w = round(w, 6);
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

