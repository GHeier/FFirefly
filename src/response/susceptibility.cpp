/**
 * @file susceptibility.cpp
 *
 * @brief Calculates the susceptibility of the system.
 *
 * @author Griffin Heier
 */

#include <iostream>
#include <math.h>
#include <complex>
#include <string>
#include <algorithm>
#include <functional>

#include <omp.h>
#include <boost/functional/hash.hpp>
#include <unordered_map>

#include "../superconductor/utilities.h"
#include "../objects/vec.h"
#include "../cfg.h"
#include "susceptibility.h"
#include "../algorithms/integration.h"
#include "../band_structure.h"

using namespace std;

//Fermi-Dirac distribution function
float fermi_dirac(float E, float T) {
    if (T == 0) {
        if (E < 0) return 1;
        return 0;
    }
    return 1 / (1 + exp(E/T));
}

float ratio(Vec k, Vec q, float w, float T) {
    float e_k = epsilon(k) - mu;
    float e_qk = epsilon(k+q) - mu;
    float dE = e_qk - e_k;
    float f_k = fermi_dirac(e_k, T);
    float f_qk = fermi_dirac(e_qk, T);

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
        vector<Vec> FS = get_FS(mu);
        float sum = 0; for (auto x : FS) sum += x.area / vp(x);
        return sum / pow(2*k_max,dim);
    }
    if (q.norm() < 0.0001) q = Vec(0.01,0.01,0.01);
    auto func = [q, w, T](Vec k) -> float {
        return ratio(k, q, w, T);
    };
    auto e_diff = [q](Vec k) -> float {
        return epsilon(k+q) - epsilon(k);
    };
    auto denom = [q, w](Vec k) -> float { 
        return w - (epsilon(k+q) - epsilon(k)); 
    };
    auto denom_diff = [q, w](Vec k) -> float { 
        return vp_diff(k, q); 
    };

    float a, b; get_surface_transformed_bounds(b, a, e_diff);
    vector<float> spacing; get_spacing_vec(spacing, w, a, b, num_points);
    //float c = tetrahedron_sum_continuous(denominator, denominator_diff, q, spacing, w, T);
    float c = surface_transform_integral(func, denom, denom_diff, spacing);
    return c / pow(2*k_max,dim);
}

complex<float> complex_susceptibility_integration(Vec q, float T, float mu, complex<float> w, int num_points) {
    auto func = [T, q, w, mu](float x, float y, float z) -> complex<float> {
        Vec k(x,y,z);
        float e_k = epsilon(k) - mu;
        float e_kq = epsilon(k+q) - mu;
        float f_kq = fermi_dirac(e_kq, T);
        float f_k = fermi_dirac(e_k, T);
        if (fabs(e_kq - e_k) < 0.0001 and fabs(w.imag()) < 0.0001 and fabs(w.real()) < 0.0001) {
            if (T == 0 or exp(e_k/T) > 1e6) return e_k < 0;
            return 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        }
        return (f_kq - f_k) / (w - (e_kq - e_k));
    };
    auto func_r = [T, q, w, mu](Vec k) -> float {
        float e_k = epsilon(k) - mu;
        float e_kq = epsilon(k+q) - mu;
        float f_kq = fermi_dirac(e_kq, T);
        float f_k = fermi_dirac(e_k, T);
        if (fabs(e_kq - e_k) < 0.0001 and fabs(w.imag()) < 0.0001 and fabs(w.real()) < 0.0001) {
            if (T == 0 or exp(e_k/T) > 1e6) return e_k < 0;
            return 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        }
        return (f_kq - f_k)*(w.real() - (e_kq - e_k)) / 
            ((w.real() - (e_kq - e_k))*(w.real() - (e_kq - e_k)) + w.imag()*w.imag());
    };
    auto func_i = [T, q, w, mu](Vec k) -> float {
        float e_k = epsilon(k) - mu;
        float e_kq = epsilon(k+q) - mu;
        float f_kq = fermi_dirac(e_kq, T);
        float f_k = fermi_dirac(e_k, T);
        if (fabs(e_kq - e_k) < 0.0001 and fabs(w.imag()) < 0.0001 and fabs(w.real()) < 0.0001) {
            if (T == 0 or exp(e_k/T) > 1e6) return e_k < 0;
            return 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
        }
        return -(f_kq - f_k) * w.imag() /
            ((w.real() - (e_kq - e_k))*(w.real() - (e_kq - e_k)) + w.imag()*w.imag());
    };
    auto denom = [q, w](Vec k) -> float { 
        return w.real() - (epsilon(k+q) - epsilon(k)); 
    };
    auto denom_diff = [q, w](Vec k) -> float { 
        return vp_diff(k, q); 
    };
    float d = pow(2*k_max,dim);
    float a, b; get_surface_transformed_bounds(b, a, denom);
    vector<float> spacing; get_spacing_vec(spacing, w.real(), a, b, num_points);
    float c_real = surface_transform_integral(func_r, denom, denom_diff, spacing);
    float c_imag = surface_transform_integral(func_i, denom, denom_diff, spacing);
    complex<float> c(c_real, c_imag);
    return c / d;

    return complex_trapezoidal_integration(func, -k_max, k_max, -k_max, k_max, 
            -k_max, k_max, num_points) / d;
}

vector<vector<vector<float>>> chi_cube(float T, float mu, float w, string message) {
    int num_points = (dim == 3) ? 50 : 1000; // Number of integral surfaces
    int m_z = m*(dim%2) + 3*((dim+1)%2);
    vector<vector<vector<float>>> cube(m, vector<vector<float>> (m, vector<float> (m_z)));
    unordered_map<string, float> map;
    float empty_val = -98214214.0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < m_z; k++) {
                Vec q((2.0*k_max*i)/(1.0*m-1), (2.0*k_max*j)/(1.0*m-1), (2.0*k_max*k)/(1.0*m_z-1));
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
        map[key] = integrate_susceptibility(string_to_vec(key), T, mu, w, num_points);
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
    //cout << q << q(2) << endl;
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
    float x = q(0), y = q(1), z = q(2);
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

