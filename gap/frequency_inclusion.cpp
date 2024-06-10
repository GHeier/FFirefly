#include <iostream>
#include <unistd.h>
#include <ctime>

#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <boost/functional/hash.hpp>

//#include <lambda_lanczos/lambda_lanczos.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/gauss.hpp>

#include "calculations.h"
#include "band_structure.h"
#include "cfg.h"
#include "fermi_surface.h"
#include "potential.h"
#include "vec.h"
#include "utilities.h"
#include "frequency_inclusion.hpp"

using std::cout;
using std::endl;
using std::sort;
using std::vector;
using std::unordered_map;
//using lambda_lanczos::LambdaLanczos;

vector<vector<Vec>> freq_tetrahedron_method(double mu) {
    assert( l % 2 != 0); // N must be odd that way frequencies are evenly spaced
    vector<vector<Vec>> basis;

    double points_0th[1] = {0}; double *p0 = points_0th;
    double points_1st[2] = {-1/pow(3,0.5), 1/pow(3,0.5)}; double *p1 = points_1st;
    double points_2nd[3] = {-pow(3/5,0.5), 0, pow(3/5,0.5)}; double *p2 = points_2nd;
    double points_3rd[4] = {-0.861136, -0.339981, 0.339981, 0.861136}; double *p3 = points_3rd;
    double points_4th[5] = {-0.90618, -0.538469, 0, 0.538469, 0.90618}; double *p4 = points_4th;

    double *points[5] = {p0, p1, p2, p3, p4};

    for (int i = 0; i < l; i++) {
        double ep = w_D * points[l-1][i] + mu;
        vector<Vec> layer = tetrahedron_method(e_base_avg, Vec(0,0,0), ep);
        basis.push_back(layer);
    }
    return basis;
}

void create_P_freq(Matrix &P, vector<vector<Vec>> &k, double T, const unordered_map<double, vector<vector<vector<double>>>> &chi_cube2) {

    int size = 0;
    for (auto x : k) size += x.size();
    cout << "Creating P Matrix with frequency\n";
    #pragma omp parallel for
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

                for (int y = 0; y < k[x].size(); y++) {
                    Vec k2 = k[x][y];
                    double d1 = pow(k1.area/vp(k1),0.5); 
                    double d2 = pow(k2.area/vp(k2),0.5); 
                    // f * d_epsilon
                    double fde1 = f_singlet(w_D * points[l-1][i], T) * weights[l-1][i];
                    double fde2 = f_singlet(w_D * points[l-1][x], T) * weights[l-1][x];
                    double w = w_D * (points[l-1][x] - points[l-1][i]);

                    P(ind1 + j,ind2 + y) = - d1 * d2 * pow(fde1*fde2,0.5) * V(k1, k2, w, T, chi_cube2); 
                }
            }
        }
    }
    cout << "P Matrix Created\n";

    //return P * 2 * w_D / (l * k_size);
    P = P * w_D * (2 / pow(2*M_PI, dim)); 
}

unordered_map <double, vector<vector<vector<double>>>> chi_cube_freq(double T, double mu, double DOS) {
    vector<double> des;
    //cout << "Desired frequencies\n";
    for (int i = 0; i < l; i++) {
        double p1 = w_D * points[l-1][i];
        for (int j = 0; j < l; j++) {
            double p2 = w_D * points[l-1][j];
            double w = p1 - p2;
            w = round(w, 6);
            if (count(des.begin(), des.end(), w) == 0)
                des.push_back(w);
        }
    }

//    cout << "Taken frequencies\n";
//    for (double w : des) {
//        cout << "?" << w << endl;
//    }

    unordered_map <double, vector<vector<vector<double>>>> cube_freq_map;
    for (int i = 0; i < des.size(); i++) {
        printf("Chi Cube %d / %d\n", i+1, des.size());
        double w = des[i];
        auto cube = chi_cube(T, mu, DOS, w);
        cube_freq_map.insert(pair<double, vector<vector<vector<double>>>>(w, cube));
    }
//    for (auto x : cube_freq_map) {
//        cout << "Frequency: " << x.first << endl;
//    }
    return cube_freq_map;
}

double calculate_chi_from_cube_map(const unordered_map<double, vector<vector<vector<double>>>> &chi_cube_map, Vec q, double w) {
    Vec v = to_IBZ_2(q);
    double d = 2*k_max/(m-1);

    double x = v.vals[0], y = v.vals[1], z = v.vals[2];
    if (dim == 2) z = 0;

    int i = floor(x / d);
    int j = floor(y / d);
    int k = floor(z / d);

    double x1 = i * d; 
    double y1 = j * d; 
    double z1 = k * d; 

    double x2 = x1 + d; 
    double y2 = y1 + d; 
    double z2 = z1 + d; 

    double dx = 0, dy = 0, dz = 0, wx = 0, wy = 0, wz = 0, w0 = 0;

    // Make sure there's no issue with indexing
    //cout << q << q.vals(2) << endl;
    //int s = chi_cube.size()-1; 
    //assert( i < s and j < s and k < chi_cube[0][0].size()-1);

    double f1 = chi_cube_map.at(w)[i][j][k], f2 = chi_cube_map.at(w)[i+1][j][k];
    double f3 = chi_cube_map.at(w)[i+1][j+1][k], f4 = chi_cube_map.at(w)[i][j+1][k];
    double f5 = chi_cube_map.at(w)[i][j][k+1], f6 = chi_cube_map.at(w)[i+1][j][k+1];
    double f7 = chi_cube_map.at(w)[i+1][j+1][k+1], f8 = chi_cube_map.at(w)[i][j+1][k+1];

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

void get_bounds3(Vec q, double &upper, double &lower, double (*func)(Vec k, Vec q)) {
    upper = 0; lower = 1000;
    int pts = 100;
    for (double i = 0; i < pts; i++) {
        double x = get_k(i, pts);
        for (double j = 0; j < pts; j++) {
            double y = get_k(j, pts);
            for (double k = 0; k < pts * (dim%2) + 1 * ((dim+1)%2); k++) {
                double z = get_k(k, pts);
                Vec k_val(x, y, z);
                double val = func(k_val, q);
                if (val > upper) upper = val;
                if (val < lower) lower = val;
            }
        }
    }
    lower *= 0.99; upper *= 0.99;
}

double denominator(Vec k, Vec q) {
    return e_diff(k, q);
}

double denominator_diff(Vec k, Vec q) {
    return vp_diff(k, q);
}

double integrand(Vec k, Vec q, double w, double T) {
    //return 1;
    //if (fabs(k.freq - w) < 0.001) return 0;
    //return 1 / (k.freq - w);

    double e_qk = epsilon(k+q) - mu;
    double e_k = epsilon(k) - mu;
    double f_k = f(e_k, T);
    double f_qk = f(e_qk, T);

    double dE = k.freq;
    //if (dE != 0) cout << dE << endl;
    //if (fabs(e_qk - e_k) > 0.1) {
    //    cout << "dE: " << dE << " " << e_qk - e_k << endl;
    //}

    if (fabs(dE) < 0.0001 and fabs(w) < 0.0001) {
        if (T == 0 or exp(e_k/T) > 1e6)
            return e_k < 0;
        return 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
    }
    if (fabs(w - dE) < 0.001) return 0;
    return (f_qk - f_k) / (w - dE);
}

void get_spacing_curve_consts(double w, double a, double b, double &A, double &upr, double &lwr) {
    A = b - w;
    lwr = (a - w) / A;
    upr = 1;
    if ( w - a > b - w) {
        A = w - a;
        lwr = -1;
        upr = (b - w) / A;
    }
    //A = (w - a);
    //upr = pow((b - w) / (w - a), 1.0/3.0);
    //if (isnan(upr)) upr = - pow(fabs(b-w)/fabs(w-a), 1.0/3.0);
    //lwr = -1;

    //if (w - a > b - w) {
    //    A = (b - w);
    //    lwr = - pow((w - a) / (b - w), 1.0/3.0);
    //    if (isnan(lwr)) lwr = - pow(fabs(b-w)/fabs(w-a), 1.0/3.0);
    //    upr = 1;
    //}
}

void get_spacing_vec(vector<double> &spacing, double w, double a, double b, int pts) {
    double A, upr, lwr;
    get_spacing_curve_consts(w, a, b, A, upr, lwr);

    auto spacing_curve = [A, w] (double i, double pts) { 
        double x = -1 + 2 * i / pts;
        return A * x + w;
        return A * pow(x,3) + w; 
    };

    double r = spacing_curve(0, pts);
    for (int i = 0; r < b; i++) {
        double t = i;
        r = spacing_curve(i, pts);
        spacing.push_back(r);
    }
}
