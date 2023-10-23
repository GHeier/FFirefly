#include <iomanip>
#include <iostream>
#include <fstream>
//#include <chrono>

#include <cstring>
#include <string>
#include <math.h>
#include <vector>

#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <gsl/gsl_integration.h>

//#include <boost/math/tools/roots.hpp>
//#include <Eigen/Dense>
//#include <gsl/gsl_integration.h>
//#include <python3.10/Python.h>

#include "cfg.h"
#include "fermi_surface.h"
#include "vec.h"
#include "analysis.h"
#include "band_structure.h"
//#include "py_port.h"
#include "save_data.h"
#include "utilities.h"
#include "potential.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
//using namespace Eigen;

//void test_chi_interpolate() {
//    Vec q(1.0,-2.2,2.7);
//    double mu = -1.0;
//    double T = 0.25;
//    auto cube = chi_cube(T, mu);
//    double chi = calculate_chi_from_cube(cube, q);
//    double chi2 = chi_trapezoidal(q, T, mu, 40);
//    Vec q2 = to_IBZ_2(q);
//    double chi3 = chi_trapezoidal(q2, T, mu, 40);
//   cout << "Chi: " << chi << " Chi 2: " << chi2 << " Chi 3: " << chi3 << endl;
//}

//void chi_interpolation_average_error() {
//    double T = 0.25;
//    double mu = 0.0;
//    srand( (unsigned)time( NULL ) );
//    vector<vector<vector<double>>> cube = chi_cube(T, mu);
//
//    double max_error = 0;
//    double min_error = 10;
//    double error_ave = 0;
//    int n = 10000;
//    for (int i = 0; i < n; i++) {
//        double x = (double) rand()/RAND_MAX * 0.5 + 2.7;
//        double y = (double) rand()/RAND_MAX * 0.5 + 2.7;
//        double z = (double) rand()/RAND_MAX * 0.5 + 2.7;
//        Vec q(x, y, z);
//        double c = calculate_chi_from_cube(cube, q);
//        double c_exact = chi_trapezoidal(q, T, mu, 40);
//        double error = abs(c - c_exact);
//        if (error > max_error) max_error = error;
//        if (error < min_error) min_error = error;
//        error_ave += error / n;
//    }
//    std::cout << std::setprecision(8);
//    cout << "Max Error: " << max_error << endl
//        << "Min Error: " << min_error << endl
//        << "Average Error: " << error_ave << endl;
//}

void test_chi_trap() {
    Vec q(3, 4.92492, 0.742601, 0.539429, false);
    double T = 1.503114;
    double mu = 0.0;
    double num_points = 20;
    cout << chi_trapezoidal(q, T, mu, num_points) << endl;
    //cout << chi(q, T, mu, num_points) << endl;
    //cout << trapezoidal_integration_3d(num_points, &test_function) << endl;
}

long double ratio_test(Vec q, Vec k, double T, double mu, double offset) {
    long double e_qk = epsilon(q+k) - mu;
    long double e_k = epsilon(k) - mu;
    if (fabs(e_qk - e_k) < offset/10) e_qk += offset;
    long double f_k = f(e_k, T);
    long double f_qk = f(e_qk, T);
    return (f_qk - f_k) / (e_k - e_qk);
}

void test_chi_limit() {
    Vec q(3, 0.0, 0.0, 0.0); Vec k(3, 1.0, 1.8, 0.7);
    long double T = 0.25, mu = -1.0;
    long double t = 1.0, b = 1/T, e = epsilon(k);
    k.to_cartesian();
    long double limit = b * exp(b*(e - mu)) / pow( exp(b*(e - mu)) + 1,2);
    for (long double i = 1; i < 15; i++) {
        long double offset = pow(10, -i);
        cout << ratio_test(q, k, T, mu, offset) << " " << limit << endl;
    }
}
 
void test_chi_q_error() {
    long double T = 0.25, mu = 0.0;
    double q_mag = M_PI - 0.001;
    double max = 0;
    while ( q_mag > 0 ) {
        Vec k1(3, q_mag, q_mag, q_mag, true); 
        Vec k2(3, q_mag+0.001, q_mag+0.001, q_mag+0.001, true);
        q_mag -= 0.001;
        double chi1 = chi_trapezoidal(k1, T, mu, 80);
        double chi2 = chi_trapezoidal(k2, T, mu, 80);
        if (max < fabs(chi2 - chi1)) max = fabs(chi2 - chi1);
    }
    cout << max << endl;
}

void test_chi_q_limits(Vec q) {
    double T = 0.25, mu = 0.0;
    double c = chi_trapezoidal(q, T, mu, 4);
}
