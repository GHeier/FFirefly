#include <iomanip>
#include <iostream>
#include <fstream>
//#include <chrono>
//#include <bits/stdc++.h>

#include <cstring>
#include <string>
//#include <math.h>
#include <vector>

#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <gsl/gsl_integration.h>
#include <omp.h>
#include <functional>

//#include <boost/math/tools/roots.hpp>
//#include <Eigen/Dense>
//#include <gsl/gsl_integration.h>
//#include <python3.10/Python.h>

#include "analysis.h"
#include "band_structure.h"
#include "calculations.h"
#include "cfg.h"
#include "eigenvec.hpp"
#include "fermi_surface.h"
#include "frequency_inclusion.hpp"
#include "matrix.hpp"
#include "potential.h"
#include "save_data.h"
#include "susceptibility.h"
#include "utilities.h"
#include "vec.h"

#include "plot.cpp"
//#include "surface_integrals.cpp"
//#include "freq_test.cpp"


using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
//using namespace Eigen;


bool test_coordinate_transforms() {
    Vec x(1.0, 0.0, 0.0, false);
    Vec y(1.5, 1.0, 0.0, false);
    Vec z(2.0, 0.0, 0.0, false);
    Vec a(1.0, 0.0, M_PI/2, false);
    Vec b(1.0, 1.4, -0.5, false);
    Vec c(6.06, 0.34, -0.96, false);

    Vec x1 = x; Vec y1 = y; Vec z1 = z; Vec a1 = a; Vec b1 = b; Vec c1 = c;

    x.to_cartesian();
    y.to_cartesian();
    z.to_cartesian();
    a.to_cartesian();
    b.to_cartesian();
    c.to_cartesian();

    if (fabs(x.vals[0] - 1.0) > 0.00001 or 
            fabs(x.vals[1] - 0.0) > 0.00001)
        cout << "Transform 1 Failed\n";
    if (fabs(y.vals[0] - 0.810453) > 0.00001 or 
            fabs(y.vals[1] - 1.26221) > 0.00001) 
        cout << "Transform 2 Failed\n";
    if (fabs(z.vals[0] - 2) > 0.00001 or 
            fabs(z.vals[1] - 0.0) > 0.00001 or 
            fabs(z.vals[2] - 0.0) > 0.00001) 
        cout << "Transform 3 Failed\n";
    if (fabs(a.vals[0] - 0) > 0.00001 or 
            fabs(a.vals[1] - 0.0) > 0.00001 or 
            fabs(a.vals[2] - 1.0) > 0.00001) 
        cout << "Transform 4 Failed\n";
    if (fabs(b.vals[0] - 0.149160) > 0.00001 or 
            fabs(b.vals[1] - 0.864813) > 0.00001 or 
            fabs(b.vals[2] - (-0.479426)) > 0.00001) 
        cout << "Transform 5 Failed\n";
    if (fabs(c.vals[0] - 3.27657) > 0.0001 or 
            fabs(c.vals[1] - 1.15904) > 0.0001 or 
            fabs(c.vals[2] - (-4.9643)) > 0.0001) 
        cout << "Transform 6 Failed\n";

    x.to_spherical();
    y.to_spherical();
    z.to_spherical();
    a.to_spherical();
    b.to_spherical();
    c.to_spherical();

    for (int i = 0; i < dim; i++) if (fabs(x.vals[i] - x1.vals[i]) > 0.00001) 
        cout << "Cartesian to spherical Transform #1 Failed";
    for (int i = 0; i < dim; i++) if (fabs(y.vals[i] - y1.vals[i]) > 0.00001) 
        cout << "Cartesian to spherical Transform #2 Failed";
    for (int i = 0; i < dim; i++) if (fabs(z.vals[i] - z1.vals[i]) > 0.00001) 
        cout << "Cartesian to spherical Transform #3 Failed";
    for (int i = 0; i < dim; i++) if (fabs(a.vals[i] - a1.vals[i]) > 0.00001) 
        cout <<"Cartesian to spherical Transform #4 Failed";
    for (int i = 0; i < dim; i++) if (fabs(b.vals[i] - b1.vals[i]) > 0.00001) 
        cout << "Cartesian to spherical Transform #5 Failed";
    for (int i = 0; i < dim; i++) if (fabs(c.vals[i] - c1.vals[i]) > 0.00001) 
        cout << "Cartesian to spherical Transform #6 Failed";
    cout << "All Transforms Completed. If no error messages, the test is passed\n";

    return true;
}

void test_IBZ2() {
    float T = 0.25;
    float w = 0.0;
    int num = 10;
    for (float i = 0; i < num; i++) {
        for (float j = 0; j < num; j++) {
            for (float k = 0; k < num; k++) {
                float x = i/num * 2 * k_max;
                float y = j/num * 2 * k_max;
                float z = k/num * 2 * k_max;
                Vec q(x, y, z); Vec input = q;
                auto f = [T, w, q] (float x, float y, float z) {
                    Vec k(x, y, z);
                    return integrand(k, q, w, T);
                };
                float chi1 = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60);
                q = to_IBZ_2(q); Vec output = q;
                float chi2 = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60);
                if (fabs(chi1 - chi2) > 0.0000001) 
                    cout << "Failed for " << input << "->" << output << endl;
            }
        }
    }
}

void test_IBZ() {
    float T = 0.25, mu = -1.2, w = 0.0;
    Vec q(1, 2, .7);
    auto f = [T, w, q] (float x, float y, float z) {
        Vec k(x, y, z);
        return integrand(k, q, w, T);
    };
    float chi1 = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60);
    q = to_IBZ_2(q);
    float chi1_t = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60);
    cout << chi1 << " " << chi1_t << endl;

    q = Vec(.1, -2, .7);
    float chi2 = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60);
    q = to_IBZ_2(q);
    float chi2_t = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60);
    cout << chi2 << " " << chi2_t << endl;

    q = Vec(-.8, -1.2, -1.1);
    float chi3 = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60);
    q = to_IBZ_2(q);
    float chi3_t = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60);
    cout << chi3 << " " << chi3_t << endl;
}


//void eigenvalue_divergence() {
//    ofstream temporary_file("eigenvalue_divergence.txt");
//    for (int i = 1; i < 20; i++) {
//        printf("Plot Progress: %i out of 50\n", i);
//
//        float cutoff = 0.03 * i;
//        init_config(mu, U, t, tn, wc, mu, U, t, tn, cutoff);
//
//        vector<vector<Vec>> freq_FS;
//        freq_FS = freq_tetrahedron_method(mu);
//        vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
//        float T = 0.25;
//
//        float DOS = get_DOS(freq_FS[(l+1)/2 - 1]);
//        unordered_map<float, vector<vector<vector<float>>>> cube;
//        if (potential_name != "test") cube = chi_cube_freq(T, mu, DOS);
//        cout << "Frequencies: " << cube.size() << endl;
//        for (auto x : cube) {
//            cout << x.first << endl;
//        }
//
//        Matrix Pf2; create_P_freq(Pf2, freq_FS, T, cube);
//        Matrix P; create_P(P, FS, T, cube);
//        float f = f_singlet_integral(T);
//
//        vector<Eigenvector> answers = power_iteration(P, 0.001);
//        vector<Eigenvector> answersf2 = power_iteration(Pf2, 0.001);
//        float eig = answers[answers.size() - 1].eigenvalue;
//        float eigf2 = answersf2[answersf2.size() - 1].eigenvalue;
//        temporary_file << wc << " " << f*eig << " " << eigf2 << endl;
//    }
//}

void integral_convergence(float T) {
    int n = 50;
    float mag = 0.04*M_PI;
    Vec q(mag, mag, mag);
    float c2 = integrate_susceptibility(q, T, mu, 0.5, 800);
    for (int i = 0; i < 10; i++) {
        float c1 = integrate_susceptibility(q, T, mu, 0.5, n);
        cout << "Points: " << n << " , Results: " << c1 << " " << c2 << endl;
        n += 50;
    }
}

void mu_to_n() {
    for (int i = 0; i < 5000; i++) {
        float newmu = -0.4 + 0.2/1000.0 * i;
        cout << "Mu: " << newmu << 
            " n: " << 2 * integrate_susceptibility(Vec(0,0,0), 0, newmu, 0, 1000) << endl;
    }
}

void eig_with_freq(float cutoff) {
    ofstream file("eig_freq.dat", std::ios_base::app);
    float T = 0.01;
    init_config(mu, U, t, tn, wc, mu, U, t, tn, cutoff);
    
    vector<vector<Vec>> freq_FS = freq_tetrahedron_method(mu);
    vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
    printf("FS created\n");

    unordered_map <float, vector<vector<vector<float>>>> cube_freq_map;
    if (potential_name.find("scalapino") != string::npos)
        cube_freq_map = chi_cube_freq(T, mu);
    printf("Cube Created\n");

    int size = 0; for (auto x : freq_FS) size += x.size();
    printf("Size: %i\n", size);

    Matrix P(size); create_P_freq(P, freq_FS, T, cube_freq_map); 
    printf("Matrix Created\n");
    Matrix P2(FS.size()); create_P(P2, FS, T, cube_freq_map);
    printf("Matrix Created\n");

    printf("Sample P values: %f %f %f %f\n", P(0,0), P(0,1), P(1,0), P(1,1));
    printf("Sample P2 values: %f %f %f %f\n", P2(0,0), P2(0,1), P2(1,0), P2(1,1));
    
    float f = f_singlet_integral(T);

    vector<Eigenvector> answers = power_iteration(P, 0.001);
    vector<Eigenvector> answers2 = power_iteration(P2, 0.001);

    float eig = answers[answers.size() - 1].eigenvalue;
    float eig2 = answers2[answers2.size() - 1].eigenvalue;

    cout << "w=0, w>0 Eigs: " << f*eig2 << " " << eig << endl;
    file << cutoff << " " << f*eig2 << " " << eig << endl;
}

void test_cube_map() {
    unordered_map<float, vector<vector<vector<float>>> > cube_freq_map;
    cube_freq_map = chi_cube_freq(0.25, -1.2);
    auto cube = chi_cube(0.25, -1.2, 0.0, "Cube 1/1");
    for (int i = 0; i < cube.size(); i++) {
        for (int j = 0; j < cube[i].size(); j++) {
            for (int k = 0; k < cube[i][j].size(); k++) {
                if (fabs(cube[i][j][k] - cube_freq_map.at(0.0)[i][j][k]) > 0.0001) {
                    printf("Cube values: %f %f\n", cube[i][j][k], cube_freq_map.at(0.0)[i][j][k]);
                }
                if (i == 0) cout << cube[i][j][k] << " " << cube_freq_map.at(0.0)[i][j][k] << endl;
            }
        }
    }
}

int main() {
    
    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs - 1);

    //test_cube_map();

    float T = 0.25, w = 0.0;
    printf("Starting Test\n");
    plot_single_chi(T, w);
    printf("Test Complete\n");
    plot_single_chi2(T, w);
    printf("Test Complete\n");
    plot_single_chi3(T, w);
    printf("Test Complete\n");
    return 0;
}
