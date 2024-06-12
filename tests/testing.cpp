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

//#include <boost/math/tools/roots.hpp>
//#include <Eigen/Dense>
//#include <gsl/gsl_integration.h>
//#include <python3.10/Python.h>

#include "cfg.h"
#include "fermi_surface.h"
#include "vec.h"
#include "matrix.hpp"
#include "eigenvec.hpp"
//#include "analysis.h"
#include "band_structure.h"
//#include "py_port.h"
#include "save_data.h"
#include "utilities.h"
#include "potential.h"
#include "frequency_inclusion.hpp"

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

void test_get_k() {
    int n = 5;
    double correct_values[n] = {-M_PI, -M_PI/2, 0, M_PI/2, M_PI};
    for (int i = 0; i < n; i++) {
        double value = get_k(i, n);
        if ( fabs(value - correct_values[i]) > 0.000001)
            cout << "get_k() failed test.\n";
    }
}

void test_set(vector<Vec> FS_vecs) {
    vector<string> q;
    for (int i = 0; i < FS_vecs.size(); i++) {
        for (int j = 0; j < FS_vecs.size(); j++) {
            Vec k_plus = to_IBZ_spherical(FS_vecs[i] + FS_vecs[j]);
            Vec k_minus = to_IBZ_spherical(FS_vecs[i] - FS_vecs[j]);
            q.push_back(vec_to_string(k_plus));
            q.push_back(vec_to_string(k_minus));
        }
    }
    unordered_set<string> s( q.begin(), q.end() );
    for (auto x: s)
        cout << x << endl;
    cout << q.size() << endl;
    cout << s.size() << endl;

}

void test_IBZ2() {
    int num = 10;
    for (double i = 0; i < num; i++) {
        for (double j = 0; j < num; j++) {
            for (double k = 0; k < num; k++) {
                double x = get_k(i, num);
                double y = get_k(j, num);
                double z = get_k(k, num);
                Vec input(x, y, z);
                Vec output = to_IBZ_2(input);
                double chi1 = chi_trapezoidal(input, 0.25, -1.2, 0, 40);
                double chi2 = chi_trapezoidal(output, 0.25, -1.2, 0, 40);
                if (fabs(chi1 - chi2) > 0.0000001) 
                    cout << "Failed for " << input << "->" << output << endl;
            }
        }
    }
}

void test_IBZ() {
    double T = 0.25, mu = -1.2;
    Vec q1(1, 2, .7);
    Vec q2(.1, -2, .7);
    Vec q3(-.8, -1.2, -1.1);
    double chi1 = chi_trapezoidal(q1, T, mu, 0, 30);
    double chi2 = chi_trapezoidal(q2, T, mu, 0, 30);
    double chi3 = chi_trapezoidal(q3, T, mu, 0, 30);
    Vec q1_t = to_IBZ_2(q1);
    Vec q2_t = to_IBZ_2(q2);
    Vec q3_t = to_IBZ_2(q3);
    double chi1_t = chi_trapezoidal(q1_t, T, mu, 0, 30);
    double chi2_t = chi_trapezoidal(q2_t, T, mu, 0, 30);
    double chi3_t = chi_trapezoidal(q3_t, T, mu, 0, 30);
    cout << chi1 << " " << chi1_t << endl;
    cout << chi2 << " " << chi2_t << endl;
    cout << chi3 << " " << chi3_t << endl;
}

//void test_pairing_interaction_parity(vector<Vec> &FS, double T, double mu) {
//    double DOS = 0; for (auto x : FS) DOS += x.area / vp(x);
//    DOS /= pow(2*M_PI, 3);
//    auto cube = chi_cube(T, mu, DOS, 0);
//    Vec k2;
//    double num = 10;
//    for (double i = 0; i < num; i++) {
//        for (double j = 0; j < num; j++) {
//            for (double k = 0; k < num; k++) {
//                double x = get_k(i, num);
//                double y = get_k(j, num);
//                double z = get_k(k, num);
//
//                Vec k1(x, y, z);
//
//                double V1 = potential_scalapino_cube(k1, k2, T, cube);
//                double V2 = potential_scalapino_cube(-1*k1, k2, T, cube);
//
//                if ( fabs(V1 - V2) > 0.00001) 
//                    cout << k1 << "->" << V1 << ", " << V2 << endl;
//            }
//        }
//    }
//
//}

//bool test_potential_test() { 
//    vector<Vec> FS = tetrahedron_method(mu); cout << "FS created\n";
//    double T = 0.25;
//    vector<vector<vector<double>>> cube;
//    MatrixXd P = create_P(FS, T, cube); cout << "Matrix Created\n";
//
//    EigenSolver<MatrixXd> s(P);
//    
//    VectorXcd vals = s.eigenvalues();// * f / FS.size();
//    //VectorXcd vals = s.eigenvalues() * f / FS.size();
//    EigenSolver<MatrixXd>::EigenvectorsType vecs;
//    vecs = s.eigenvectors(); 
//    
//    cout << vals << endl;
//
//    //if (vals[0] == 1.0 and vals[1] == 0.5) return true;
//    return false;
//}

//void eigenvalue_divergence() {
//    ofstream temporary_file("eigenvalue_divergence.txt");
//    for (int i = 1; i < 20; i++) {
//        printf("Plot Progress: %i out of 50\n", i);
//
//        double cutoff = 0.03 * i;
//        init_config(mu, U, t, tn, w_D, mu, U, t, tn, cutoff);
//
//        vector<vector<Vec>> freq_FS;
//        freq_FS = freq_tetrahedron_method(mu);
//        vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
//        double T = 0.25;
//
//        double DOS = get_DOS(freq_FS[(l+1)/2 - 1]);
//        unordered_map<double, vector<vector<vector<double>>>> cube;
//        if (potential_name != "test") cube = chi_cube_freq(T, mu, DOS);
//        cout << "Frequencies: " << cube.size() << endl;
//        for (auto x : cube) {
//            cout << x.first << endl;
//        }
//
//        Matrix Pf2; create_P_freq(Pf2, freq_FS, T, cube);
//        Matrix P; create_P(P, FS, T, cube);
//        double f = f_singlet_integral(T);
//
//        vector<Eigenvector> answers = power_iteration(P, 0.001);
//        vector<Eigenvector> answersf2 = power_iteration(Pf2, 0.001);
//        double eig = answers[answers.size() - 1].eigenvalue;
//        double eigf2 = answersf2[answersf2.size() - 1].eigenvalue;
//        temporary_file << w_D << " " << f*eig << " " << eigf2 << endl;
//    }
//}

void integral_convergence(double T) {
    int n = 50;
    double mag = 0.04*M_PI;
    Vec q(mag, mag, mag);
    double c2 = integrate_susceptibility(q, T, mu, 0.5, 800);
    for (int i = 0; i < 10; i++) {
        double c1 = integrate_susceptibility(q, T, mu, 0.5, n);
        cout << "Points: " << n << " , Results: " << c1 << " " << c2 << endl;
        n += 50;
    }
}

void imaginary_integral_convergence(double T, double w) {
    double mag = M_PI;
    Vec q(mag, mag, mag);
    //cout << "Answer: " << integrate_susceptibility(q, T, mu, w, 500) << endl;
    for (int i = 50; i < 500; i+=50) {
        double c1 = imaginary_integration(q, T, mu, w, i, 0.0001);
        double c2 = chi_trapezoidal(q, T, mu, w, i);
        double c3 = integrate_susceptibility(q, T, mu, w, i);
        //double c3 = 0;
        cout << "Points: " << i << ", Results: " << c1 << " " << c2 << " " << c3 << endl;
    }
}

void split_surface_vals(Vec q) {
    vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
    double min_val = 1000;
    double max_val = -1000;
    for (auto k : FS) {
        cout << k << " " << e_split(k, q) << endl;
        if (e_split(k, q) < min_val) min_val = e_split(k, q);
        if (e_split(k, q) > max_val) max_val = e_split(k, q);
        if (e_split(k, q) > max_val) max_val = e_split(k, q);
        if (e_split(k, q) > max_val) max_val = e_split(k, q);
    }
    cout << min_val << " " << max_val << endl;
}

void mu_to_n() {
    for (int i = 0; i < 5000; i++) {
        double newmu = -0.4 + 0.2/1000.0 * i;
        cout << "Mu: " << newmu << 
            " n: " << 2 * integrate_susceptibility(Vec(0,0,0), 0, newmu, 0, 1000) << endl;
    }
}

void eig_with_freq(double cutoff) {
    ofstream file("eig_freq.dat", std::ios_base::app);
    double T = 0.01;
    init_config(mu, U, t, tn, w_D, mu, U, t, tn, cutoff);
    
    vector<vector<Vec>> freq_FS = freq_tetrahedron_method(mu);
    vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
    printf("FS created\n");

    unordered_map <double, vector<vector<vector<double>>>> cube_freq_map;
    if (potential_name.find("scalapino") != string::npos)
        cube_freq_map = chi_cube_freq(T, mu, 0.0);
    printf("Cube Created\n");

    int size = 0; for (auto x : freq_FS) size += x.size();
    printf("Size: %i\n", size);

    Matrix P(size); create_P_freq(P, freq_FS, T, cube_freq_map); 
    printf("Matrix Created\n");
    Matrix P2(FS.size()); create_P(P2, FS, T, cube_freq_map);
    printf("Matrix Created\n");

    printf("Sample P values: %f %f %f %f\n", P(0,0), P(0,1), P(1,0), P(1,1));
    printf("Sample P2 values: %f %f %f %f\n", P2(0,0), P2(0,1), P2(1,0), P2(1,1));
    
    double f = f_singlet_integral(T);

    vector<Eigenvector> answers = power_iteration(P, 0.001);
    vector<Eigenvector> answers2 = power_iteration(P2, 0.001);

    double eig = answers[answers.size() - 1].eigenvalue;
    double eig2 = answers2[answers2.size() - 1].eigenvalue;

    cout << "w=0, w>0 Eigs: " << f*eig2 << " " << eig << endl;
    file << cutoff << " " << f*eig2 << " " << eig << endl;
}

int main() {
    
    int num_procs = omp_get_num_procs();
    omp_set_num_threads(num_procs - 1);

    double T = 0.25, w = 0.0;
    //plot_single_chi(T, w);
    //plot_single_chi2(T, w);
    //plot_single_chi3(T, w);
    //return 0;

    vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
    unordered_map <double, vector<vector<vector<double>>>> cube_freq_map;
    Matrix P(FS.size()); 
    printf("P size: %d, val: %f\n", P.size, P(0,0));
    save_potential_vs_q(FS, P, "potential.dat");
    printf("P size: %d, val: %f\n", P.size, P(0,0));
    create_P(P, FS, T, cube_freq_map);
    printf("P size: %d, val: %f\n", P.size, P(0,0));
    vector<Eigenvector> answers = power_iteration(P, 0.001);
    printf("P size: %d, val: %f\n", P.size, P(0,0));

    for (int i = 0; i < answers.size(); i++) {
        printf("Eigenvalue: %f\n", answers[i].eigenvalue);
    }

    Eigenvector y(FS.size());
    for (int i = 0; i < y.size; i++) {
        y.eigenvector[i] = 1;
    }
    printf("Eigenvalue: %f\n", dot(y, P * y) / dot(y, y));

    Eigenvector r(FS.size(), true);
    printf("Random Eigenvector:\n");
    for (int i = 0; i < r.size; i++) {
        cout << r.eigenvector[i] << " ";
    }
    cout << endl;

    Eigenvector after = P * r;
    printf("Multiplied Eigenvector:\n");
    for (int i = 0; i < after.size; i++) {
        cout << after.eigenvector[i] << " ";
    }
    after.normalize();
    cout << endl;

    vector<double> x2(FS.size(), 0);
        srand(time(0)); 
    for (int i = 0; i < x2.size(); i++) {
        x2[i] = ((double)rand()) / RAND_MAX;
    }
    printf("Random Eigenvector:\n");
    for (int i = 0; i < x2.size(); i++) {
        cout << x2[i] << " ";
    }
    //Eigenvector x(FS.size());
    //for (int i = 0; i < x.size; i++) {
    //    x.eigenvector[i] = 1;
    //}
    //printf("\nFirst Eigenvector: %f %f %f %f\n", x.eigenvector[0], x.eigenvector[1], x.eigenvector[2], x.eigenvector[3]);
    //printf("Multiplied Eigenvector:\n");
    //cout << P * x << endl;

    //double result = 0;
    //for (int i = 0; i < x.size; i++) {
    //    result += x.eigenvector[i] * P(i, 0);
    //}
    //printf("Result: %f\n", result);





    return 0;



    for (int i = 1; i < 30; i++) {
        double cutoff = 0.03 * i;
        cout << "Cutoff: " << cutoff << endl;
        init_config(mu, U, t, tn, w_D, mu, U, t, tn, cutoff);
        eig_with_freq(cutoff);
    }
    return 0;
}
