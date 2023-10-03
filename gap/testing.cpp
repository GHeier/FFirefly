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

    if (fabs(x.vals(0) - 1.0) > 0.00001 or 
            fabs(x.vals(1) - 0.0) > 0.00001)
        cout << "Transform 1 Failed\n";
    if (fabs(y.vals(0) - 0.810453) > 0.00001 or 
            fabs(y.vals(1) - 1.26221) > 0.00001) 
        cout << "Transform 2 Failed\n";
    if (fabs(z.vals(0) - 2) > 0.00001 or 
            fabs(z.vals(1) - 0.0) > 0.00001 or 
            fabs(z.vals(2) - 0.0) > 0.00001) 
        cout << "Transform 3 Failed\n";
    if (fabs(a.vals(0) - 0) > 0.00001 or 
            fabs(a.vals(1) - 0.0) > 0.00001 or 
            fabs(a.vals(2) - 1.0) > 0.00001) 
        cout << "Transform 4 Failed\n";
    if (fabs(b.vals(0) - 0.149160) > 0.00001 or 
            fabs(b.vals(1) - 0.864813) > 0.00001 or 
            fabs(b.vals(2) - (-0.479426)) > 0.00001) 
        cout << "Transform 5 Failed\n";
    if (fabs(c.vals(0) - 3.27657) > 0.0001 or 
            fabs(c.vals(1) - 1.15904) > 0.0001 or 
            fabs(c.vals(2) - (-4.9643)) > 0.0001) 
        cout << "Transform 6 Failed\n";

    x.to_spherical();
    y.to_spherical();
    z.to_spherical();
    a.to_spherical();
    b.to_spherical();
    c.to_spherical();

    for (int i = 0; i < dim; i++) if (fabs(x.vals(i) - x1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #1 Failed";
    for (int i = 0; i < dim; i++) if (fabs(y.vals(i) - y1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #2 Failed";
    for (int i = 0; i < dim; i++) if (fabs(z.vals(i) - z1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #3 Failed";
    for (int i = 0; i < dim; i++) if (fabs(a.vals(i) - a1.vals(i)) > 0.00001) 
        cout <<"Cartesian to spherical Transform #4 Failed";
    for (int i = 0; i < dim; i++) if (fabs(b.vals(i) - b1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #5 Failed";
    for (int i = 0; i < dim; i++) if (fabs(c.vals(i) - c1.vals(i)) > 0.00001) 
        cout << "Cartesian to spherical Transform #6 Failed";
    cout << "All Transforms Completed. If no error messages, the test is passed\n";

    return true;
}

void plot_chi(double T) {
    ofstream file("chi_plot.dat");
    file << "q chi " << endl;
    double n = 80.0;
    for (double mu = 0.0; mu > -4.0; mu--) {
        for (double i = 0; i < n; i++) {
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            double c = chi_trapezoidal(q, T, mu, 60);
            //double c2 = varied_chi_trapezoidal(q, T, mu, 40);
            file << q_mag << " " << c << endl;
        }
        file << endl;
    }
}

void plot_potential(double T) {
    ofstream file("info.log");
    file << "q V " << endl;
    double n = 80.0;
    for (double mu = 0.0; mu > -2.0; mu--) {
        vector<vector<vector<double>>> cube = chi_cube(T, mu);
        for (double i = 0; i < n; i++) {
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, 0);
            double c = calculate_chi_from_cube(cube, q);
            //double c = chi_trapezoidal(q, T, mu, 60);
            Vec zero;
            double potential = V(zero, q, T, cube);
            //double potential = U*U * c / ( 1 - U * c) + U*U*U * c*c / ( 1 - U*U * c*c);
            file << q_mag << " " << potential << endl;
        }
        file << endl;
    }
}

void plot_chi2(double T) {
    ofstream file("chi_plot2.dat");
    file << "q chi " << endl;
    double n = 80.0;
    for (double mu = -0.6; mu > -4.0; mu--) {
        vector<vector<vector<double>>> cube = chi_cube(T, mu);
        for (double i = 0; i < n; i++) {
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            double c = calculate_chi_from_cube(cube, q);
            file << q_mag << " " << c << endl;
        }
        file << endl;
    }
}

void plot_chi4(double T) {
    ofstream file("chi_plot4.dat");
    file << "q chi " << endl;
    double n = 80.0;
    for (double mu = 0.0; mu > -4.0; mu--) {
        for (double i = 0; i < n; i++) {
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            double c = integrate_susceptibility(q, T, mu);
            file << q_mag << " " << c << endl;
        }
        file << endl;
    }
}

void test_chi_interpolate() {
    Vec q(1.0,-2.2,2.7);
    double mu = -1.0;
    double T = 0.25;
    auto cube = chi_cube(T, mu);
    double chi = calculate_chi_from_cube(cube, q);
    double chi2 = chi_trapezoidal(q, T, mu, 40);
    Vec q2 = to_IBZ_2(q);
    double chi3 = chi_trapezoidal(q2, T, mu, 40);
   cout << "Chi: " << chi << " Chi 2: " << chi2 << " Chi 3: " << chi3 << endl;
}

void chi_interpolation_average_error() {
    double T = 0.25;
    double mu = 0.0;
    srand( (unsigned)time( NULL ) );
    vector<vector<vector<double>>> cube = chi_cube(T, mu);

    double max_error = 0;
    double min_error = 10;
    double error_ave = 0;
    int n = 10000;
    for (int i = 0; i < n; i++) {
        double x = (double) rand()/RAND_MAX * 0.5 + 2.7;
        double y = (double) rand()/RAND_MAX * 0.5 + 2.7;
        double z = (double) rand()/RAND_MAX * 0.5 + 2.7;
        Vec q(x, y, z);
        double c = calculate_chi_from_cube(cube, q);
        double c_exact = chi_trapezoidal(q, T, mu, 40);
        double error = abs(c - c_exact);
        if (error > max_error) max_error = error;
        if (error < min_error) min_error = error;
        error_ave += error / n;
    }
    std::cout << std::setprecision(8);
    cout << "Max Error: " << max_error << endl
        << "Min Error: " << min_error << endl
        << "Average Error: " << error_ave << endl;
}

void plot_potential_q() {
    ofstream file("potential_q.dat");
    double T = 1.0/4.0;
    file << "q V " << endl;
    double mu = -1.0;
    double n = 60.0;
    for (double i = 0; i < n; i++) {
        double q_mag = i/(n-1) * M_PI;
        Vec q_sub(q_mag, q_mag, q_mag);
        double chi_sub = chi_trapezoidal(q_sub, T, mu, 20);
        double Vs1 = U*U * chi_sub / ( 1 - U * chi_sub);
        double Vs2 = U*U*U * chi_sub*chi_sub / ( 1 - U*U * chi_sub*chi_sub);
        file << q_mag*pow(3,0.5) << " " << Vs1+Vs2 << endl;
    }
    file << endl;
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

void test_chi_trap() {
    Vec q(3, 4.92492, 0.742601, 0.539429, false);
    double T = 1.503114;
    double mu = 0.0;
    double num_points = 20;
    cout << chi_trapezoidal(q, T, mu, num_points) << endl;
    //cout << chi(q, T, mu, num_points) << endl;
    //cout << trapezoidal_integration_3d(num_points, &test_function) << endl;
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

double trapezoid_area_from_points(Vec k1, Vec k2, Vec k3, Vec k4) {
    // Define distances 
    Vec k12 = k1 - k2; if (k12.cartesian) k12.to_spherical();
    Vec k23 = k2 - k3; if (k23.cartesian) k23.to_spherical();
    Vec k34 = k3 - k4; if (k34.cartesian) k34.to_spherical();
    Vec k14 = k1 - k4; if (k14.cartesian) k14.to_spherical();
    Vec k13 = k1 - k3; if (k13.cartesian) k13.to_spherical();
    double d12 = k12.vals(0);
    double d23 = k23.vals(0);
    double d34 = k34.vals(0);
    double d14 = k14.vals(0);
    double d13 = k13.vals(0);

    double s = (d12 + d23 + d34 + d14)/2; 
    double theta1 = acos( (pow(d12,2) + pow(d23,2) - pow(d13,2)) / (2*d12*d23));
    double theta2 = acos( (pow(d14,2) + pow(d34,2) - pow(d13,2)) / (2*d14*d34));
    double theta = theta1 + theta2;
    double area = pow((s-d12)*(s-d23)*(s-d34)*(s-d14) - d12*d23*d34*d14*pow(cos(theta/2),2), 0.5);

    return area;
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
                double chi1 = chi_trapezoidal(input, 0.25, -1.2, 40);
                double chi2 = chi_trapezoidal(output, 0.25, -1.2, 40);
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
    double chi1 = chi_trapezoidal(q1, T, mu, 30);
    double chi2 = chi_trapezoidal(q2, T, mu, 30);
    double chi3 = chi_trapezoidal(q3, T, mu, 30);
    Vec q1_t = to_IBZ_2(q1);
    Vec q2_t = to_IBZ_2(q2);
    Vec q3_t = to_IBZ_2(q3);
    double chi1_t = chi_trapezoidal(q1_t, T, mu, 30);
    double chi2_t = chi_trapezoidal(q2_t, T, mu, 30);
    double chi3_t = chi_trapezoidal(q3_t, T, mu, 30);
    cout << chi1 << " " << chi1_t << endl;
    cout << chi2 << " " << chi2_t << endl;
    cout << chi3 << " " << chi3_t << endl;
}

void test_diagonalization() {
    int n = 20;
    auto func = [](double x, double y, double a, double b) {
        return cos(x)*cos(y) + 0.0;
        return (cos(x)-cos(y))*(cos(a)-cos(b)) + 1.0;
    };
    MatrixXd V(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double theta1 = 6.28*i/n, phi1 = 3.14*j/n;
            double x = cos(theta1)*sin(phi1), y = sin(theta1)*sin(phi1);
            V(i,j) = func(theta1, 2*phi1, 0, 0);
            //for (int k = 0; k < n; k++) {
            //    for (int l = 0; l < n; l++) {
            //        double theta1 = 6.28*i/n, phi1 = 3.14*j/n;
            //        double theta2 = 6.28*k/n, phi2 = 3.14*l/n;
            //        double x = cos(theta1)*sin(phi1), y = sin(theta1)*sin(phi1);
            //        double a = cos(theta2)*sin(phi2), b = sin(theta2)*sin(phi2);
            //        int ind1 = i*n + j, ind2 = k*n + l;
            //        V(ind1,ind2) = func(x, y, a, b);
            //        //V(ind1,ind2) = func(theta1, phi1, theta2, phi2);
            //    }
            //}
        }
    }
    EigenSolver<MatrixXd> s(V);
    VectorXcd vals = s.eigenvalues() / (n);
    EigenSolver<MatrixXd>::EigenvectorsType vecs = s.eigenvectors();

    vector<EigAndVec> solutions = combine_eigs_and_vecs(vals.real(), vecs.real());
    sort(solutions.rbegin(), solutions.rend());

    ofstream file("1val.dat");
    ofstream file2("2val.dat");
    cout << solutions[0].eig << endl << solutions[1].eig << endl;
    file << solutions[0].vec;
    file2 << solutions[1].vec;
}

void test_pairing_interaction_parity(vector<Vec> &FS, double T, double mu) {
    auto cube = chi_cube(T, mu);
    Vec k2;
    double num = 10;
    for (double i = 0; i < num; i++) {
        for (double j = 0; j < num; j++) {
            for (double k = 0; k < num; k++) {
                double x = get_k(i, num);
                double y = get_k(j, num);
                double z = get_k(k, num);

                Vec k1(x, y, z);

                double V1 = potential_scalapino_cube(k1, k2, T, cube);
                double V2 = potential_scalapino_cube(-1*k1, k2, T, cube);

                if ( fabs(V1 - V2) > 0.00001) 
                    cout << k1 << "->" << V1 << ", " << V2 << endl;
            }
        }
    }

}

vector<double> get_wave_from_file() {
    ifstream file("../data/const3D_mu=-1.2_U=4.0_wD=1.0_n=6.dat");
    string line;
    getline(file, line);
    vector<double> wave;
    while (file) {
        double x, y, z, c1;
        file >> x >> y >> z >> c1;
        wave.push_back(c1);
        string line;
        getline(file, line);
    }
    return wave;
}

double test_surface_integral_discrete(vector<Vec> &FS_vecs, int i, double T, const vector<vector<vector<double>>> &chi_cube, auto &V, bool norm = false) {
    auto ind_func = [&V](auto &FS_vecs, auto &dumb_wave, int i, int j, double T, const vector<vector<vector<double>>> &chi_cube, bool norm = false) {
        Vec k1 = FS_vecs[i]; Vec k2 = FS_vecs[j];
        if (not norm) {
            return dumb_wave[i] * V(k1, k2, T, chi_cube) * dumb_wave[j] / (vp(k1) * vp(k2));
        }
        return dumb_wave[j] * dumb_wave[j];// / (vp(k2)*vp(k2));
    };
    vector<double> dumb_wave = get_wave_from_file();
    double sum = 0;
    for (int j = 0; j < FS_vecs.size(); j++) {
        Vec k1 = FS_vecs[j];
        double f = ind_func(FS_vecs, dumb_wave, i, j, T, chi_cube, norm);
        sum += f*k1.area;
    }
    return sum;
}

double test_4d_surface_integral_discrete(vector<Vec> &FS_vecs, const vector<vector<vector<double>>> &chi_cube, double T, auto &V) {
    
    double normalized = test_surface_integral_discrete(FS_vecs, 0, T, chi_cube, V, true);
    cout << "Discrete Normalization: " << normalized << endl;
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < FS_vecs.size(); i++) {
        Vec k = FS_vecs[i];
        double f = test_surface_integral_discrete(FS_vecs, i, T, chi_cube, V);
        sum += f*k.area;
    }
    return -sum / normalized;

}


double test_surface_integral4(vector<Vec> &FS_vecs, Vec p, double T, double mu, const vector<vector<vector<double>>> &chi_cube, auto &func) {
    double sum = 0;
    for (int i = 0; i < FS_vecs.size(); i++) {
        Vec k1 = FS_vecs[i];
        double f = func(k1, p, T, mu, chi_cube);
        sum += f*k1.area;
    }
    return sum * (2 / pow(2*M_PI,3));
}

double test_4d_surface_integral4(vector<Vec> &FS_vecs, const vector<vector<vector<double>>> &chi_cube, double T, double mu, auto &wave, auto &V) {

    auto norm_func = [&wave](Vec k, Vec extra, double T, double mu, const vector<vector<vector<double>>> &chi_cube) {
        return pow(wave(k, mu),2) / vp(k);
    };
    auto func = [&wave, &V](Vec p1, Vec p2, double T, double mu, const vector<vector<vector<double>>> &chi_cube) { 
        return wave(p1, mu) * V(p1,p2,T,chi_cube) * wave(p2, mu) / (vp(p1) * vp(p2)); 
    };

    double normalized = test_surface_integral4(FS_vecs, FS_vecs[0], T, mu, chi_cube, norm_func);
    //cout << "Normalized: " << normalized << endl;

    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < FS_vecs.size(); i++) {
        Vec k = FS_vecs[i];
        double f = test_surface_integral4(FS_vecs, k, T, mu, chi_cube, func);
        sum += f*k.area;
    }
    return -sum / normalized * (2 / pow(2*M_PI,3));
}


void plot_coupling() {

    auto dx2_y2 = [](Vec k, double mu) { 
        Vec q = k; 
        if (q.cartesian == false) q.to_cartesian();
        return cos(q.vals(1))-cos(q.vals(0)); 
    };
    auto dz2_1 = [](Vec k, double mu) { 
        Vec q = k; 
        if (q.cartesian == false) q.to_cartesian();
        return 2*cos(q.vals(2))-cos(q.vals(0))-cos(q.vals(1)); 
    };
    auto dxy = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(0))*sin(q.vals(1));
    };
    auto dxz = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(0))*sin(q.vals(2));
    };
    auto dyz = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(1))*sin(q.vals(2));
    };
    auto px = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(0));
    };
    auto py = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(1));
    };
    auto pz = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return sin(q.vals(2));
    };
    auto p_wave2 = [](Vec k, double mu) { 
        Vec q = k; if (not q.cartesian) q.to_cartesian();
        return (sin(q.vals(0)) - sin(q.vals(1)));
    };
    auto extended_s_wave = [](Vec k, double mu) { 
        return -mu/2 + cos(k.vals(0)) + cos(k.vals(1)) + cos(k.vals(2)); 
    };
    auto wave_norm = [](Vec k, double mu) {
        return 1.0; 
    };
    auto Vs = [](Vec p1, Vec p2, double T, const auto &chi_cube) { 
        return V(p1, p2, T, chi_cube);
    };
    auto Vt = [](Vec p1, Vec p2, double T, const auto &chi_cube) { 
        return 1.0;
        Vec q = to_IBZ_2(p1 - p2);
        double chi = calculate_chi_from_cube(chi_cube, q);
        return -pow(U,2)*chi / (1 - pow(U*chi,2));
    };
    auto V_norm = [](Vec p1, Vec p2, double T, const auto &chi_cube) { 
        Vec q = to_IBZ_2(p1 - p2);
        double chi = calculate_chi_from_cube(chi_cube, q);
        double V = pow(U,3)*pow(chi,2) / (1-U*chi) + pow(U,2)*chi / (1-pow(U*chi,2));
        return -V;
    };

    double T = 0.0065;
    T = 0.25;
    double mu = -0.9;
    ofstream file("coupling.dat");

    while (mu > -2.3) {
        cout << "Mu: " << mu << endl;
        vector<Vec> FS = tetrahedron_method(mu); 
        cout << "Fermi Surface Size: " << FS.size() << endl;
        double area = 0; for (auto x : FS) area += x.area;

        vector<vector<vector<double>>> cube = chi_cube(T, mu);
        //MatrixXd P = create_P(FS, T, cube);
        //vector<EigAndVec> solutions = power_iteration(P, 0.0001);
        cout << "Power Iteration Done\n";
        //cout << P.eigenvalues().real().maxCoeff();

        //double eig = solutions[solutions.size()-1].eig;
        double eig = 0;
        cout << "Eig Found: " << eig << " \n";

        // Renormalization
        //double renormalization = test_4d_surface_integral4(FS, cube, T, mu, wave_norm, V_norm)+1.0;
        double renormalization = 1;
        //cout << "\nEigs: " << endl;
        //for (auto x : solutions) {
        //    cout << x.eig / renormalization << endl;
        //}
        //cout << endl;
        // D-waves
        double dx2_y2_coupling = test_4d_surface_integral4(FS, cube, T, mu, dx2_y2, Vs);
        cout << dx2_y2_coupling << endl;
        //double dz2_1_coupling, dxy_coupling, dxz_coupling, dyz_coupling, px_coupling, py_coupling, pz_coupling, s_coupling, extended_s_coupling, dx_coupling = 0;
        double dz2_1_coupling = test_4d_surface_integral4(FS, cube, T, mu, dz2_1, Vs);
        double dxy_coupling = test_4d_surface_integral4(FS, cube, T, mu, dxy, Vs);
        double dxz_coupling = test_4d_surface_integral4(FS, cube, T, mu, dxz, Vs);
        double dyz_coupling = test_4d_surface_integral4(FS, cube, T, mu, dyz, Vs);
        // P-waves
        double px_coupling = test_4d_surface_integral4(FS, cube, T, mu, px, Vt);
        double py_coupling = test_4d_surface_integral4(FS, cube, T, mu, py, Vt);
        double pz_coupling = test_4d_surface_integral4(FS, cube, T, mu, pz, Vt);

        //// P-wave for singlet potential
        double dx_coupling = test_4d_surface_integral4(FS, cube, T, mu, px, Vs);
        ////double dy_coupling = test_4d_surface_integral4(FS, cube, T, mu, py, Vs);
        ////double dz_coupling = test_4d_surface_integral4(FS, cube, T, mu, pz, Vs);
        ////double p_coupling2 = test_4d_surface_integral4(FS, cube, T, mu, p_wave2, Vs);
        double s_coupling = test_4d_surface_integral4(FS, cube, T, mu, wave_norm, Vs);
        double extended_s_coupling = test_4d_surface_integral4(FS, cube, T, mu, extended_s_wave, Vs);
        //double dumb_wave_coupling = test_4d_surface_integral_discrete(FS, cube, T, Vs);
        //double dz2_1_coupling = 0, dxz_coupling = 0, dyz_coupling = 0, px_coupling = 0, py_coupling = 0, pz_coupling = 0, s_coupling = 0; 
        file << setprecision(4);
        file << mu << " ";
        file << renormalization << " ";

        file << dx2_y2_coupling / renormalization << " "<< " ";
        file << dz2_1_coupling / renormalization << " "<< " ";
        file << dxy_coupling / renormalization << " "<< " ";
        file << dxz_coupling / renormalization << " "<< " ";
        file << dyz_coupling / renormalization << " "<< " ";

        file << px_coupling / renormalization << " "<< " ";
        file << py_coupling / renormalization << " "<< " ";
        file << pz_coupling / renormalization << " "<< " ";

        //cout << dy_coupling / renormalization << " ";
        //cout << dz_coupling / renormalization << " ";
        //cout << "P-coupling2:  p_coupling2 / renormalization << " ";
        cout << s_coupling << " " << s_coupling / renormalization << endl;
        file << s_coupling / renormalization << " ";
        file << extended_s_coupling / renormalization  << " ";
        //cout << "Dumb Wave:  dumb_wave_coupling / renormalization  << " ";

        //double eig = 0;
        file << eig / renormalization << " " << " ";
        file << endl;
        mu -= 0.1;
    }

}

void test_integral_pairing_interaction_parity(vector<Vec> &FS, double T, double mu) {

    auto wave = [](Vec k, double mu) { 
        Vec q = k; 
        if (q.cartesian == false) q.to_cartesian();
        return cos(q.vals(1))-cos(q.vals(0)); 
    };
    auto Vs = [](Vec p1, Vec p2, double T, const auto &chi_cube) { 
        return V(p1, p2, T, chi_cube);
    };
    auto func = [&wave, &Vs](Vec p1, Vec p2, double T, double mu, const vector<vector<vector<double>>> &chi_cube) { 
        return wave(p1, mu) * Vs(p1,p2,T,chi_cube) * wave(p2, mu);// / (vp(p1) * vp(p2)); 
    };

    auto cube = chi_cube(T, mu);
    double num = 10;
    for (double i = 0; i < num; i++) {
        for (double j = 0; j < num; j++) {
            for (double k = 0; k < num; k++) {
                double x = 2*k_max*(i/num-1.0);
                double y = 2*k_max*(j/num-1.0);
                double z = 2*k_max*(k/num-1.0);

                Vec p(x, y, z);

                double val1 = test_surface_integral4(FS, p, T, mu, cube, func);
                double val2 = test_surface_integral4(FS, -1*p, T, mu, cube, func);
                if ( fabs(val1 - val2) > 0.0000001) {
                    cout << setprecision(10);
                    cout << "P: " << p << endl;
                    cout << "Integral: " << val1 << " " << val2 << endl;
                    cout << "Chi: " << calculate_chi_from_cube(cube, p) << endl;
                }
            }
        }
    }

}

double heavy_step(Vec k) {
    return k.vals.norm() < 1;
}

double ratio_1D(Vec q, Vec k) {
    double e_k = epsilon_sphere(k);
    double e_qk = epsilon_sphere(k+q);

    double f_k = heavy_step(k);
    double f_qk = heavy_step(k+q);

    return (f_qk - f_k) / (e_k - e_qk);
}

double explicit_formula_1d(Vec q, Vec k) {
    return 1 / (q.vals.squaredNorm() - 4*k.vals.squaredNorm());
}

double example_formula_1D(Vec q, Vec k) {
    //return (k - q).vals.norm();
    return 1/(k.vals(0) - q.vals(0));
}

bool k_in_dense_box(Vec k, vector<Vec> singularities, vector<Vec> singularity_widths) {
    for (int i = 0; i < singularities.size(); i++) {
        Vec singularity = singularities[i];
        Vec singularity_width = singularity_widths[i];
        if (k.vals(0) >= singularity.vals(0) - singularity_width.vals(0) 
            and k.vals(0) <= singularity.vals(0) + singularity_width.vals(0)
            and k.vals(1) >= singularity.vals(1) - singularity_width.vals(1) 
            and k.vals(1) <= singularity.vals(1) + singularity_width.vals(1) 
            and k.vals(2) >= singularity.vals(2) - singularity_width.vals(2) 
            and k.vals(2) <= singularity.vals(2) + singularity_width.vals(2)) {
            return true;
        }
    }
    return false;
}

int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    //double sigma = *((double *) fdata); // we can pass σ via fdata argument
    Vec q = *((Vec *) fdata);
    //double sum = 0;
    //unsigned i;
    //for (i = 0; i < ndim; ++i) sum += x[i] * x[i];
    //// compute the output value: note that fdim should == 1 from below
    //fval[0] = exp(-sigma * sum);
    Vec k_val(x[0], x[1], x[2]);
    double T = 0.25;
    fval[0] = ratio(q, k_val, T, mu);
    //cout << k_val << fval[0] << endl;
    //fval[0] = x[0]*x[0] * x[1]*x[1] * x[2]*x[2];
    return 0; // success*
}

void plot_DOS() {
    ofstream file("info.log");
    double num = 100;
    vector<double> DOS(num);
    double mu_max = 5.9;
    for (int i = 0; i < num; i++) {
        double new_mu = -mu_max + 2.0*mu_max*i/(num-1);
        init_config(mu, U, t, tn, new_mu, U, t, tn);

        auto wave = [](Vec k, double mu) { 
            return 1.0;
        };
        auto norm_func = [&wave](Vec k, Vec extra, double T, double mu, const vector<vector<vector<double>>> &chi_cube) {
            //return 1.0;
            return pow(wave(k, mu),2) / vp(k);
        };

        vector<Vec> FS = tetrahedron_method(mu);
        filter_FS(FS);
        double T = 25;
        vector<vector<vector<double>>> cube;
        double test = test_surface_integral4(FS, FS[0], T, mu, cube, norm_func);
        DOS[i] = test;
        file << new_mu << " " << DOS[i] << endl;
    }
}

int main() {
    plot_chi2(0.1); return 0;
    //plot_DOS();
    //return 0;
    double T = 0.25;
    //Vec q(3.14,3.14,3.14);
    //cout << integrate_susceptibility(q, T, mu);
    //return 0;
    //plot_potential(0.25);
    plot_coupling();
    //return 0;
    vector<Vec> FS = tetrahedron_method(mu);
    cout << FS.size() << endl;
    //filter_FS(FS);
    cout << FS.size() << endl;

    save_FS(FS);
    vector<vector<vector<double>>> cube;
    auto wave = [](Vec k, double mu) { 
        return 1.0;
    };
    auto norm_func = [&wave](Vec k, Vec extra, double T, double mu, const vector<vector<vector<double>>> &chi_cube) {
        return 1.0;
        return pow(wave(k, mu),2) / vp(k);
    };
    double test = test_surface_integral4(FS, FS[0], T, mu, cube, norm_func);
    cout << test << endl;
    return 0;
}
