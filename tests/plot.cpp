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
#include "calculations.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
//using namespace Eigen;

void plot_chi(double T) {
    ofstream file("chi_plot.dat");
    file << "q chi " << endl;
    double n = 100.0;
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
    auto FS = tetrahedron_method(mu);
    double DOS = 0;
    for (int i = 0; i < FS.size(); i++)
        DOS += FS[i].area / vp(FS[i]);
    double n = 80.0;
    for (double mu = 0.0; mu > -2.0; mu--) {
        vector<vector<vector<double>>> cube = chi_cube(T, mu, DOS, 0);
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
    double n = 100.0;
    for (double mu = -0.0; mu > -4.0; mu-=1.0) {
        vector<Vec> FS = tetrahedron_method(mu);
        double DOS = get_DOS(FS);
        vector<vector<vector<double>>> cube = chi_cube(T, mu, DOS, 0);
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
    double n = 10.0;
    for (double mu = 0.0; mu > -4.0; mu--) {
        for (double i = 0; i < n; i++) {
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            double c = integrate_susceptibility(q, T, mu, 0);
            file << q_mag << " " << c << endl;
        }
        file << endl;
    }
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


