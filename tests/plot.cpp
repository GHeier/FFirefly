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

void plot_chi(double T, double w) {
    ofstream file("chi_plot.dat");
    file << "q chi " << endl;
    double n = 50.0;
    for (double new_mu = 0; new_mu > -4.0; new_mu--) {
        for (double i = 0; i < n; i++) {
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            double c = chi_trapezoidal(q, T, new_mu, w, 100);
            //double c2 = varied_chi_trapezoidal(q, T, mu, 40);
            file << q_mag << " " << c << endl;
        }
        file << endl;
    }
}

void plot_potential(double T) {
    ofstream file("info.log");
    file << "q V " << endl;
    auto FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
    double DOS = 0;
    for (int i = 0; i < FS.size(); i++)
        DOS += FS[i].area / vp(FS[i]);
    double n = 80.0;
    for (double mu = 0.0; mu > -2.0; mu--) {
        unordered_map<double, vector<vector<vector<double>>>> cube = chi_cube_freq(T, mu, DOS);
        for (double i = 0; i < n; i++) {
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, 0);
            Vec zero;
            double potential = V(zero, q, 0, T, cube);
            //double potential = U*U * c / ( 1 - U * c) + U*U*U * c*c / ( 1 - U*U * c*c);
            file << q_mag << " " << potential << endl;
        }
        file << endl;
    }
}

void plot_chi2(double T) {
    ofstream file("chi_plot2.dat");
    file << "q chi " << endl;
    double n = 50.0;
    for (double mu = -0.0; mu > -4.0; mu-=1.0) {
        vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
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

void plot_chi3(double T, double w, double eta, int pts) {
    ofstream file("chi_plot3.dat");
    file << "q chi " << endl;
    double n = 50.0;
    for (double mu = 0.0; mu < 4.0; mu++) {
        for (double i = 0; i < n; i++) {
            printf("\r Mu %.1f: %.3f" , mu, i/(n-1));
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            double c = imaginary_integration(q, T, mu, w, pts, eta);
            //double c = integrate_susceptibility(q, T, mu, 0, 3*pts);
            //double c = modified_integral_wrapper(q, T, mu, 0, delta, pts);
            file << q_mag << " " << c << endl;
        }
        file << endl;
    }
}

void plot_chi4(double T, double w, int pts) {
    ofstream file("chi_plot4.dat");
    file << "q chi " << endl;
    Vec q_temp(M_PI, M_PI, M_PI);
    double n = 50.0;
    for (double mu = 0.0; mu > -4.0; mu--) {
        for (double i = 0; i < n; i++) {
            printf("\r Mu %.1f: %.3f" , mu, 100.0*i/(n-1));
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            double c = integrate_susceptibility(q, T, mu, w, pts);
            //cout << c << endl;
            //cout << integrate_susceptibility(q_temp,T,mu,0,30) << endl;
            file << q_mag << " " << c << endl;
            fflush(stdout);
        }
        file << endl;
    }
}

void plot_chi5(double T, double w) {
    ofstream file("chi_plot5.dat");
    file << "q chi " << endl;
    double n = 50.0;
    for (double new_mu = 0.0; new_mu > -4.0; new_mu--) {
        init_config(mu, U, t, tn, w_D, new_mu, U, t, tn, w_D);
        vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
        for (double i = 0; i < n; i++) {
            printf("\r Mu %.1f: %.3f" , new_mu, 100.0*i/(n-1));
            fflush(stdout);
            double q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            double c = chi_ep_integrate(q, w, T);
            //cout << c << endl;
            //cout << integrate_susceptibility(q_temp,T,mu,0,30) << endl;
            file << q_mag << " " << c << endl;
        }
        printf("\n");
        file << endl;
    }
}

void plot_single_chi(double T, double w) {
    ofstream file("single_chi_plot.dat");
    file << "q chi " << endl;
    double n = 50.0;
    double new_mu = 1.0;
    init_config(mu, U, t, tn, w_D, new_mu, U, t, tn, w_D);
    for (double i = 0; i < n; i++) {
        printf("\r Mu %.1f: %.3f" , new_mu, 100.0*i/(n-1));
        fflush(stdout);
        double q_mag = i/(n-1) * M_PI;
        Vec q(q_mag, q_mag, q_mag);
        double c = chi_trapezoidal(q, T, new_mu, w, 100);
        //double c2 = varied_chi_trapezoidal(q, T, mu, 40);
        file << q_mag << " " << c << endl;
    }

}

void plot_single_chi2(double T, double w) {
    ofstream file("single_chi_plot2.dat");
    file << "q chi " << endl;
    double n = 50.0;
    double new_mu = 1.0;
    init_config(mu, U, t, tn, w_D, new_mu, U, t, tn, w_D);
    for (double i = 0; i < n; i++) {
        printf("\r Mu %.1f: %.3f" , new_mu, 100.0*i/(n-1));
        fflush(stdout);
        double q_mag = i/(n-1) * M_PI;
        Vec q(q_mag, q_mag, q_mag);
        double c = chi_ep_integrate(q, w, T);
        file << q_mag << " " << c << endl;
    }

}

void plot_surfaces(Vec q, double T, double w) {
    double a, b;
    get_bounds3(q, b, a, denominator);
    a *= 0.99; b*= 0.99;
    b = (b - a) * 0.5;
    a = -b;
    printf("Bounds: %.3f, %.3f\n", a, b);
    double A, upr, lwr;
    get_spacing_curve_consts(w, a, b, A, upr, lwr);
    auto spacing = [A, w, lwr, upr] (double i, double pts) { 
        double x = lwr + (upr- lwr) * i / pts;
        printf("X: %.3f\n", x);
        return A * x + w;
        return A * pow(x,3) + w; 
    };

    int pts = 40;
    for (int i = 0; i <= pts; i++) {
        double s = spacing(i, pts);
        ofstream file("test_surface" + to_string(i) + ".dat");
        vector<Vec> surface = tetrahedron_method(denominator, q, s);
        printf("S: %.3f, size: %d\n", s, surface.size());
        for (auto x: surface) {
            file << x << endl;
        }
    }

}

void plot_surfaces2(Vec q, double T, double w) {
    double a, b;
    get_bounds3(q, b, a, denominator);

    vector<double> spacing;
    get_spacing_vec(spacing, w, a, b, 20);

    for (int i = 0; i < spacing.size(); i++) {
        ofstream file("test_surface" + to_string(i) + ".dat");
        vector<Vec> surface = tetrahedron_method(denominator, q, spacing[i]);
        //printf("S: %.3f, size: %d\n", spacing[i], surface.size());
        for (auto x: surface) {
            file << x << endl;
        }
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
        double chi_sub = chi_trapezoidal(q_sub, T, mu, 0, 20);
        double Vs1 = U*U * chi_sub / ( 1 - U * chi_sub);
        double Vs2 = U*U*U * chi_sub*chi_sub / ( 1 - U*U * chi_sub*chi_sub);
        file << q_mag*pow(3,0.5) << " " << Vs1+Vs2 << endl;
    }
    file << endl;
}


