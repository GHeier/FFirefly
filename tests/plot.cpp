#include <iomanip>
#include <iostream>
#include <fstream>
//#include <chrono>

#include <cstring>
#include <string>
#include <math.h>
#include <vector>
#include <complex>

#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <gsl/gsl_integration.h>
#include <functional>

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

void plot_chi(float T, float w) {
    ofstream file("chi_plot.dat");
    float n = 50.0;
    for (float new_mu = 0; new_mu > -4.0; new_mu--) {
        for (float i = 0; i < n; i++) {
            float q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            auto f = [T, w, q] (float x, float y, float z) -> float {
                Vec k(x, y, z);
                return ratio(k, q, w, T);
            };
            float c = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60) / pow(2*M_PI, dim);
            //float c2 = varied_chi_trapezoidal(q, T, mu, 40);
            file << q_mag << " " << c << endl;
        }
        file << endl;
    }
}

void plot_complex_chi(float T, complex<float> w) {
    ofstream file("chi_complex_plot.dat");
    float n = 500.0;
    for (float new_mu = 0; new_mu > -4.0; new_mu--) {
        for (float i = 0; i < n; i++) {
            float q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            complex<float> c = complex_susceptibility_integration(q, T, mu, w, 100);
            //float c2 = varied_chi_trapezoidal(q, T, mu, 40);
            file << q_mag << " " << c.real() << " " << c.imag() << endl;
            cout << q_mag << " " << c.real() << " " << c.imag() << endl;
        }
        file << endl;
    }
}

void plot_potential(float T) {
    ofstream file("info.log");
    file << "q V " << endl;
    auto FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
    float DOS = 0;
    for (int i = 0; i < FS.size(); i++)
        DOS += FS[i].area / vp(FS[i]);
    float n = 80.0;
    for (float mu = 0.0; mu > -2.0; mu--) {
        unordered_map<float, vector<vector<vector<float>>>> cube = chi_cube_freq(T, mu);
        for (float i = 0; i < n; i++) {
            float q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, 0);
            Vec zero;
            float potential = V(zero, q, 0, T, cube);
            //float potential = U*U * c / ( 1 - U * c) + U*U*U * c*c / ( 1 - U*U * c*c);
            file << q_mag << " " << potential << endl;
        }
        file << endl;
    }
}

void plot_chi2(float T) {
    ofstream file("chi_plot2.dat");
    file << "q chi " << endl;
    float n = 50.0;
    for (float mu = -0.0; mu > -4.0; mu-=1.0) {
        vector<Vec> FS = tetrahedron_method(e_base_avg, Vec(0,0,0), mu);
        float DOS = get_DOS(FS);
        vector<vector<vector<float>>> cube = chi_cube(T, mu, DOS, 0);
        for (float i = 0; i < n; i++) {
            float q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            float c = calculate_chi_from_cube(cube, q);
            file << q_mag << " " << c << endl;
        }
        file << endl;
    }
}

void plot_chi3(float T, float w, float eta, int pts) {
    ofstream file("chi_plot3.dat");
    file << "q chi " << endl;
    float n = 50.0;
    for (float mu = 0.0; mu < 4.0; mu++) {
        for (float i = 0; i < n; i++) {
            printf("\r Mu %.1f: %.3f" , mu, i/(n-1));
            float q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            float c = integrate_susceptibility(q, T, mu, 0, 3*pts);
            //float c = modified_integral_wrapper(q, T, mu, 0, delta, pts);
            file << q_mag << " " << c << endl;
        }
        file << endl;
    }
}

void plot_chi4(float T, float w, int pts) {
    ofstream file("chi_plot4.dat");
    file << "q chi " << endl;
    Vec q_temp(M_PI, M_PI, M_PI);
    float n = 50.0;
    for (float mu = 0.0; mu > -4.0; mu--) {
        for (float i = 0; i < n; i++) {
            printf("\r Mu %.1f: %.3f" , mu, 100.0*i/(n-1));
            float q_mag = i/(n-1) * M_PI;
            Vec q(q_mag, q_mag, q_mag);
            float c = integrate_susceptibility(q, T, mu, w, pts);
            //cout << c << endl;
            //cout << integrate_susceptibility(q_temp,T,mu,0,30) << endl;
            file << q_mag << " " << c << endl;
            fflush(stdout);
        }
        file << endl;
    }
}

void plot_single_chi(float T, float w) {
    ofstream file("single_chi_plot.dat");
    file << "q chi " << endl;
    float n = 50.0;
    float new_mu = 1.0;
    init_config(mu, U, t, tn, wc, new_mu, U, t, tn, wc);
    for (float i = 0; i < n; i++) {
        printf("\r Mu %.1f: %.3f" , new_mu, 100.0*i/(n-1));
        fflush(stdout);
        float q_mag = i/(n-1) * M_PI;
        Vec q(q_mag, q_mag, q_mag);
        auto f = [T, w, q] (float x, float y, float z) {
            Vec k(x, y, z);
            return integrand(k, q, w, T);
        };
        float c = trapezoidal_integration(f, -M_PI, M_PI, -M_PI, M_PI, -M_PI, M_PI, 60);
        file << q_mag << " " << c << endl;
    }
    cout << endl;
}

void plot_single_chi2(float T, float w) {
    ofstream file("single_chi_plot2.dat");
    file << "q chi " << endl;
    float n = 50.0;
    float new_mu = 1.0;
    init_config(mu, U, t, tn, wc, new_mu, U, t, tn, wc);
    vector<vector<vector<float>>> cube = chi_cube(T, mu, w, "Cube 1/1");
    for (float i = 0; i < n; i++) {
        printf("\r Mu %.1f: %.3f" , new_mu, 100.0*i/(n-1));
        fflush(stdout);
        float q_mag = i/(n-1) * M_PI;
        Vec q(q_mag, q_mag, q_mag);
        float c = calculate_chi_from_cube(cube, q);
        file << q_mag << " " << c << endl;
    }
    cout << endl;
}

void plot_single_chi3(float T, float w) {
    ofstream file("single_chi_plot3.dat");
    file << "q chi " << endl;
    float n = 50.0;
    float new_mu = 1.0;
    init_config(mu, U, t, tn, wc, new_mu, U, t, tn, wc);
    for (float i = 0; i < n; i++) {
        printf("\r Mu %.1f: %.3f" , new_mu, 100.0*i/(n-1));
        fflush(stdout);
        float q_mag = i/(n-1) * M_PI;
        Vec q(q_mag, q_mag, q_mag);
        float c = integrate_susceptibility(q, T, mu, w, s_pts);
        if (isnan(c)) {
            cout << "NAN: " << q << endl;
            assert(false);
        }
        file << q_mag << " " << c << endl;
    }
    cout << endl;
}

void plot_surfaces(Vec q, float T, float w) {
    float a, b;
    get_bounds(q, b, a, denominator);
    a *= 0.99; b*= 0.99;
    b = (b - a) * 0.5;
    a = -b;
    printf("Bounds: %.3f, %.3f\n", a, b);
    float A, upr, lwr;
    get_spacing_curve_consts(w, a, b, A, upr, lwr);
    auto spacing = [A, w, lwr, upr] (int i, int pts) { 
        float x = lwr + (upr- lwr) * i / pts;
        printf("X: %.3f\n", x);
        return A * x + w;
    };

    int pts = 40;
    for (int i = 0; i <= pts; i++) {
        float s = spacing(i, pts);
        ofstream file("test_surface" + to_string(i) + ".dat");
        vector<Vec> surface = tetrahedron_method(denominator, q, s);
        printf("S: %.3f, size: %d\n", s, surface.size());
        for (auto x: surface) {
            file << x << endl;
        }
    }

}

void plot_surfaces2(Vec q, float T, float w) {
    float a, b;
    get_bounds(q, b, a, denominator);

    vector<float> spacing;
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


