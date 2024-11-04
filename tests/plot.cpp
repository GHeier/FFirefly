#include <iostream>
#include <fstream>

#include <cstring>
#include <string>
#include <math.h>
#include <complex>

#include <boost/functional/hash.hpp>
#include <gsl/gsl_integration.h>

#include "../gap/cfg.h"
#include "../gap/fermi_surface.h"
#include "../gap/vec.h"
#include "../gap/analysis.h"
#include "../gap/save_data.h"
#include "../gap/calculations.h"
#include "../gap/susceptibility.h"
#include "../gap/save_data.h"
#include "../gap/utilities.h"

using namespace std;

void plot_real_susceptibility_integration(float w) {
    float T = 0.25;
    float mu = 0.0;
    int num_points = 100;
    ofstream file("susceptibility_integration.dat");
    for (int i = 0; i < 50; i++) {
        float mag = M_PI * i / 49.0;
        Vec q(mag, mag, mag);
        float c = integrate_susceptibility(q, T, mu, w, num_points);
        file << mag << " " << c << endl;
    }
    printf("Output written to susceptibility_integration.dat\n");
}

void plot_complex_susceptibility_integration(complex<float> w) {
    float T = 0.25;
    float mu = 0.0;
    int num_points = 100;
    ofstream file("complex_susceptibility_integration.dat");
    for (int i = 0; i < 50; i++) {
        float mag = M_PI * i / 49.0;
        Vec q(mag, mag, mag);
        complex<float> c = complex_susceptibility_integration(q, T, mu, w, num_points);
        file << mag << " " << c.real() << " " << c.imag() << endl;
    }
    printf("Output written to complex_susceptibility_integration.dat\n");
}

void plot_complex_susceptibility_integration_v_w(Vec q) {
    float T = 0.005;
    float mu = 1.0;
    int num_points = 2000;
    ofstream file("complex_susceptibility_integration.dat");
    for (int i = 0; i < 50; i++) {
        float mag = max_freq * i / 49.0;
        complex<float> w(mag, 0.01);
        complex<float> c = complex_susceptibility_integration(q, T, mu, w, num_points);
        file << mag << " " << c.real() << " " << c.imag() << endl;
    }
    printf("Output written to complex_susceptibility_integration.dat\n");
}

void plot_analytic_susceptibility_integration(float w) {
    float mu = 1.0;
    int num_points = 300;
    ofstream file("analytic_susceptibility_integration.dat");
    for (int i = 1; i < 50; i++) {
        float mag = M_PI * i / 49.0;
        Vec q(mag, mag, mag);
        float c = analytic_tetrahedron_sum(q, w, num_points);
        file << mag << " " << c << endl;
    }
    printf("Output written to analytic_susceptibility_integration.dat\n");
}

void plot_real_trapezoidal_susceptibility_integration(float w) {
    float T = 0.25;
    float mu = 0.0;
    int num_points = 100;
    ofstream file("trap_susceptibility_integration.dat");
    for (int i = 0; i < 50; i++) {
        float mag = M_PI * i / 49.0;
        Vec q(mag, mag, mag);

        auto func = [T, q, w, mu](float x, float y, float z) -> float {
            Vec k(x,y,z);
            float e_k = epsilon(k) - mu;
            float e_kq = epsilon(k+q) - mu;
            float f_kq = f(e_kq, T);
            float f_k = f(e_k, T);
            if (fabs(e_kq - e_k) < 0.0001 and fabs(w) < 0.0001) {
                if (T == 0 or exp(e_k/T) > 1e6) return e_k < 0;
                return 1/T * exp(e_k/T) / pow( exp(e_k/T) + 1,2);
            }
            return (f_kq - f_k) / (w - (e_kq - e_k));
        };


        float c = trapezoidal_integration(func, -k_max, k_max, -k_max, k_max, -k_max, k_max, num_points);
        c /= pow(2*M_PI, dim);
        file << mag << " " << c << endl;
    }
    printf("Output written to trap_susceptibility_integration.dat\n");
}

void plot_matsubara_cube_v_q(complex<float> w) {
    float T = 0.005;
    float mu = 1.0;
    MatCube matsubara_cube = create_matsubara_cube(T, mu, m, w_pts, -max_freq, max_freq, 100);
    save_matsubara_cube(matsubara_cube, matsubara_cube.w_min, matsubara_cube.w_max, "matsubara_cube.dat");
    ofstream file("matsubara_test_plot.dat");

    // Plot along q = (pi, pi, pi)
    for (int i = 0; i < 100; i++) {
        float mag = M_PI * i / 99.0;
        Vec q(mag, mag, mag);
        complex<float> c = matsubara_cube(q, w);
        file << mag << " " << c.real() << " " << c.imag() << endl;
    }
    printf("Output written to matsubara_test_plot.dat\n");
}

void plot_matsubara_cube_v_w(Vec q) {
    float T = 0.005;
    float mu = 1.0;
    MatCube matsubara_cube = create_matsubara_cube(T, mu, m, w_pts, -max_freq, max_freq, 100);
    save_matsubara_cube(matsubara_cube, matsubara_cube.w_min, matsubara_cube.w_max, "matsubara_cube.dat");
    ofstream file("matsubara_test_w_plot.dat");
    for (int i = 0; i < 100; i++) {
        float mag = max_freq * i / 99.0;
        complex<float> w(mag, 0.01);
        complex<float> c = matsubara_cube(q, w);
        file << mag << " " << c.real() << " " << c.imag() << endl;
    }
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


