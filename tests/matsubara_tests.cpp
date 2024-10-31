#include <fstream>
#include <complex>
#include <cassert>
#include <openblas/lapacke.h>
#include "../gap/frequency_inclusion.hpp"
#include "../gap/cfg.h"
#include "../gap/susceptibility.h"


void plot_matsubara_cube(complex<float> w) {
    MatCube matsubara_cube = create_matsubara_cube(0.25, 0.0, 10, 10, -10, 10, 100);
    ofstream file("matsubara_test_plot.dat");

    // Plot along q = (pi, pi, pi)
    for (int i = 0; i < 100; i++) {
        float mag = M_PI * i / 99.0;
        Vec q(mag, mag, mag);
        complex<float> c = matsubara_cube(q, w);
        file << mag << " " << c.real() << " " << c.imag() << endl;
    }
}

bool compare_real_vs_complex_susceptibility() {
    float T = 0.25;
    float mu = 0.0;
    int num_points = 100;
    MatCube matsubara_cube = create_matsubara_cube(T, mu, 10, 10, -10, 10, num_points);
    complex<float> w(0.0, 0.0);
    ofstream file("real_vs_complex_integration.dat");
    for (int i = 0; i < 50; i++) {
        float mag = M_PI * i / 99.0;
        Vec q(mag, mag, mag);
        complex<float> c_complex = matsubara_cube(q, w);
        float c_real = integrate_susceptibility(q, T, mu, w.real(), num_points); 
        file << mag << " " << c_real << " " << c_complex.real() << endl; 
        //cout << "Real: " << c_real << " Imaginary: " << c_complex.real() << endl;
        if (abs(c_real - c_complex.real()) > 0.001) {
            return false;
        }
    }
    return true;
}

bool test_real_integration() {
    auto func = [](float x, float y, float z) -> double {
        double c = 1.0;
        return c / (x + c);
    };
    printf("Real Integration Test\n");
    float c = trapezoidal_integration(func, 1.0, 2.0, 0.0, 1.0, 0.0, 1.0, 1000);
    cout << "Real Integration Test: " << c << endl;
    return (abs(c - 0.69316) < 0.001);
}

bool test_complex_integration() {
    auto func = [](float x, float y, float z) -> complex<float> {
        complex<float> c = {1.0, 0.0};
        return (float)1 / (x + c);
    };
    complex<float> c = complex_trapezoidal_integration(func, 1.0, 2.0, 0.0, 1.0, 0.0, 1.0, 1000);
    cout << "Complex Integration Test: " << c.real() << " " << c.imag() << endl;
    return (abs(c.real() - 0.45815) < 0.001 and abs(c.imag() - (-0.32175)) < 0.001);
}
