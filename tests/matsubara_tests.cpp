#include <fstream>
#include <complex>
#include <cassert>
#include <openblas/lapacke.h>
#include "../gap/frequency_inclusion.hpp"
#include "../gap/cfg.h"
#include "../gap/susceptibility.h"

using namespace std;

bool check_matcube_interpolation() {
    float T = 0.25;
    float mu = 0.0;
    MatCube matsubara_cube = create_matsubara_cube(T, mu, m, w_pts, -max_freq, max_freq, 100);
    complex<float> w1(0.0, 0.0);
    complex<float> w2(0.0, 0.5);
    complex<float> w3(0.0, 1.0);
    complex<float> w4(0.0, M_PI);

    Vec q1(0.0, 0.0, 0.0); complex<float> c1 = matsubara_cube(q1, w1);
    Vec q2(0.5, 0.5, 0.5); complex<float> c2 = matsubara_cube(q2, w2);
    //Vec q3(1.0, 1.0, 1.0); complex<float> c3 = matsubara_cube(q3, w3);
    //Vec q4(M_PI, M_PI, M_PI); complex<float> c4 = matsubara_cube(q4, w4);

    complex<float> i1 = complex_susceptibility_integration(q1, T, mu, w1.imag(), 100);
    complex<float> i2 = complex_susceptibility_integration(q2, T, mu, w2.imag(), 100);
    //complex<float> i3 = complex_susceptibility_integration(q3, T, mu, w3.imag(), 100);
    //complex<float> i4 = complex_susceptibility_integration(q4, T, mu, w4.imag(), 100);

    printf("Matcube: %f %f \n", c1.real(), c2.real());
    printf("Integration: %f  %f\n", i1.real(), i2.real());
    return (abs(c1 - i1) < 0.001 and abs(c2 - i2) < 0.001 );
    //        and abs(c3 - i3) < 0.001 and abs(c4 - i4) < 0.001);

    return true;
}

bool check_susceptibility_integration_methods_are_equivalent(Vec q, float T, float mu, float w, int num_points) {
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


    float c1 = trapezoidal_integration(func, -k_max, k_max, -k_max, k_max, -k_max, k_max, num_points);
    c1 /= pow(2*M_PI, dim);

    float c2 = integrate_susceptibility(q, T, mu, w, num_points);
    MatCube matsubara_cube = create_matsubara_cube(T, mu, m, w_pts, -max_freq, max_freq, 100);
    complex<float> c3 = matsubara_cube(q, w);
    return (abs(c1 - c2) < 0.001 and abs(c1 - c3.real()) < 0.001);
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
    printf("Output written to real_vs_complex_integration.dat\n");
    return true;
}

bool compare_real_vs_complex_susceptibility_integration(Vec q, float T, float mu, float w, float num_points) {
    float c_real = integrate_susceptibility(q, T, mu, w, num_points);
    complex<float> c_imag = complex_susceptibility_integration(q, T, mu, w, num_points);
    if (abs(c_real - c_imag.real()) > 0.001) {
        return false;
    }
    return true;
}

bool test_real_integration() {
    auto func = [](float x, float y, float z) -> double {
        double c = 0.0;
        return 1 / (x + c);
    };
    float c = trapezoidal_integration(func, 1.0, 2.0, 0.0, 1.0, 0.0, 1.0, 300);
    return (abs(c - 0.69316) < 0.001);
}

bool test_complex_integration() {
    auto func = [](float x, float y, float z) -> complex<float> {
        complex<float> c = {0.0, 1.0};
        return (float)1.0/(x+c);
    };
    complex<float> c = complex_trapezoidal_integration(func, 1.0, 2.0, 0.0, 1.0, 0.0, 1.0, 300);
    return (abs(c.real() - 0.45815) < 0.001 and abs(c.imag() - (-0.32175)) < 0.001);
}
