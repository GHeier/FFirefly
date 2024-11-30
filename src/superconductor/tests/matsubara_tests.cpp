#include <fstream>
#include <complex>
#include <cassert>
#include <openblas/lapacke.h>
#include "../gap/frequency_inclusion.hpp"
#include "../gap/cfg.h"
#include "../gap/susceptibility.h"
#include "../gap/save_data.h"
#include "../gap/band_structure.h"
#include "../gap/integration.h"
#include "../gap/susceptibility.h"

using namespace std;

bool check_matcube_interpolation() {
    float T = 0.25;
    float mu = 0.0;
    MatCube matsubara_cube = create_matsubara_cube(T, mu, m, w_pts, -max_freq, max_freq, 100);
    save_matsubara_cube(matsubara_cube, matsubara_cube.w_min, matsubara_cube.w_max, "matsubara_cube.dat");

    complex<float> w1(0.0, 0.0);
    complex<float> w2(0.0, 0.5);
    complex<float> w3(0.0, 1.0);

    Vec q1(0.0, 0.0, 0.0); complex<float> c1 = matsubara_cube(q1, w1);
    Vec q2(0.5, 0.5, 0.5); complex<float> c2 = matsubara_cube(q2, w2);
    Vec q3(1.0, 1.0, 1.0); complex<float> c3 = matsubara_cube(q3, w3);

    complex<float> i1 = complex_susceptibility_integration(q1, T, mu, w1, 100);
    complex<float> i2 = complex_susceptibility_integration(q2, T, mu, w2, 100);
    complex<float> i3 = complex_susceptibility_integration(q3, T, mu, w3, 100);
    printf("Real and imaginary parts of Matsubara cube and integration methods\n");
    printf("i1: (%f, %f) c1: (%f, %f)\n", i1.real(), i1.imag(), c1.real(), c1.imag());
    printf("i2: (%f, %f) c2: (%f, %f)\n", i2.real(), i2.imag(), c2.real(), c2.imag());
    printf("i3: (%f, %f) c3: (%f, %f)\n", i3.real(), i3.imag(), c3.real(), c3.imag());

    return (abs(c1 - i1) < 0.001 and abs(c2 - i2) < 0.001 and abs(c3 - i3) < 0.001); 
}

bool check_susceptibility_integration_methods_are_equivalent(Vec q, float T, float mu, float w, int num_points) {
    auto func = [T, q, w, mu](float x, float y, float z) -> float {
        Vec k(x,y,z);
        float e_k = epsilon(k) - mu;
        float e_kq = epsilon(k+q) - mu;
        float f_kq = fermi_dirac(e_kq, T);
        float f_k = fermi_dirac(e_k, T);
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

