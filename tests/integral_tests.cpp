#include "../gap/susceptibility.h"
#include "../gap/cfg.h"
#include "../gap/integration.h"


using namespace std;



bool compare_real_vs_complex_susceptibility_integration(Vec q, float T, float mu, float w, float num_points) {
    float c_real = integrate_susceptibility(q, T, mu, w, num_points);
    complex<float> c_imag = complex_susceptibility_integration(q, T, mu, w, num_points);
    if (abs(c_real - c_imag.real()) > 0.001) {
        return false;
    }
    return true;
}

bool check_nonzero_imaginary_part() {
    Vec q(1.0, 1.0, 1.0);
    float T = 0.25;
    float mu = 1.0;
    complex<float> w(0.0, 0.5);
    int num_points = 100;
    complex<float> c = complex_susceptibility_integration(q, T, mu, w, num_points);
    return (abs(c.imag()) > 0.001);
}

void complex_integration_convergence_test() {
    Vec q(1.0, 1.0, 1.0);
    float T = 0.005;
    float mu = 1.0;
    complex<float> w(1.0, 0.5);
    int num_points = 100;
    for (int i = 0; i < 100; i++) {
        int num_points = 100 + i * 10;
        complex<float> c = complex_susceptibility_integration(q, T, mu, w, num_points);
        cout << num_points << " " << c.real() << " " << c.imag() << endl;
    }
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
