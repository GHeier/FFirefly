#include "src/algorithms/fft.hpp"
#include "src/config/load/c_config.h"
#include "fft_tests.hpp"
#include <cmath>
#include <iostream>

using namespace std;

const double TEST_TOLERANCE = 1e-10;

/**
 * @brief Test 1D FFT with a delta function
 * Delta function in time domain -> constant in frequency domain
 */
bool test_1d_delta() {
    int n = 8;
    vector<complex<double>> data(n, 0.0);
    data[0] = 1.0; // Delta at t=0

    vector<complex<double>> result = fft_1d(data, false, false);

    // All frequency components should be 1.0
    for (int i = 0; i < n; i++) {
        if (abs(result[i] - 1.0) > TEST_TOLERANCE) {
            printf("1D Delta test failed at index %d: expected 1.0, got (%f, %f)\n",
                   i, result[i].real(), result[i].imag());
            return false;
        }
    }

    return true;
}

/**
 * @brief Test 1D FFT with a constant function
 * Constant in time domain -> delta in frequency domain
 */
bool test_1d_constant() {
    int n = 8;
    vector<complex<double>> data(n, 1.0);

    vector<complex<double>> result = fft_1d(data, false, false);

    // Only DC component should be non-zero (= n)
    if (abs(result[0] - complex<double>(n, 0.0)) > TEST_TOLERANCE) {
        printf("1D Constant test failed at DC: expected %d, got (%f, %f)\n",
               n, result[0].real(), result[0].imag());
        return false;
    }

    for (int i = 1; i < n; i++) {
        if (abs(result[i]) > TEST_TOLERANCE) {
            printf("1D Constant test failed at index %d: expected 0, got (%f, %f)\n",
                   i, result[i].real(), result[i].imag());
            return false;
        }
    }

    return true;
}

/**
 * @brief Test 1D FFT inverse transform
 * FFT followed by IFFT should return original data
 */
bool test_1d_inverse() {
    int n = 16;
    vector<complex<double>> data(n);

    // Create some arbitrary data
    for (int i = 0; i < n; i++) {
        data[i] = complex<double>(sin(2.0 * M_PI * i / n), cos(2.0 * M_PI * i / n));
    }

    vector<complex<double>> original = data;
    vector<complex<double>> fft_result = fft_1d(data, false, false);
    vector<complex<double>> ifft_result = fft_1d(fft_result, true, false);

    for (int i = 0; i < n; i++) {
        if (abs(ifft_result[i] - original[i]) > TEST_TOLERANCE) {
            printf("1D Inverse test failed at index %d: expected (%f, %f), got (%f, %f)\n",
                   i, original[i].real(), original[i].imag(),
                   ifft_result[i].real(), ifft_result[i].imag());
            return false;
        }
    }

    return true;
}

/**
 * @brief Test 1D FFT with a cosine function
 * cos(2πkt/n) should give peaks at k and n-k
 */
bool test_1d_cosine() {
    int n = 16;
    int k = 3; // Frequency
    vector<complex<double>> data(n);

    for (int i = 0; i < n; i++) {
        data[i] = cos(2.0 * M_PI * k * i / n);
    }

    vector<complex<double>> result = fft_1d(data, false, false);

    // Should have peaks at k and n-k, each with amplitude n/2
    double expected_amplitude = n / 2.0;

    for (int i = 0; i < n; i++) {
        if (i == k || i == n - k) {
            if (abs(abs(result[i]) - expected_amplitude) > TEST_TOLERANCE) {
                printf("1D Cosine test failed at index %d: expected amplitude %f, got %f\n",
                       i, expected_amplitude, abs(result[i]));
                return false;
            }
        } else {
            if (abs(result[i]) > TEST_TOLERANCE) {
                printf("1D Cosine test failed at index %d: expected 0, got (%f, %f)\n",
                       i, result[i].real(), result[i].imag());
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Test in-place vs copy modes
 */
bool test_1d_in_place() {
    int n = 8;
    vector<complex<double>> data(n);

    for (int i = 0; i < n; i++) {
        data[i] = complex<double>(i, 0);
    }

    vector<complex<double>> data_copy = data;

    // In-place transform
    vector<complex<double>> result_in_place = fft_1d(data, false, true);

    // Copy transform
    vector<complex<double>> result_copy = fft_1d(data_copy, false, false);

    // Original data_copy should be unchanged
    for (int i = 0; i < n; i++) {
        if (data_copy[i] != complex<double>(i, 0)) {
            printf("In-place test failed: original data was modified in copy mode\n");
            return false;
        }
    }

    // Results should be identical
    for (int i = 0; i < n; i++) {
        if (abs(result_in_place[i] - result_copy[i]) > TEST_TOLERANCE) {
            printf("In-place test failed at index %d: in-place and copy results differ\n", i);
            return false;
        }
    }

    return true;
}

/**
 * @brief Test 2D FFT with separable functions
 * f(x,y) = cos(2πkx/nx) * cos(2πky/ny)
 */
bool test_2d_separable() {
    int nx = 8, ny = 8;
    int kx = 2, ky = 3;
    vector<complex<double>> data(nx * ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double val = cos(2.0 * M_PI * kx * i / nx) * cos(2.0 * M_PI * ky * j / ny);
            data[i * ny + j] = val;
        }
    }

    vector<complex<double>> result = fft_2d(data, nx, ny, false, false);

    // Should have peaks at (kx, ky), (kx, ny-ky), (nx-kx, ky), (nx-kx, ny-ky)
    double expected_amplitude = (nx * ny) / 4.0;

    int peak_count = 0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            bool is_peak = (i == kx || i == nx - kx) && (j == ky || j == ny - ky);

            if (is_peak) {
                if (abs(abs(result[i * ny + j]) - expected_amplitude) > TEST_TOLERANCE) {
                    printf("2D Separable test failed at (%d,%d): expected amplitude %f, got %f\n",
                           i, j, expected_amplitude, abs(result[i * ny + j]));
                    return false;
                }
                peak_count++;
            } else {
                if (abs(result[i * ny + j]) > TEST_TOLERANCE) {
                    printf("2D Separable test failed at (%d,%d): expected 0, got (%f, %f)\n",
                           i, j, result[i * ny + j].real(), result[i * ny + j].imag());
                    return false;
                }
            }
        }
    }

    if (peak_count != 4) {
        printf("2D Separable test failed: expected 4 peaks, found %d\n", peak_count);
        return false;
    }

    return true;
}

/**
 * @brief Test 2D FFT inverse
 */
bool test_2d_inverse() {
    int nx = 8, ny = 8;
    vector<complex<double>> data(nx * ny);

    for (int i = 0; i < nx * ny; i++) {
        data[i] = complex<double>(sin(i * 0.5), cos(i * 0.3));
    }

    vector<complex<double>> original = data;
    vector<complex<double>> fft_result = fft_2d(data, nx, ny, false, false);
    vector<complex<double>> ifft_result = fft_2d(fft_result, nx, ny, true, false);

    for (int i = 0; i < nx * ny; i++) {
        if (abs(ifft_result[i] - original[i]) > TEST_TOLERANCE) {
            printf("2D Inverse test failed at index %d\n", i);
            return false;
        }
    }

    return true;
}

/**
 * @brief Test 3D FFT with delta function
 */
bool test_3d_delta() {
    int nx = 4, ny = 4, nz = 4;
    vector<complex<double>> data(nx * ny * nz, 0.0);
    data[0] = 1.0; // Delta at origin

    vector<complex<double>> result = fft_3d(data, nx, ny, nz, false, false);

    // All frequency components should be 1.0
    for (int i = 0; i < nx * ny * nz; i++) {
        if (abs(result[i] - 1.0) > TEST_TOLERANCE) {
            printf("3D Delta test failed at index %d\n", i);
            return false;
        }
    }

    return true;
}

/**
 * @brief Test 3D FFT inverse
 */
bool test_3d_inverse() {
    int nx = 4, ny = 4, nz = 4;
    vector<complex<double>> data(nx * ny * nz);

    for (int i = 0; i < nx * ny * nz; i++) {
        data[i] = complex<double>(sin(i * 0.2), cos(i * 0.1));
    }

    vector<complex<double>> original = data;
    vector<complex<double>> fft_result = fft_3d(data, nx, ny, nz, false, false);
    vector<complex<double>> ifft_result = fft_3d(fft_result, nx, ny, nz, true, false);

    for (int i = 0; i < nx * ny * nz; i++) {
        if (abs(ifft_result[i] - original[i]) > TEST_TOLERANCE) {
            printf("3D Inverse test failed at index %d\n", i);
            return false;
        }
    }

    return true;
}

/**
 * @brief Test nD FFT with 4D data
 */
bool test_4d_inverse() {
    vector<int> dims = {4, 4, 4, 2};
    int total_size = 4 * 4 * 4 * 2;
    vector<complex<double>> data(total_size);

    for (int i = 0; i < total_size; i++) {
        data[i] = complex<double>(sin(i * 0.3), cos(i * 0.2));
    }

    vector<complex<double>> original = data;
    vector<complex<double>> fft_result = fft_nd(data, dims, false, false);
    vector<complex<double>> ifft_result = fft_nd(fft_result, dims, true, false);

    for (int i = 0; i < total_size; i++) {
        if (abs(ifft_result[i] - original[i]) > TEST_TOLERANCE) {
            printf("4D Inverse test failed at index %d\n", i);
            return false;
        }
    }

    return true;
}

/**
 * @brief Test real input FFT
 */
bool test_1d_real_input() {
    int n = 8;
    vector<double> data(n);

    for (int i = 0; i < n; i++) {
        data[i] = cos(2.0 * M_PI * 2 * i / n);
    }

    vector<complex<double>> result = fft_1d_real(data, false);

    // Should have peaks at k=2 and k=n-2
    double expected_amplitude = n / 2.0;

    for (int i = 0; i < n; i++) {
        if (i == 2 || i == n - 2) {
            if (abs(abs(result[i]) - expected_amplitude) > TEST_TOLERANCE) {
                printf("1D Real input test failed at index %d\n", i);
                return false;
            }
        } else {
            if (abs(result[i]) > TEST_TOLERANCE) {
                printf("1D Real input test failed at index %d: expected 0, got (%f, %f)\n",
                       i, result[i].real(), result[i].imag());
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Test Gaussian function
 * Gaussian in space -> Gaussian in frequency
 * For discrete case, we test the inverse property
 */
bool test_1d_gaussian() {
    int n = 32;
    vector<complex<double>> data(n);

    double sigma = 2.0;
    double center = n / 2.0;

    // Create Gaussian centered at n/2
    for (int i = 0; i < n; i++) {
        double x = i - center;
        data[i] = exp(-x * x / (2.0 * sigma * sigma));
    }

    vector<complex<double>> original = data;
    vector<complex<double>> fft_result = fft_1d(data, false, false);
    vector<complex<double>> ifft_result = fft_1d(fft_result, true, false);

    // Test roundtrip
    for (int i = 0; i < n; i++) {
        if (abs(ifft_result[i] - original[i]) > TEST_TOLERANCE) {
            printf("1D Gaussian test failed at index %d\n", i);
            return false;
        }
    }

    return true;
}

// ============================================================================
// FLOAT FFT TESTS
// ============================================================================

const float TEST_TOLERANCE_F = 1e-5f;

/**
 * @brief Test 1D float FFT with delta function
 */
bool test_1d_delta_f() {
    int n = 8;
    vector<complex<float>> data(n, 0.0f);
    data[0] = 1.0f;

    vector<complex<float>> result = fft_1d_f(data, false, false);

    for (int i = 0; i < n; i++) {
        if (abs(result[i] - complex<float>(1.0f, 0.0f)) > TEST_TOLERANCE_F) {
            printf("1D Delta (float) test failed at index %d\n", i);
            return false;
        }
    }
    return true;
}

/**
 * @brief Test 1D float FFT inverse
 */
bool test_1d_inverse_f() {
    int n = 16;
    vector<complex<float>> data(n);

    for (int i = 0; i < n; i++) {
        data[i] = complex<float>(sinf(2.0f * M_PI * i / n), cosf(2.0f * M_PI * i / n));
    }

    vector<complex<float>> original = data;
    vector<complex<float>> fft_result = fft_1d_f(data, false, false);
    vector<complex<float>> ifft_result = fft_1d_f(fft_result, true, false);

    for (int i = 0; i < n; i++) {
        if (abs(ifft_result[i] - original[i]) > TEST_TOLERANCE_F) {
            printf("1D Inverse (float) test failed at index %d\n", i);
            return false;
        }
    }
    return true;
}

/**
 * @brief Test 1D float FFT with real input
 */
bool test_1d_real_input_f() {
    int n = 8;
    vector<float> data(n);

    for (int i = 0; i < n; i++) {
        data[i] = cosf(2.0f * M_PI * 2 * i / n);
    }

    vector<complex<float>> result = fft_1d_real_f(data, false);

    float expected_amplitude = n / 2.0f;

    for (int i = 0; i < n; i++) {
        if (i == 2 || i == n - 2) {
            if (abs(abs(result[i]) - expected_amplitude) > TEST_TOLERANCE_F) {
                printf("1D Real input (float) test failed at index %d\n", i);
                return false;
            }
        } else {
            if (abs(result[i]) > TEST_TOLERANCE_F) {
                printf("1D Real input (float) test failed at index %d: expected 0\n", i);
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Test 2D float FFT inverse
 */
bool test_2d_inverse_f() {
    int nx = 8, ny = 8;
    vector<complex<float>> data(nx * ny);

    for (int i = 0; i < nx * ny; i++) {
        data[i] = complex<float>(sinf(i * 0.5f), cosf(i * 0.3f));
    }

    vector<complex<float>> original = data;
    vector<complex<float>> fft_result = fft_2d_f(data, nx, ny, false, false);
    vector<complex<float>> ifft_result = fft_2d_f(fft_result, nx, ny, true, false);

    for (int i = 0; i < nx * ny; i++) {
        if (abs(ifft_result[i] - original[i]) > TEST_TOLERANCE_F) {
            printf("2D Inverse (float) test failed at index %d\n", i);
            return false;
        }
    }
    return true;
}

/**
 * @brief Test 3D float FFT inverse
 */
bool test_3d_inverse_f() {
    int nx = 4, ny = 4, nz = 4;
    vector<complex<float>> data(nx * ny * nz);

    for (int i = 0; i < nx * ny * nz; i++) {
        data[i] = complex<float>(sinf(i * 0.2f), cosf(i * 0.1f));
    }

    vector<complex<float>> original = data;
    vector<complex<float>> fft_result = fft_3d_f(data, nx, ny, nz, false, false);
    vector<complex<float>> ifft_result = fft_3d_f(fft_result, nx, ny, nz, true, false);

    for (int i = 0; i < nx * ny * nz; i++) {
        if (abs(ifft_result[i] - original[i]) > TEST_TOLERANCE_F) {
            printf("3D Inverse (float) test failed at index %d\n", i);
            return false;
        }
    }
    return true;
}

/**
 * @brief Test nD float FFT with 4D data
 */
bool test_4d_inverse_f() {
    vector<int> dims = {4, 4, 4, 2};
    int total_size = 4 * 4 * 4 * 2;
    vector<complex<float>> data(total_size);

    for (int i = 0; i < total_size; i++) {
        data[i] = complex<float>(sinf(i * 0.3f), cosf(i * 0.2f));
    }

    vector<complex<float>> original = data;
    vector<complex<float>> fft_result = fft_nd_f(data, dims, false, false);
    vector<complex<float>> ifft_result = fft_nd_f(fft_result, dims, true, false);

    for (int i = 0; i < total_size; i++) {
        if (abs(ifft_result[i] - original[i]) > TEST_TOLERANCE_F) {
            printf("4D Inverse (float) test failed at index %d\n", i);
            return false;
        }
    }
    return true;
}

bool fft_tests() {
    int num_tests = 18;
    bool all_tests[num_tests] = {
        test_1d_delta(),
        test_1d_constant(),
        test_1d_inverse(),
        test_1d_cosine(),
        test_1d_in_place(),
        test_1d_real_input(),
        test_1d_gaussian(),
        test_2d_separable(),
        test_2d_inverse(),
        test_3d_delta(),
        test_3d_inverse(),
        test_4d_inverse(),
        test_1d_delta_f(),
        test_1d_inverse_f(),
        test_1d_real_input_f(),
        test_2d_inverse_f(),
        test_3d_inverse_f(),
        test_4d_inverse_f(),
    };
    return print_test_results(all_tests, num_tests, "FFT tests");
}
