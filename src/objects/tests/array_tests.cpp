#include "src/config/load/c_config.h"
#include "src/objects/array.hpp"
#include "array_tests.hpp"
#include <cmath>
#include <iostream>

using namespace std;

const float ARRAY_TEST_TOL = 1e-5f;

/**
 * @brief Test 1D array creation and indexing
 */
bool test_array_1d() {
    Array arr({8});

    // Test initialization
    if (arr.size() != 8) {
        printf("Array 1D test failed: wrong size\n");
        return false;
    }

    // Test indexing
    for (int i = 0; i < 8; i++) {
        arr(i) = i * 2.0f;
    }

    for (int i = 0; i < 8; i++) {
        if (abs(arr(i) - i * 2.0f) > ARRAY_TEST_TOL) {
            printf("Array 1D test failed: indexing error\n");
            return false;
        }
    }

    return true;
}

/**
 * @brief Test 2D array creation and indexing
 */
bool test_array_2d() {
    Array arr({4, 4});

    // Test size
    if (arr.size() != 16) {
        printf("Array 2D test failed: wrong size\n");
        return false;
    }

    // Test indexing
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            arr(i, j) = i * 4.0f + j;
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (abs(arr(i, j) - (i * 4.0f + j)) > ARRAY_TEST_TOL) {
                printf("Array 2D test failed: indexing error\n");
                return false;
            }
        }
    }

    return true;
}

/**
 * @brief Test 3D array creation and indexing
 */
bool test_array_3d() {
    Array arr({2, 3, 4});

    // Test size
    if (arr.size() != 24) {
        printf("Array 3D test failed: wrong size\n");
        return false;
    }

    // Test indexing
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                arr(i, j, k) = i * 12.0f + j * 4.0f + k;
            }
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 4; k++) {
                float expected = i * 12.0f + j * 4.0f + k;
                if (abs(arr(i, j, k) - expected) > ARRAY_TEST_TOL) {
                    printf("Array 3D test failed: indexing error\n");
                    return false;
                }
            }
        }
    }

    return true;
}

/**
 * @brief Test array addition
 */
bool test_array_addition() {
    Array arr1({4});
    Array arr2({4});

    for (int i = 0; i < 4; i++) {
        arr1(i) = i;
        arr2(i) = i * 2.0f;
    }

    Array result = arr1 + arr2;

    for (int i = 0; i < 4; i++) {
        if (abs(result(i) - i * 3.0f) > ARRAY_TEST_TOL) {
            printf("Array addition test failed\n");
            return false;
        }
    }

    return true;
}

/**
 * @brief Test array subtraction
 */
bool test_array_subtraction() {
    Array arr1({4});
    Array arr2({4});

    for (int i = 0; i < 4; i++) {
        arr1(i) = i * 5.0f;
        arr2(i) = i * 2.0f;
    }

    Array result = arr1 - arr2;

    for (int i = 0; i < 4; i++) {
        if (abs(result(i) - i * 3.0f) > ARRAY_TEST_TOL) {
            printf("Array subtraction test failed\n");
            return false;
        }
    }

    return true;
}

/**
 * @brief Test scalar multiplication
 */
bool test_array_scalar_mult() {
    Array arr({4});

    for (int i = 0; i < 4; i++) {
        arr(i) = i + 1.0f;
    }

    Array result = arr * 3.0f;

    for (int i = 0; i < 4; i++) {
        if (abs(result(i) - (i + 1.0f) * 3.0f) > ARRAY_TEST_TOL) {
            printf("Array scalar multiplication test failed\n");
            return false;
        }
    }

    return true;
}

/**
 * @brief Test dot product
 */
bool test_array_dot() {
    Array arr1({4});
    Array arr2({4});

    for (int i = 0; i < 4; i++) {
        arr1(i) = i + 1.0f;
        arr2(i) = i + 1.0f;
    }

    // Dot product of [1,2,3,4] with itself = 1+4+9+16 = 30
    float result = arr1.dot(arr2);
    float expected = 30.0f;

    if (abs(result - expected) > ARRAY_TEST_TOL) {
        printf("Array dot product test failed: expected %f, got %f\n", expected, result);
        return false;
    }

    return true;
}

/**
 * @brief Test assignment operator returns a copy
 */
bool test_array_copy() {
    Array arr1({4});

    for (int i = 0; i < 4; i++) {
        arr1(i) = i * 2.0f;
    }

    Array arr2 = arr1; // Should create a copy

    // Modify original
    arr1(0) = 999.0f;

    // Check that copy is unchanged
    if (abs(arr2(0) - 0.0f) > ARRAY_TEST_TOL) {
        printf("Array copy test failed: copy was modified\n");
        return false;
    }

    // Check other elements match
    for (int i = 1; i < 4; i++) {
        if (abs(arr2(i) - i * 2.0f) > ARRAY_TEST_TOL) {
            printf("Array copy test failed: elements don't match\n");
            return false;
        }
    }

    return true;
}

/**
 * @brief Test 1D FFT on array
 */
bool test_array_fft_1d() {
    Array arr({8});

    // Create a simple cosine wave
    for (int i = 0; i < 8; i++) {
        arr(i) = cosf(2.0f * M_PI * 2 * i / 8);
    }

    vector<complex<float>> result = arr.fft_real(false);

    // Should have peaks at frequency 2 and 6 (8-2)
    float expected_amplitude = 4.0f; // n/2 for real cosine

    bool found_peaks = false;
    if (abs(abs(result[2]) - expected_amplitude) < ARRAY_TEST_TOL &&
        abs(abs(result[6]) - expected_amplitude) < ARRAY_TEST_TOL) {
        found_peaks = true;
    }

    if (!found_peaks) {
        printf("Array FFT 1D test failed: peaks not found\n");
        return false;
    }

    return true;
}

/**
 * @brief Test 2D FFT on array
 */
bool test_array_fft_2d() {
    Array arr({4, 4});

    // Fill with some data
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            arr(i, j) = sinf(i * 0.5f) + cosf(j * 0.3f);
        }
    }

    vector<complex<float>> fft_result = arr.fft(false, false);

    // Verify we got the right size
    if (fft_result.size() != 16) {
        printf("Array FFT 2D test failed: wrong output size\n");
        return false;
    }

    // Test that DC component has non-zero magnitude
    if (abs(fft_result[0]) < ARRAY_TEST_TOL) {
        printf("Array FFT 2D test failed: DC component is zero\n");
        return false;
    }

    return true;
}

/**
 * @brief Test 3D FFT on array
 */
bool test_array_fft_3d() {
    Array arr({2, 2, 2});

    // Fill with data
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                arr(i, j, k) = i + j * 2.0f + k * 4.0f;
            }
        }
    }

    vector<complex<float>> result = arr.fft(false, false);

    // Verify size
    if (result.size() != 8) {
        printf("Array FFT 3D test failed: wrong output size\n");
        return false;
    }

    return true;
}

bool array_tests() {
    int num_tests = 11;
    bool all_tests[num_tests] = {
        test_array_1d(),
        test_array_2d(),
        test_array_3d(),
        test_array_addition(),
        test_array_subtraction(),
        test_array_scalar_mult(),
        test_array_dot(),
        test_array_copy(),
        test_array_fft_1d(),
        test_array_fft_2d(),
        test_array_fft_3d(),
    };
    return print_test_results(all_tests, num_tests, "Array tests");
}
