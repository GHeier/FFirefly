/**
 * Array class - nD array with FFT support
 *
 * Author: Claude Code (with FFTW integration)
 */

#include "array.hpp"
#include "../config/load/cpp_config.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <numeric>

// Default constructor
Array::Array() {
    total_size = 0;
}

// Constructor with dimensions
Array::Array(vector<int> dimensions) : dims(dimensions) {
    // Validate dimensions
    for (int dim : dims) {
        if (dim <= 0) {
            throw invalid_argument("All dimensions must be positive");
        }
    }

    // Calculate total size
    total_size = 1;
    for (int dim : dims) {
        total_size *= dim;
    }

    printv("Array dimensions: [");
    for (size_t i = 0; i < dims.size(); i++) {
        printv("%d", dims[i]);
        if (i < dims.size() - 1) printv(", ");
    }
    printv("]\n");
    printv("Total size: %zu\n", total_size);

    float MB = (float)(total_size * sizeof(float)) / 1000000.0f;
    if (MB > 100) {
        printv("Allocating %.1f GB\n", MB / 1024);
    } else if (MB > 1) {
        printv("Allocating %.0f MB\n", MB);
    } else {
        printv("Allocating %.3f MB\n", MB);
    }

    try {
        vals.resize(total_size, 0.0f);
    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed: " << e.what() << '\n';
        exit(EXIT_FAILURE);
    }

    printv("Array allocated\n");
}

// Constructor with dimensions and fill value
Array::Array(vector<int> dimensions, float fill_value) : Array(dimensions) {
    fill(fill_value);
}

// Helper to flatten multi-dimensional index
size_t Array::flatten_index(const vector<int>& indices) const {
    if (indices.size() != dims.size()) {
        throw invalid_argument("Index dimensionality must match array dimensionality");
    }

    size_t flat_idx = 0;
    size_t stride = 1;

    // Row-major order (C-style)
    for (int i = dims.size() - 1; i >= 0; i--) {
        if (indices[i] < 0 || indices[i] >= dims[i]) {
            throw out_of_range("Index out of range");
        }
        flat_idx += indices[i] * stride;
        stride *= dims[i];
    }

    return flat_idx;
}

// Helper for 1D indexing
size_t Array::flatten_index(int i) const {
    if (dims.size() != 1) {
        throw invalid_argument("1D indexing requires 1D array");
    }
    if (i < 0 || i >= dims[0]) {
        throw out_of_range("Index out of range");
    }
    return i;
}

// Helper for 2D indexing
size_t Array::flatten_index(int i, int j) const {
    if (dims.size() != 2) {
        throw invalid_argument("2D indexing requires 2D array");
    }
    if (i < 0 || i >= dims[0] || j < 0 || j >= dims[1]) {
        throw out_of_range("Index out of range");
    }
    return i * dims[1] + j;
}

// Helper for 3D indexing
size_t Array::flatten_index(int i, int j, int k) const {
    if (dims.size() != 3) {
        throw invalid_argument("3D indexing requires 3D array");
    }
    if (i < 0 || i >= dims[0] || j < 0 || j >= dims[1] || k < 0 || k >= dims[2]) {
        throw out_of_range("Index out of range");
    }
    return i * dims[1] * dims[2] + j * dims[2] + k;
}

// Accessors
float &Array::operator()(vector<int> indices) {
    return vals[flatten_index(indices)];
}

float Array::operator()(vector<int> indices) const {
    return vals[flatten_index(indices)];
}

float &Array::operator()(int i) {
    return vals[flatten_index(i)];
}

float &Array::operator()(int i, int j) {
    return vals[flatten_index(i, j)];
}

float &Array::operator()(int i, int j, int k) {
    return vals[flatten_index(i, j, k)];
}

// Assignment operator (returns a copy)
Array &Array::operator=(const Array &other) {
    if (this != &other) {
        dims = other.dims;
        total_size = other.total_size;
        vals = other.vals; // This creates a copy
    }
    return *this;
}

// Addition
Array Array::operator+(const Array &other) const {
    if (dims != other.dims) {
        throw invalid_argument("Array dimensions must match for addition");
    }

    Array result(dims);
    for (size_t i = 0; i < total_size; i++) {
        result.vals[i] = vals[i] + other.vals[i];
    }
    return result;
}

Array &Array::operator+=(const Array &other) {
    if (dims != other.dims) {
        throw invalid_argument("Array dimensions must match for addition");
    }

    for (size_t i = 0; i < total_size; i++) {
        vals[i] += other.vals[i];
    }
    return *this;
}

// Subtraction
Array Array::operator-(const Array &other) const {
    if (dims != other.dims) {
        throw invalid_argument("Array dimensions must match for subtraction");
    }

    Array result(dims);
    for (size_t i = 0; i < total_size; i++) {
        result.vals[i] = vals[i] - other.vals[i];
    }
    return result;
}

Array &Array::operator-=(const Array &other) {
    if (dims != other.dims) {
        throw invalid_argument("Array dimensions must match for subtraction");
    }

    for (size_t i = 0; i < total_size; i++) {
        vals[i] -= other.vals[i];
    }
    return *this;
}

// Scalar multiplication
Array Array::operator*(float scalar) const {
    Array result(dims);
    for (size_t i = 0; i < total_size; i++) {
        result.vals[i] = vals[i] * scalar;
    }
    return result;
}

Array &Array::operator*=(float scalar) {
    for (size_t i = 0; i < total_size; i++) {
        vals[i] *= scalar;
    }
    return *this;
}

// Scalar division
Array Array::operator/(float scalar) const {
    if (scalar == 0.0f) {
        throw invalid_argument("Division by zero");
    }

    Array result(dims);
    for (size_t i = 0; i < total_size; i++) {
        result.vals[i] = vals[i] / scalar;
    }
    return result;
}

Array &Array::operator/=(float scalar) {
    if (scalar == 0.0f) {
        throw invalid_argument("Division by zero");
    }

    for (size_t i = 0; i < total_size; i++) {
        vals[i] /= scalar;
    }
    return *this;
}

// Comparison
bool Array::operator==(const Array &other) const {
    if (dims != other.dims) {
        return false;
    }

    for (size_t i = 0; i < total_size; i++) {
        if (abs(vals[i] - other.vals[i]) > 1e-6f) {
            return false;
        }
    }

    return true;
}

// Dot product (sum of element-wise products)
float Array::dot(const Array &other) const {
    if (dims != other.dims) {
        throw invalid_argument("Array dimensions must match for dot product");
    }

    float result = 0.0f;
    for (size_t i = 0; i < total_size; i++) {
        result += vals[i] * other.vals[i];
    }

    return result;
}

// FFT operations
vector<complex<float>> Array::fft(bool inverse, bool in_place) {
    return fft_call(inverse, in_place);
}

vector<complex<float>> Array::ifft(bool in_place) {
    return fft_call(true, in_place);
}

vector<complex<float>> Array::fft_call(bool inverse, bool in_place) {
    // Convert float data to complex
    vector<complex<float>> complex_data(total_size);
    for (size_t i = 0; i < total_size; i++) {
        complex_data[i] = complex<float>(vals[i], 0.0f);
    }

    // Perform FFT based on dimensionality
    if (dims.size() == 1) {
        return fft_1d_f(complex_data, inverse, in_place);
    } else if (dims.size() == 2) {
        return fft_2d_f(complex_data, dims[0], dims[1], inverse, in_place);
    } else if (dims.size() == 3) {
        return fft_3d_f(complex_data, dims[0], dims[1], dims[2], inverse, in_place);
    } else {
        return fft_nd_f(complex_data, dims, inverse, in_place);
    }
}

// FFT for real input (optimized)
vector<complex<float>> Array::fft_real(bool inverse) {
    if (dims.size() == 1) {
        return fft_1d_real_f(vals, inverse);
    } else {
        // For multi-dimensional, use general FFT
        return fft(inverse, false);
    }
}

// Fill array with a value
void Array::fill(float value) {
    for (size_t i = 0; i < total_size; i++) {
        vals[i] = value;
    }
}

// Print array (limited output for large arrays)
void Array::print() const {
    printf("Array with dimensions: [");
    for (size_t i = 0; i < dims.size(); i++) {
        printf("%d", dims[i]);
        if (i < dims.size() - 1) printf(", ");
    }
    printf("]\n");

    if (total_size <= 100) {
        // Print all elements for small arrays
        for (size_t i = 0; i < total_size; i++) {
            printf("%.4f ", vals[i]);
            if ((i + 1) % 10 == 0) printf("\n");
        }
        if (total_size % 10 != 0) printf("\n");
    } else {
        // Print summary for large arrays
        printf("First 10 elements: ");
        for (size_t i = 0; i < 10 && i < total_size; i++) {
            printf("%.4f ", vals[i]);
        }
        printf("...\n");
    }
}
