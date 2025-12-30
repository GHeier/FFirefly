#pragma once

#include "src/algorithms/fft.hpp"
#include <vector>
#include <complex>
#include <iostream>

using namespace std;

/**
 * Array class
 *
 * This class represents an n-dimensional array of floats with FFT capabilities.
 * Data is stored in a flattened 1D vector and accessed via multi-dimensional indexing.
 * Supports mathematical operations and Fourier transforms.
 */
class Array {
  public:
    vector<int> dims;           // Dimensions of the array
    vector<float> vals;         // Flattened data storage
    size_t total_size;          // Total number of elements

    // Constructors
    Array();
    Array(vector<int> dimensions);
    Array(vector<int> dimensions, float fill_value);

    // Accessors
    float &operator()(vector<int> indices);
    float operator()(vector<int> indices) const;

    // For convenience, overload for 1D, 2D, 3D access
    float &operator()(int i);
    float &operator()(int i, int j);
    float &operator()(int i, int j, int k);

    // Arithmetic operators (element-wise)
    Array &operator=(const Array &other);
    Array operator+(const Array &other) const;
    Array operator-(const Array &other) const;
    Array operator*(float scalar) const;
    Array operator/(float scalar) const;
    Array &operator+=(const Array &other);
    Array &operator-=(const Array &other);
    Array &operator*=(float scalar);
    Array &operator/=(float scalar);

    // Comparison
    bool operator==(const Array &other) const;

    // Dot product (flattened element-wise product sum)
    float dot(const Array &other) const;

    // FFT operations (returns complex array as vector)
    vector<complex<float>> fft(bool inverse = false, bool in_place = true);
    vector<complex<float>> fft_real(bool inverse = false);
    vector<complex<float>> ifft(bool in_place = true);

    // Utility methods
    int ndim() const { return dims.size(); }
    size_t size() const { return total_size; }
    void fill(float value);
    void print() const;

  private:
    size_t flatten_index(const vector<int>& indices) const;
    size_t flatten_index(int i) const;
    size_t flatten_index(int i, int j) const;
    size_t flatten_index(int i, int j, int k) const;
    vector<complex<float>> fft_call(bool inverse, bool in_place);
};
