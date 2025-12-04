#pragma once

#include <vector>
#include <complex>
#include <fftw3.h>

using namespace std;

// ============================================================================
// DOUBLE PRECISION FFT FUNCTIONS
// ============================================================================

/**
 * @brief Performs a 1D Fast Fourier Transform using FFTW (double precision).
 *
 * @param data Input vector of complex values
 * @param inverse If true, performs inverse FFT
 * @param in_place If true, modifies input data; if false, returns new data
 *
 * @return Vector of complex FFT coefficients
 */
vector<complex<double>> fft_1d(vector<complex<double>>& data, bool inverse = false, bool in_place = false);

/**
 * @brief Performs a 1D Fast Fourier Transform on real input using FFTW (double precision).
 *
 * @param data Input vector of real values
 * @param inverse If true, performs inverse FFT
 *
 * @return Vector of complex FFT coefficients
 */
vector<complex<double>> fft_1d_real(vector<double>& data, bool inverse = false);

/**
 * @brief Performs a 2D Fast Fourier Transform using FFTW (double precision).
 *
 * @param data Input 2D array of complex values (flattened: [nx * ny])
 * @param nx Number of points in x dimension
 * @param ny Number of points in y dimension
 * @param inverse If true, performs inverse FFT
 * @param in_place If true, modifies input data; if false, returns new data
 *
 * @return Flattened vector of complex FFT coefficients
 */
vector<complex<double>> fft_2d(vector<complex<double>>& data, int nx, int ny,
                                bool inverse = false, bool in_place = false);

/**
 * @brief Performs a 3D Fast Fourier Transform using FFTW (double precision).
 *
 * @param data Input 3D array of complex values (flattened: [nx * ny * nz])
 * @param nx Number of points in x dimension
 * @param ny Number of points in y dimension
 * @param nz Number of points in z dimension
 * @param inverse If true, performs inverse FFT
 * @param in_place If true, modifies input data; if false, returns new data
 *
 * @return Flattened vector of complex FFT coefficients
 */
vector<complex<double>> fft_3d(vector<complex<double>>& data, int nx, int ny, int nz,
                                bool inverse = false, bool in_place = false);

/**
 * @brief Performs an nD Fast Fourier Transform using FFTW (double precision).
 *
 * @param data Input nD array of complex values (flattened)
 * @param dims Vector of dimensions [n0, n1, ..., n_{d-1}]
 * @param inverse If true, performs inverse FFT
 * @param in_place If true, modifies input data; if false, returns new data
 *
 * @return Flattened vector of complex FFT coefficients
 */
vector<complex<double>> fft_nd(vector<complex<double>>& data, vector<int> dims,
                                bool inverse = false, bool in_place = false);

// ============================================================================
// SINGLE PRECISION (FLOAT) FFT FUNCTIONS
// ============================================================================

/**
 * @brief Performs a 1D Fast Fourier Transform using FFTW (single precision).
 *
 * @param data Input vector of complex values
 * @param inverse If true, performs inverse FFT
 * @param in_place If true, modifies input data; if false, returns new data
 *
 * @return Vector of complex FFT coefficients
 */
vector<complex<float>> fft_1d_f(vector<complex<float>>& data, bool inverse = false, bool in_place = false);

/**
 * @brief Performs a 1D Fast Fourier Transform on real input using FFTW (single precision).
 *
 * @param data Input vector of real values
 * @param inverse If true, performs inverse FFT
 *
 * @return Vector of complex FFT coefficients
 */
vector<complex<float>> fft_1d_real_f(vector<float>& data, bool inverse = false);

/**
 * @brief Performs a 2D Fast Fourier Transform using FFTW (single precision).
 *
 * @param data Input 2D array of complex values (flattened: [nx * ny])
 * @param nx Number of points in x dimension
 * @param ny Number of points in y dimension
 * @param inverse If true, performs inverse FFT
 * @param in_place If true, modifies input data; if false, returns new data
 *
 * @return Flattened vector of complex FFT coefficients
 */
vector<complex<float>> fft_2d_f(vector<complex<float>>& data, int nx, int ny,
                                 bool inverse = false, bool in_place = false);

/**
 * @brief Performs a 3D Fast Fourier Transform using FFTW (single precision).
 *
 * @param data Input 3D array of complex values (flattened: [nx * ny * nz])
 * @param nx Number of points in x dimension
 * @param ny Number of points in y dimension
 * @param nz Number of points in z dimension
 * @param inverse If true, performs inverse FFT
 * @param in_place If true, modifies input data; if false, returns new data
 *
 * @return Flattened vector of complex FFT coefficients
 */
vector<complex<float>> fft_3d_f(vector<complex<float>>& data, int nx, int ny, int nz,
                                 bool inverse = false, bool in_place = false);

/**
 * @brief Performs an nD Fast Fourier Transform using FFTW (single precision).
 *
 * @param data Input nD array of complex values (flattened)
 * @param dims Vector of dimensions [n0, n1, ..., n_{d-1}]
 * @param inverse If true, performs inverse FFT
 * @param in_place If true, modifies input data; if false, returns new data
 *
 * @return Flattened vector of complex FFT coefficients
 */
vector<complex<float>> fft_nd_f(vector<complex<float>>& data, vector<int> dims,
                                 bool inverse = false, bool in_place = false);
