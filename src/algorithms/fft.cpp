#include "fft.hpp"
#include <cstring>
#include <stdexcept>

using namespace std;

// Helper function to convert std::complex to fftw_complex
static void copy_to_fftw(const vector<complex<double>>& src, fftw_complex* dst, size_t n) {
    for (size_t i = 0; i < n; i++) {
        dst[i][0] = src[i].real();
        dst[i][1] = src[i].imag();
    }
}

// Helper function to convert fftw_complex to std::complex
static void copy_from_fftw(const fftw_complex* src, vector<complex<double>>& dst, size_t n) {
    for (size_t i = 0; i < n; i++) {
        dst[i] = complex<double>(src[i][0], src[i][1]);
    }
}

vector<complex<double>> fft_1d(vector<complex<double>>& data, bool inverse, bool in_place) {
    int n = data.size();

    if (in_place) {
        // Allocate FFTW arrays
        fftw_complex* fftw_data = fftw_alloc_complex(n);

        // Copy input data
        copy_to_fftw(data, fftw_data, n);

        // Create plan and execute
        fftw_plan plan = fftw_plan_dft_1d(n, fftw_data, fftw_data,
                                          inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                          FFTW_ESTIMATE);
        fftw_execute(plan);

        // Copy back to input
        copy_from_fftw(fftw_data, data, n);

        // Normalize for inverse transform
        if (inverse) {
            for (int i = 0; i < n; i++) {
                data[i] /= n;
            }
        }

        // Clean up
        fftw_destroy_plan(plan);
        fftw_free(fftw_data);

        return data;
    } else {
        // Allocate FFTW arrays
        fftw_complex* fftw_in = fftw_alloc_complex(n);
        fftw_complex* fftw_out = fftw_alloc_complex(n);

        // Copy input data
        copy_to_fftw(data, fftw_in, n);

        // Create plan and execute
        fftw_plan plan = fftw_plan_dft_1d(n, fftw_in, fftw_out,
                                          inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                          FFTW_ESTIMATE);
        fftw_execute(plan);

        // Copy to result
        vector<complex<double>> result(n);
        copy_from_fftw(fftw_out, result, n);

        // Normalize for inverse transform
        if (inverse) {
            for (int i = 0; i < n; i++) {
                result[i] /= n;
            }
        }

        // Clean up
        fftw_destroy_plan(plan);
        fftw_free(fftw_in);
        fftw_free(fftw_out);

        return result;
    }
}

vector<complex<double>> fft_1d_real(vector<double>& data, bool inverse) {
    int n = data.size();

    if (!inverse) {
        // Real to complex transform
        double* fftw_in = fftw_alloc_real(n);
        fftw_complex* fftw_out = fftw_alloc_complex(n);

        // Copy input data
        for (int i = 0; i < n; i++) {
            fftw_in[i] = data[i];
        }

        // Create plan and execute
        fftw_plan plan = fftw_plan_dft_r2c_1d(n, fftw_in, fftw_out, FFTW_ESTIMATE);
        fftw_execute(plan);

        // Copy to result (r2c produces n/2+1 complex values, but we expand to full n)
        vector<complex<double>> result(n);
        int n_half = n / 2 + 1;

        // Copy positive frequencies
        for (int i = 0; i < n_half; i++) {
            result[i] = complex<double>(fftw_out[i][0], fftw_out[i][1]);
        }

        // Fill negative frequencies using Hermitian symmetry
        for (int i = n_half; i < n; i++) {
            result[i] = conj(result[n - i]);
        }

        // Clean up
        fftw_destroy_plan(plan);
        fftw_free(fftw_in);
        fftw_free(fftw_out);

        return result;
    } else {
        // Convert real to complex and use complex FFT
        vector<complex<double>> complex_data(n);
        for (int i = 0; i < n; i++) {
            complex_data[i] = complex<double>(data[i], 0.0);
        }
        return fft_1d(complex_data, inverse, false);
    }
}

vector<complex<double>> fft_2d(vector<complex<double>>& data, int nx, int ny,
                                bool inverse, bool in_place) {
    if (data.size() != (size_t)(nx * ny)) {
        throw invalid_argument("Data size must match nx * ny");
    }

    if (in_place) {
        // Allocate FFTW arrays
        fftw_complex* fftw_data = fftw_alloc_complex(nx * ny);

        // Copy input data
        copy_to_fftw(data, fftw_data, nx * ny);

        // Create plan and execute
        fftw_plan plan = fftw_plan_dft_2d(nx, ny, fftw_data, fftw_data,
                                          inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                          FFTW_ESTIMATE);
        fftw_execute(plan);

        // Copy back to input
        copy_from_fftw(fftw_data, data, nx * ny);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < data.size(); i++) {
                data[i] /= (nx * ny);
            }
        }

        // Clean up
        fftw_destroy_plan(plan);
        fftw_free(fftw_data);

        return data;
    } else {
        // Allocate FFTW arrays
        fftw_complex* fftw_in = fftw_alloc_complex(nx * ny);
        fftw_complex* fftw_out = fftw_alloc_complex(nx * ny);

        // Copy input data
        copy_to_fftw(data, fftw_in, nx * ny);

        // Create plan and execute
        fftw_plan plan = fftw_plan_dft_2d(nx, ny, fftw_in, fftw_out,
                                          inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                          FFTW_ESTIMATE);
        fftw_execute(plan);

        // Copy to result
        vector<complex<double>> result(nx * ny);
        copy_from_fftw(fftw_out, result, nx * ny);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < result.size(); i++) {
                result[i] /= (nx * ny);
            }
        }

        // Clean up
        fftw_destroy_plan(plan);
        fftw_free(fftw_in);
        fftw_free(fftw_out);

        return result;
    }
}

vector<complex<double>> fft_3d(vector<complex<double>>& data, int nx, int ny, int nz,
                                bool inverse, bool in_place) {
    if (data.size() != (size_t)(nx * ny * nz)) {
        throw invalid_argument("Data size must match nx * ny * nz");
    }

    if (in_place) {
        // Allocate FFTW arrays
        fftw_complex* fftw_data = fftw_alloc_complex(nx * ny * nz);

        // Copy input data
        copy_to_fftw(data, fftw_data, nx * ny * nz);

        // Create plan and execute
        fftw_plan plan = fftw_plan_dft_3d(nx, ny, nz, fftw_data, fftw_data,
                                          inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                          FFTW_ESTIMATE);
        fftw_execute(plan);

        // Copy back to input
        copy_from_fftw(fftw_data, data, nx * ny * nz);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < data.size(); i++) {
                data[i] /= (nx * ny * nz);
            }
        }

        // Clean up
        fftw_destroy_plan(plan);
        fftw_free(fftw_data);

        return data;
    } else {
        // Allocate FFTW arrays
        fftw_complex* fftw_in = fftw_alloc_complex(nx * ny * nz);
        fftw_complex* fftw_out = fftw_alloc_complex(nx * ny * nz);

        // Copy input data
        copy_to_fftw(data, fftw_in, nx * ny * nz);

        // Create plan and execute
        fftw_plan plan = fftw_plan_dft_3d(nx, ny, nz, fftw_in, fftw_out,
                                          inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                          FFTW_ESTIMATE);
        fftw_execute(plan);

        // Copy to result
        vector<complex<double>> result(nx * ny * nz);
        copy_from_fftw(fftw_out, result, nx * ny * nz);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < result.size(); i++) {
                result[i] /= (nx * ny * nz);
            }
        }

        // Clean up
        fftw_destroy_plan(plan);
        fftw_free(fftw_in);
        fftw_free(fftw_out);

        return result;
    }
}

vector<complex<double>> fft_nd(vector<complex<double>>& data, vector<int> dims,
                                bool inverse, bool in_place) {
    // Validate dimensions
    size_t total_size = 1;
    for (int dim : dims) {
        if (dim <= 0) {
            throw invalid_argument("All dimensions must be positive");
        }
        total_size *= dim;
    }

    if (data.size() != total_size) {
        throw invalid_argument("Data size must match product of dimensions");
    }

    int rank = dims.size();
    int* n = new int[rank];
    for (int i = 0; i < rank; i++) {
        n[i] = dims[i];
    }

    if (in_place) {
        // Allocate FFTW arrays
        fftw_complex* fftw_data = fftw_alloc_complex(total_size);

        // Copy input data
        copy_to_fftw(data, fftw_data, total_size);

        // Create plan and execute
        fftw_plan plan = fftw_plan_dft(rank, n, fftw_data, fftw_data,
                                       inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                       FFTW_ESTIMATE);
        fftw_execute(plan);

        // Copy back to input
        copy_from_fftw(fftw_data, data, total_size);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < data.size(); i++) {
                data[i] /= total_size;
            }
        }

        // Clean up
        fftw_destroy_plan(plan);
        fftw_free(fftw_data);
        delete[] n;

        return data;
    } else {
        // Allocate FFTW arrays
        fftw_complex* fftw_in = fftw_alloc_complex(total_size);
        fftw_complex* fftw_out = fftw_alloc_complex(total_size);

        // Copy input data
        copy_to_fftw(data, fftw_in, total_size);

        // Create plan and execute
        fftw_plan plan = fftw_plan_dft(rank, n, fftw_in, fftw_out,
                                       inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                       FFTW_ESTIMATE);
        fftw_execute(plan);

        // Copy to result
        vector<complex<double>> result(total_size);
        copy_from_fftw(fftw_out, result, total_size);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < result.size(); i++) {
                result[i] /= total_size;
            }
        }

        // Clean up
        fftw_destroy_plan(plan);
        fftw_free(fftw_in);
        fftw_free(fftw_out);
        delete[] n;

        return result;
    }
}

// ============================================================================
// SINGLE PRECISION (FLOAT) FFT FUNCTIONS
// ============================================================================

// Helper function to convert std::complex<float> to fftwf_complex
static void copy_to_fftwf(const vector<complex<float>>& src, fftwf_complex* dst, size_t n) {
    for (size_t i = 0; i < n; i++) {
        dst[i][0] = src[i].real();
        dst[i][1] = src[i].imag();
    }
}

// Helper function to convert fftwf_complex to std::complex<float>
static void copy_from_fftwf(const fftwf_complex* src, vector<complex<float>>& dst, size_t n) {
    for (size_t i = 0; i < n; i++) {
        dst[i] = complex<float>(src[i][0], src[i][1]);
    }
}

vector<complex<float>> fft_1d_f(vector<complex<float>>& data, bool inverse, bool in_place) {
    int n = data.size();

    if (in_place) {
        // Allocate FFTW arrays
        fftwf_complex* fftw_data = fftwf_alloc_complex(n);

        // Copy input data
        copy_to_fftwf(data, fftw_data, n);

        // Create plan and execute
        fftwf_plan plan = fftwf_plan_dft_1d(n, fftw_data, fftw_data,
                                            inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                            FFTW_ESTIMATE);
        fftwf_execute(plan);

        // Copy back to input
        copy_from_fftwf(fftw_data, data, n);

        // Normalize for inverse transform
        if (inverse) {
            for (int i = 0; i < n; i++) {
                data[i] /= n;
            }
        }

        // Clean up
        fftwf_destroy_plan(plan);
        fftwf_free(fftw_data);

        return data;
    } else {
        // Allocate FFTW arrays
        fftwf_complex* fftw_in = fftwf_alloc_complex(n);
        fftwf_complex* fftw_out = fftwf_alloc_complex(n);

        // Copy input data
        copy_to_fftwf(data, fftw_in, n);

        // Create plan and execute
        fftwf_plan plan = fftwf_plan_dft_1d(n, fftw_in, fftw_out,
                                            inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                            FFTW_ESTIMATE);
        fftwf_execute(plan);

        // Copy to result
        vector<complex<float>> result(n);
        copy_from_fftwf(fftw_out, result, n);

        // Normalize for inverse transform
        if (inverse) {
            for (int i = 0; i < n; i++) {
                result[i] /= n;
            }
        }

        // Clean up
        fftwf_destroy_plan(plan);
        fftwf_free(fftw_in);
        fftwf_free(fftw_out);

        return result;
    }
}

vector<complex<float>> fft_1d_real_f(vector<float>& data, bool inverse) {
    int n = data.size();

    if (!inverse) {
        // Real to complex transform
        float* fftw_in = fftwf_alloc_real(n);
        fftwf_complex* fftw_out = fftwf_alloc_complex(n);

        // Copy input data
        for (int i = 0; i < n; i++) {
            fftw_in[i] = data[i];
        }

        // Create plan and execute
        fftwf_plan plan = fftwf_plan_dft_r2c_1d(n, fftw_in, fftw_out, FFTW_ESTIMATE);
        fftwf_execute(plan);

        // Copy to result (r2c produces n/2+1 complex values, but we expand to full n)
        vector<complex<float>> result(n);
        int n_half = n / 2 + 1;

        // Copy positive frequencies
        for (int i = 0; i < n_half; i++) {
            result[i] = complex<float>(fftw_out[i][0], fftw_out[i][1]);
        }

        // Fill negative frequencies using Hermitian symmetry
        for (int i = n_half; i < n; i++) {
            result[i] = conj(result[n - i]);
        }

        // Clean up
        fftwf_destroy_plan(plan);
        fftwf_free(fftw_in);
        fftwf_free(fftw_out);

        return result;
    } else {
        // Convert real to complex and use complex FFT
        vector<complex<float>> complex_data(n);
        for (int i = 0; i < n; i++) {
            complex_data[i] = complex<float>(data[i], 0.0f);
        }
        return fft_1d_f(complex_data, inverse, false);
    }
}

vector<complex<float>> fft_2d_f(vector<complex<float>>& data, int nx, int ny,
                                 bool inverse, bool in_place) {
    if (data.size() != (size_t)(nx * ny)) {
        throw invalid_argument("Data size must match nx * ny");
    }

    if (in_place) {
        // Allocate FFTW arrays
        fftwf_complex* fftw_data = fftwf_alloc_complex(nx * ny);

        // Copy input data
        copy_to_fftwf(data, fftw_data, nx * ny);

        // Create plan and execute
        fftwf_plan plan = fftwf_plan_dft_2d(nx, ny, fftw_data, fftw_data,
                                            inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                            FFTW_ESTIMATE);
        fftwf_execute(plan);

        // Copy back to input
        copy_from_fftwf(fftw_data, data, nx * ny);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < data.size(); i++) {
                data[i] /= (nx * ny);
            }
        }

        // Clean up
        fftwf_destroy_plan(plan);
        fftwf_free(fftw_data);

        return data;
    } else {
        // Allocate FFTW arrays
        fftwf_complex* fftw_in = fftwf_alloc_complex(nx * ny);
        fftwf_complex* fftw_out = fftwf_alloc_complex(nx * ny);

        // Copy input data
        copy_to_fftwf(data, fftw_in, nx * ny);

        // Create plan and execute
        fftwf_plan plan = fftwf_plan_dft_2d(nx, ny, fftw_in, fftw_out,
                                            inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                            FFTW_ESTIMATE);
        fftwf_execute(plan);

        // Copy to result
        vector<complex<float>> result(nx * ny);
        copy_from_fftwf(fftw_out, result, nx * ny);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < result.size(); i++) {
                result[i] /= (nx * ny);
            }
        }

        // Clean up
        fftwf_destroy_plan(plan);
        fftwf_free(fftw_in);
        fftwf_free(fftw_out);

        return result;
    }
}

vector<complex<float>> fft_3d_f(vector<complex<float>>& data, int nx, int ny, int nz,
                                 bool inverse, bool in_place) {
    if (data.size() != (size_t)(nx * ny * nz)) {
        throw invalid_argument("Data size must match nx * ny * nz");
    }

    if (in_place) {
        // Allocate FFTW arrays
        fftwf_complex* fftw_data = fftwf_alloc_complex(nx * ny * nz);

        // Copy input data
        copy_to_fftwf(data, fftw_data, nx * ny * nz);

        // Create plan and execute
        fftwf_plan plan = fftwf_plan_dft_3d(nx, ny, nz, fftw_data, fftw_data,
                                            inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                            FFTW_ESTIMATE);
        fftwf_execute(plan);

        // Copy back to input
        copy_from_fftwf(fftw_data, data, nx * ny * nz);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < data.size(); i++) {
                data[i] /= (nx * ny * nz);
            }
        }

        // Clean up
        fftwf_destroy_plan(plan);
        fftwf_free(fftw_data);

        return data;
    } else {
        // Allocate FFTW arrays
        fftwf_complex* fftw_in = fftwf_alloc_complex(nx * ny * nz);
        fftwf_complex* fftw_out = fftwf_alloc_complex(nx * ny * nz);

        // Copy input data
        copy_to_fftwf(data, fftw_in, nx * ny * nz);

        // Create plan and execute
        fftwf_plan plan = fftwf_plan_dft_3d(nx, ny, nz, fftw_in, fftw_out,
                                            inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                            FFTW_ESTIMATE);
        fftwf_execute(plan);

        // Copy to result
        vector<complex<float>> result(nx * ny * nz);
        copy_from_fftwf(fftw_out, result, nx * ny * nz);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < result.size(); i++) {
                result[i] /= (nx * ny * nz);
            }
        }

        // Clean up
        fftwf_destroy_plan(plan);
        fftwf_free(fftw_in);
        fftwf_free(fftw_out);

        return result;
    }
}

vector<complex<float>> fft_nd_f(vector<complex<float>>& data, vector<int> dims,
                                 bool inverse, bool in_place) {
    // Validate dimensions
    size_t total_size = 1;
    for (int dim : dims) {
        if (dim <= 0) {
            throw invalid_argument("All dimensions must be positive");
        }
        total_size *= dim;
    }

    if (data.size() != total_size) {
        throw invalid_argument("Data size must match product of dimensions");
    }

    int rank = dims.size();
    int* n = new int[rank];
    for (int i = 0; i < rank; i++) {
        n[i] = dims[i];
    }

    if (in_place) {
        // Allocate FFTW arrays
        fftwf_complex* fftw_data = fftwf_alloc_complex(total_size);

        // Copy input data
        copy_to_fftwf(data, fftw_data, total_size);

        // Create plan and execute
        fftwf_plan plan = fftwf_plan_dft(rank, n, fftw_data, fftw_data,
                                         inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                         FFTW_ESTIMATE);
        fftwf_execute(plan);

        // Copy back to input
        copy_from_fftwf(fftw_data, data, total_size);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < data.size(); i++) {
                data[i] /= (float)total_size;
            }
        }

        // Clean up
        fftwf_destroy_plan(plan);
        fftwf_free(fftw_data);
        delete[] n;

        return data;
    } else {
        // Allocate FFTW arrays
        fftwf_complex* fftw_in = fftwf_alloc_complex(total_size);
        fftwf_complex* fftw_out = fftwf_alloc_complex(total_size);

        // Copy input data
        copy_to_fftwf(data, fftw_in, total_size);

        // Create plan and execute
        fftwf_plan plan = fftwf_plan_dft(rank, n, fftw_in, fftw_out,
                                         inverse ? FFTW_BACKWARD : FFTW_FORWARD,
                                         FFTW_ESTIMATE);
        fftwf_execute(plan);

        // Copy to result
        vector<complex<float>> result(total_size);
        copy_from_fftwf(fftw_out, result, total_size);

        // Normalize for inverse transform
        if (inverse) {
            for (size_t i = 0; i < result.size(); i++) {
                result[i] /= (float)total_size;
            }
        }

        // Clean up
        fftwf_destroy_plan(plan);
        fftwf_free(fftw_in);
        fftwf_free(fftw_out);
        delete[] n;

        return result;
    }
}
