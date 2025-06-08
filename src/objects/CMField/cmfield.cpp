/*
 * This file contains the implementation of the class CMField
 * This class interpolates BZ(n+1 dimension) fields. It has been implemented
 * to allow for input via a file or a list of points and values. They are
 * constructed to exist over a mesh, so it is ideal for Brillouin Zone
 * calculations.
 *
 * Author: Griffin Heier
 */
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../../algorithms/interpolate.hpp"
#include "../CMData/cmdata.hpp"
#include "../vec.hpp"
#include "cmfield.hpp"

using namespace std;

CMField::CMField() {
    data = CMData();
    values = vector<vector<complex<Vec>>>();

    domain = vector<Vec>();
    inv_domain = vector<Vec>();
    first = Vec();
    nx = 0, ny = 0, nz = 0, nw = 0;
    wmax = 0, wmin = 0;
}

CMField::~CMField() {
    // No manual cleanup is necessary here unless you allocate memory with new.
    // Standard containers (like vector) and member objects (like CMData) clean up automatically.

    // Optional: clear vectors explicitly (only needed for large memory footprint or debugging)
    values.clear();
    domain.clear();
    inv_domain.clear();
}

CMField::CMField(CMData &data) {
    this->data = data;
    if (data.dimension == 1 or (data.dimension == 0 and data.with_w))
        find_domain_1d(data);
    else
        find_domain(data);
    if (data.with_w) {
        wmin = *min_element(data.w_points.begin(), data.w_points.end());
        wmax = *max_element(data.w_points.begin(), data.w_points.end());
    }
    inv_domain = invertMatrix(domain, data.dimension);

    make_values_2d(values, data);
    empty_CMData(data);
}

CMField::CMField(vector<Vec> points, vector<complex<Vec>> values, int dimension,
                 bool with_w, bool with_n, bool is_complex, bool is_vector) {
    data = CMData(points, values, dimension, with_w, with_n, is_complex,
                  is_vector);
    if (data.dimension == 1 or (data.dimension == 0 and data.with_w))
        find_domain_1d(data);
    else
        find_domain(data);
    if (data.with_w) {
        wmin = *min_element(data.w_points.begin(), data.w_points.end());
        wmax = *max_element(data.w_points.begin(), data.w_points.end());
    }
    inv_domain = invertMatrix(domain, data.dimension);

    make_values_2d(this->values, data);
    empty_CMData(data);
}

void empty_CMData(CMData &data) {
    data.points = vector<Vec>();
    // data.w_points = vector<float>();
    data.values = vector<complex<Vec>>();
    data.n_inds = vector<int>();
}

void make_values_2d(vector<vector<complex<Vec>>> &values, CMData &data) {
    if (!data.with_n) {
        values.push_back(data.values);
        return;
    }
    for (int i = 0; i < data.n_inds.size(); i++) {
        int start = data.n_inds[i];
        int end = data.n_inds[i + 1];
        if (i + 1 >= data.n_inds.size())
            end = data.values.size();
        if (end == 0)
            break;
        values.emplace_back(data.values.begin() + start,
                            data.values.begin() + end);
    }
}

float sizeof_values(vector<vector<complex<Vec>>> &values) {
    float MB = 0;
    for (auto x : values) {
        MB += 1.0 * x.size() * 64 / (1024 * 1024);
    }
    return MB;
}

CMField load_CMField(string filename) {
    CMData data = load(filename);
    CMField temp(data);
    return temp;
}

void save_CMField(string filename, CMField &field) {
    save(field.data, filename);
}

// Function to invert a matrix represented as vector<Vec>
vector<Vec> invertMatrix(vector<Vec> &matrix, int n) {
    // Create augmented matrix [A|I]
    vector<vector<float>> augmented(n, std::vector<float>(2 * n, 0.0f));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmented[i][j] = matrix[i](j);
            if (fabs(augmented[i][j]) < 1e-5)
                augmented[i][j] = 0.0f;
        }
        augmented[i][n + i] = 1.0f; // Identity matrix
    }
    // Gaussian Elimination
    for (size_t i = 0; i < n; ++i) {
        // Partial Pivoting
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::fabs(augmented[k][i]) > std::fabs(augmented[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(augmented[i], augmented[maxRow]);

        // Check for singular matrix
        if (std::fabs(augmented[i][i]) < 1e-6) {
            throw std::runtime_error(
                "vector<vector<float>> is singular and cannot be inverted.");
        }

        // Normalize the pivot row
        float pivot = augmented[i][i];
        for (size_t j = 0; j < 2 * n; ++j) {
            augmented[i][j] /= pivot;
        }

        // Eliminate other rows
        for (size_t k = 0; k < n; ++k) {
            if (k != i) {
                float factor = augmented[k][i];
                for (size_t j = 0; j < 2 * n; ++j) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }

    // Extract the inverted matrix
    vector<Vec> inv(n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            inv[i](j) = augmented[i][n + j];
        }
    }

    return inv;
}

vector<float> matrix_multiplication(vector<Vec> &matrix, Vec &vec, int n) {
    vector<float> result(n, 0.0f);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += matrix[i](j) * vec(j);
        }
    }
    return result;
}

Vec vec_matrix_multiplication(vector<Vec> &matrix, Vec &vec, int n) {
    Vec result;
    result.dimension = n;
    for (int i = 0; i < n; i++) {
        result(i) = 0;
        for (int j = 0; j < n; j++) {
            result(i) += matrix[i](j) * vec(j);
        }
    }
    return result;
}

Vec cross_product(Vec a, Vec b) {
    if (a.dimension == 3) {
        Vec result;
        result(0) = a(1) * b(2) - a(2) * b(1);
        result(1) = a(2) * b(0) - a(0) * b(2);
        result(2) = a(0) * b(1) - a(1) * b(0);
        return result;
    }
    if (a.dimension == 2) {
        Vec result;
        result(0) = a(0) * b(1) - a(1) * b(0);
        return result;
    }
    if (a.dimension == 1) {
        return a / a.norm();
    }
    throw runtime_error("Cross product only defined for 2D and 3D vectors.");
}

void CMField::find_domain_1d(CMData &data) {
    vector<Vec> domain;
    Vec first = data.points[0];
    int nw = data.w_points.size();
    int nn = data.n_inds.size();
    if (!data.with_w)
        nw = 1;
    if (!data.with_n)
        nn = 1;
    nx = data.points.size() / (nw * nn);
    Vec newpoint = (data.points[data.points.size() - 1] - first);
    newpoint.dimension = data.dimension;
    domain.push_back(newpoint);
    this->domain = domain;
    this->first = first;
    inv_domain = invertMatrix(domain, data.dimension);
}

void CMField::find_domain(CMData &data) {
    // Define initial number of points
    int section = 0;
    vector<int> section_sizes = {1, 1, 1};

    vector<Vec> domain;

    // Set initial point, vector, and hopping
    Vec first = data.points[0];
    int jump = 1;
    int w_size = data.w_points.size();
    if (w_size > 1)
        jump = w_size;

    int first_n = -1;
    if (data.n_inds.size() > 0)
        first_n = data.n_inds[1];

    Vec lattice_vec = (data.points[jump] - first);
    lattice_vec.dimension = data.dimension;

    for (int i = jump; i < data.points.size(); i += jump) {
        if (data.with_n and i >= first_n)
            break;
        Vec current_vec = (data.points[i] - first);
        if (cross_product(lattice_vec, current_vec).norm() <
            1e-6) { // If the same lattice vec
            section_sizes[section]++;
        } else { // Trigger if deviation from previous lattice vector
            Vec current_point = data.points[i];
            Vec p = (data.points[i - jump] - first);
            p.dimension = data.dimension;
            domain.push_back(p);
            jump = i;
            i -= jump;
            section++;
            if (section >= data.dimension)
                break;
            Vec temp = current_point - first;
            lattice_vec = (current_point - first);
            lattice_vec.dimension = data.dimension;
        }
    }
    // Save last lattice vector
    Vec last = data.points[data.points.size() - jump];
    last.dimension = data.dimension;
    domain.push_back((last - first));
    reverse(domain.begin(), domain.end());
    this->domain = domain;
    this->first = first;

    nx = section_sizes[0];
    // if (data.with_w) nx--;
    if (data.dimension > 1)
        ny = section_sizes[1];
    if (data.dimension > 2)
        nz = section_sizes[2];
    if (data.dimension > 3)
        nw = section_sizes[3];
}

double fold(double x) {
    // return x;
    double decimal = x - std::floor(x);
    return decimal;
}

void fold_to_first_BZ(Vec &p) {
    p.x = fold(p.x);
    p.y = fold(p.y);
    p.z = fold(p.z);
}

complex<Vec> CMField::operator()(float w) {
    if (!data.with_w)
        throw runtime_error("This field does not have a w dimension.");
    if (data.with_n)
        throw runtime_error("This field requires an index to be specified.");
    return CMF_search_1d(w, data.w_points, values[0]);
}

complex<Vec> CMField::operator()(int n, float w) {
    if (!data.with_w)
        throw runtime_error("This field does not have a w dimension.");
    if (!data.with_n)
        throw runtime_error("This field forbids an index to be specified.");
    if (n < 1 || n > values.size())
        throw out_of_range("Index starts at 1 and is out of bounds.");
    return CMF_search_1d(w, data.w_points, values[n - 1]);
}

complex<Vec> CMField::operator()(Vec point, float w) {
    if (w != 0 && !data.with_w)
        throw runtime_error("This field does not have a w dimension.");
    if (data.with_n)
        throw runtime_error("This field requires an index to be specified.");
    Vec shifted = point - first;
    Vec p = vec_matrix_multiplication(inv_domain, shifted, data.dimension);
    fold_to_first_BZ(p);
    p.w = w;
    if (!data.with_w) {
        if (data.dimension == 1)
            return interpolate_1D(p.x, 0, 1, values[0]);
        if (data.dimension == 2)
            return interpolate_2D(p.x, p.y, 0, 1, 0, 1, nx, ny, values[0]);
        if (data.dimension == 3)
            return interpolate_3D(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, nx, ny, nz,
                                  values[0]);
    } else {
        if (data.dimension == 1)
            return CMF_search_2d(p.x, w, nx, data.w_points, values[0]);
        if (data.dimension == 2)
            return CMF_search_3d(p.x, p.y, w, nx, ny, data.w_points, values[0]);
        if (data.dimension == 3)
            return CMF_search_4d(p.x, p.y, p.z, w, nx, ny, nz, data.w_points,
                                 values[0]);
    }
    throw runtime_error("Invalid dimension.");
}

complex<Vec> CMField::operator()(int n, Vec point, float w) {
    if (w != 0 && !data.with_w)
        throw runtime_error("This field does not have a w dimension.");
    if (!data.with_n)
        throw runtime_error("This field forbids an index to be specified.");
    if (n < 1 || n > values.size())
        throw out_of_range("Index starts at 1 and is out of bounds.");
    Vec shifted = point - first;
    Vec p = vec_matrix_multiplication(inv_domain, shifted, data.dimension);
    fold_to_first_BZ(p);
    p.w = w;
    if (!data.with_w) {
        if (data.dimension == 1)
            return interpolate_1D(p.x, 0, 1, values[n - 1]);
        if (data.dimension == 2)
            return interpolate_2D(p.x, p.y, 0, 1, 0, 1, nx, ny, values[n - 1]);
        if (data.dimension == 3)
            return interpolate_3D(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, nx, ny, nz,
                                  values[n - 1]);
    } else {
        if (data.dimension == 1)
            return CMF_search_2d(p.x, w, nx, data.w_points, values[n - 1]);
        if (data.dimension == 2)
            return CMF_search_3d(p.x, p.y, w, nx, ny, data.w_points,
                                 values[n - 1]);
        if (data.dimension == 3)
            return CMF_search_4d(p.x, p.y, p.z, w, nx, ny, nz, data.w_points,
                                 values[n - 1]);
    }
    throw runtime_error("Invalid dimension.");
}

complex<Vec> CMF_search_1d(float w_val, vector<float> &w_points,
                           vector<complex<Vec>> &f) {

    // Interpolates a 1D function f(w) given on a grid between w_min and w_max
    // at the point w_val using linear interpolation
    // w_val: value at which to interpolate
    // w_min, w_max: minimum and maximum values of w
    // f: vector of function values at the grid points
    // returns: interpolated value of f(w_val)
    float w_min = w_points[0], w_max = w_points[w_points.size() - 1];
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (w_val < w_min || w_val > w_max)
        throw out_of_range("w_val out of bounds");

    int i = binary_search(w_val, w_points);
    if (i == w_points.size() - 1)
        i--;

    if (w_val < w_points[i] || w_val > w_points[i + 1])
        throw out_of_range("Binary search failed\n");
    float w_rel = (w_val - w_points[i]) / (w_points[i + 1] - w_points[i]);
    if (w_rel < 0 || w_rel > 1)
        throw out_of_range("w_rel out of bounds");

    complex<Vec> result = (1 - w_rel) * f[i] + w_rel * f[i + 1];
    return result;
}

complex<Vec> CMF_search_2d(float x_val, float w_val, int nx,
                           vector<float> &w_points, vector<complex<Vec>> &f) {
    // Interpolates a 2D function f(x, y) given on a grid between x_min and
    // x_max and y_min and y_max at the point (x_val, y_val) using bilinear
    // interpolation x_val, y_val: values at which to interpolate x_min, x_max:
    // minimum and maximum values of x y_min, y_max: minimum and maximum values
    // of y f: vector of function values at the grid points returns:
    // interpolated value of f(x_val)
    float x_min = 0, x_max = 1;
    float w_min = w_points[0], w_max = w_points[w_points.size() - 1];
    int nw = w_points.size();

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (x_val < x_min || x_val > x_max)
        throw out_of_range("x_val out of bounds");
    if (w_val < w_min || w_val > w_max)
        throw out_of_range("w_val out of bounds");

    float dx = (x_max - x_min) / (nx - 1);

    int i = (x_val - x_min) / dx;
    int j = binary_search(w_val, w_points);
    if (i < 0 || i >= nx)
        throw out_of_range("i out of bounds");
    if (j < 0 || j >= w_points.size())
        throw out_of_range("j out of bounds");
    if (i == nx - 1)
        i--;
    if (j == w_points.size() - 1)
        j--;

    float x_rel = (x_val - x_min) / dx - i;
    float w_rel = (w_val - w_points[j]) / (w_points[j + 1] - w_points[j]);
    if (x_rel < 0 || x_rel > 1)
        throw out_of_range("x_rel out of bounds");
    if (w_rel < 0 || w_rel > 1)
        throw out_of_range("w_rel out of bounds");

    complex<Vec> result = (1 - x_rel) * (1 - w_rel) * f[i * nw + j] +
                          x_rel * (1 - w_rel) * f[(i + 1) * nw + j] +
                          (1 - x_rel) * w_rel * f[i * nw + j + 1] +
                          x_rel * w_rel * f[(i + 1) * nw + j + 1];

    return result;
}

complex<Vec> CMF_search_3d(float x_val, float y_val, float w_val, int nx,
                           int ny, vector<float> &w_points,
                           vector<complex<Vec>> &f) {
    // Interpolates a 3D function f(x, y, w) given on a grid between x_min and
    // x_max and y_min and y_max and w_min and w_max at the point (x_val, y_val,
    // w_val) using trilinear interpolation x_val, y_val, w_val: values at which
    // to interpolate x_min, x_max: minimum and maximum values of x y_min,
    // y_max: minimum and maximum values of y w_min, w_max: minimum and maximum
    // values of w f: vector of function values at the grid points returns:
    // interpolated value of f(x_val, y_val, w_val)
    float x_min = 0, x_max = 1;
    float y_min = 0, y_max = 1;
    float w_min = w_points[0], w_max = w_points[w_points.size() - 1];
    int nw = w_points.size();

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (x_val < x_min || x_val > x_max)
        throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max)
        throw out_of_range("y_val out of bounds");
    if (w_val < w_min || w_val > w_max)
        throw out_of_range("w_val out of bounds");

    float dx = (x_max - x_min) / (nx - 1);
    float dy = (y_max - y_min) / (ny - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = binary_search(w_val, w_points);
    if (i < 0 || i >= nx)
        throw out_of_range("i out of bounds");
    if (j < 0 || j >= ny)
        throw out_of_range("j out of bounds");
    if (i == nx - 1)
        i--;
    if (j == ny - 1)
        j--;
    if (k == nw - 1)
        k--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float w_rel = (w_val - w_points[k]) / (w_points[k + 1] - w_points[k]);
    if (x_rel < 0 || x_rel > 1)
        throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1)
        throw out_of_range("y_rel out of bounds");
    if (w_rel < 0 || w_rel > 1)
        throw out_of_range("w_rel out of bounds");

    complex<Vec> result =
        (1 - x_rel) * (1 - y_rel) * (1 - w_rel) * f[i * ny * nw + j * nw + k] +
        x_rel * (1 - y_rel) * (1 - w_rel) * f[(i + 1) * ny * nw + j * nw + k] +
        (1 - x_rel) * y_rel * (1 - w_rel) * f[i * ny * nw + (j + 1) * nw + k] +
        x_rel * y_rel * (1 - w_rel) * f[(i + 1) * ny * nw + (j + 1) * nw + k] +
        (1 - x_rel) * (1 - y_rel) * w_rel * f[i * ny * nw + j * nw + k + 1] +
        x_rel * (1 - y_rel) * w_rel * f[(i + 1) * ny * nw + j * nw + k + 1] +
        (1 - x_rel) * y_rel * w_rel * f[i * ny * nw + (j + 1) * nw + k + 1] +
        x_rel * y_rel * w_rel * f[(i + 1) * ny * nw + (j + 1) * nw + k + 1];

    return result;
}

complex<Vec> CMF_search_4d(float x_val, float y_val, float z_val, float w_val,
                           int nx, int ny, int nz, vector<float> &w_points,
                           vector<complex<Vec>> &f) {
    // Interpolates a 4D function f(x, y, z, w) given on a grid between x_min
    // and x_max and y_min and y_max and z_min and z_max and w_min and w_max at
    // the point (x_val, y_val, z_val, w_val) using trilinear interpolation
    // x_val, y_val, z_val, w_val: values at which to interpolate x_min, x_max:
    // minimum and maximum values of x y_min, y_max: minimum and maximum values
    // of y z_min, z_max: minimum and maximum values of z w_min, w_max: minimum
    // and maximum values of w f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val, z_val, w_val)
    float x_min = 0, x_max = 1;
    float y_min = 0, y_max = 1;
    float z_min = 0, z_max = 1;
    float w_min = w_points[0], w_max = w_points[w_points.size() - 1];
    int nw = w_points.size();

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    z_val = sanitize_within_bounds(z_val, z_min, z_max);
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (x_val < x_min || x_val > x_max)
        throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max)
        throw out_of_range("y_val out of bounds");
    if (z_val < z_min || z_val > z_max)
        throw out_of_range("z_val out of bounds");
    if (w_val < w_min || w_val > w_max)
        throw out_of_range("w_val out of bounds");

    float dx = (x_max - x_min) / (nx - 1);
    float dy = (y_max - y_min) / (ny - 1);
    float dz = (z_max - z_min) / (nz - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    int l = binary_search(w_val, w_points);
    if (i < 0 || i > nx)
        throw out_of_range("i out of bounds");
    if (j < 0 || j > ny)
        throw out_of_range("j out of bounds");
    if (k < 0 || k > nz)
        throw out_of_range("k out of bounds");
    if (i == nx - 1)
        i--;
    if (j == ny - 1)
        j--;
    if (k == nz - 1)
        k--;
    if (l == w_points.size() - 1)
        l--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    float w_rel = (w_val - w_points[l]) / (w_points[l + 1] - w_points[l]);
    if (x_rel < 0 || x_rel > 1)
        throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1)
        throw out_of_range("y_rel out of bounds");
    if (z_rel < 0 || z_rel > 1)
        throw out_of_range("z_rel out of bounds");
    if (w_rel < 0 || w_rel > 1)
        throw out_of_range("w_rel out of bounds");

    complex<Vec> result =
        (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * (1 - w_rel) *
            f[i * ny * nz * nw + j * nz * nw + k * nw + l] +
        x_rel * (1 - y_rel) * (1 - z_rel) * (1 - w_rel) *
            f[(i + 1) * ny * nz * nw + j * nz * nw + k * nw + l] +
        (1 - x_rel) * y_rel * (1 - z_rel) * (1 - w_rel) *
            f[i * ny * nz * nw + (j + 1) * nz * nw + k * nw + l] +
        x_rel * y_rel * (1 - z_rel) * (1 - w_rel) *
            f[(i + 1) * ny * nz * nw + (j + 1) * nz * nw + k * nw + l] +
        (1 - x_rel) * (1 - y_rel) * z_rel * (1 - w_rel) *
            f[i * ny * nz * nw + j * nz * nw + (k + 1) * nw + l] +
        x_rel * (1 - y_rel) * z_rel * (1 - w_rel) *
            f[(i + 1) * ny * nz * nw + j * nz * nw + (k + 1) * nw + l] +
        (1 - x_rel) * y_rel * z_rel * (1 - w_rel) *
            f[i * ny * nz * nw + (j + 1) * nz * nw + (k + 1) * nw + l] +
        x_rel * y_rel * z_rel * (1 - w_rel) *
            f[(i + 1) * ny * nz * nw + (j + 1) * nz * nw + (k + 1) * nw + l] +
        (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * w_rel *
            f[i * ny * nz * nw + j * nz * nw + k * nw + l + 1] +
        x_rel * (1 - y_rel) * (1 - z_rel) * w_rel *
            f[(i + 1) * ny * nz * nw + j * nz * nw + k * nw + l + 1] +
        (1 - x_rel) * y_rel * (1 - z_rel) * w_rel *
            f[i * ny * nz * nw + (j + 1) * nz * nw + k * nw + l + 1] +
        x_rel * y_rel * (1 - z_rel) * w_rel *
            f[(i + 1) * ny * nz * nw + (j + 1) * nz * nw + k * nw + l + 1] +
        (1 - x_rel) * (1 - y_rel) * z_rel * w_rel *
            f[i * ny * nz * nw + j * nz * nw + (k + 1) * nw + l + 1] +
        x_rel * (1 - y_rel) * z_rel * w_rel *
            f[(i + 1) * ny * nz * nw + j * nz * nw + (k + 1) * nw + l + 1] +
        (1 - x_rel) * y_rel * z_rel * w_rel *
            f[i * ny * nz * nw + (j + 1) * nz * nw + (k + 1) * nw + l + 1] +
        x_rel * y_rel * z_rel * w_rel *
            f[(i + 1) * ny * nz * nw + (j + 1) * nz * nw + (k + 1) * nw + l +
              1];

    return result;
}
