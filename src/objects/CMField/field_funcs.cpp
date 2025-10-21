/*
 * Field interpolation functions used by DataEvaluator and Field classes.
 * These functions perform multi-dimensional interpolation with frequency support.
 *
 * Author: Griffin Heier
 */
#include <algorithm>
#include <complex>
#include <stdexcept>
#include <vector>

#include "../../algorithms/interpolate.hpp"
#include "../vec.hpp"
#include "field_funcs.hpp"

using namespace std;

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

    //printf("i, j, k: %d, %d, %d\n", i, j, k);
    //printf("nx, ny, nw: %d, %d, %d\n", nx, ny, nw);
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

    // Special case: if nz=1, reduce to 3D interpolation
    if (nz == 1) {
        return CMF_search_3d(x_val, y_val, w_val, nx, ny, w_points, f);
    }

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
