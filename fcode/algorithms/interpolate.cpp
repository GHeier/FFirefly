#include <iostream>
#include <vector>
#include <cassert>
#include <complex>


using namespace std;

float sanitize_within_bounds(float value, float min_bound, float max_bound, float tolerance = 1e-5) {
    if (value > max_bound && value - max_bound < tolerance) {
        return max_bound;  // Adjust down to max bound
    } else if (value < min_bound && min_bound - value < tolerance) {
        return min_bound;  // Adjust up to min bound
    }
    return value;  // Return value as-is if within bounds
}

float interpolate_1D(float x_val, float x_min, float x_max, vector<float> &f) {
    // Interpolates a 1D function f(x) given on a grid between x_min and x_max
    // at the point x_val using linear interpolation
    // x_val: value at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (f.size() < 2) throw invalid_argument("f size too small");

    float dx = (x_max - x_min) / (f.size() - 1);
    int i = (x_val - x_min) / dx;
    if (i < 0 || i >= f.size()) throw out_of_range("i out of bounds");
    if (i == f.size() - 1) i--;

    float x_rel = (x_val - x_min) / dx - i;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    float result = f[i] + x_rel * (f[i + 1] - f[i]);

    return result;
}

float interpolate_2D(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, vector<vector<float>> &f) {
    // Interpolates a 2D function f(x, y) given on a grid between x_min and x_max
    // and y_min and y_max at the point (x_val, y_val) using bilinear interpolation
    // x_val, y_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (f.size() < 2) throw invalid_argument("f size too small");
    if (f[0].size() < 2) throw invalid_argument("f[0] size too small");

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    if (i < 0 || i >= f.size()) throw out_of_range("i out of bounds");
    if (j < 0 || j >= f[0].size()) throw out_of_range("j out of bounds");

    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;


    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");

    float result = (1 - x_rel) * (1 - y_rel) * f[i][j] +
                   x_rel * (1 - y_rel) * f[i + 1][j] +
                   (1 - x_rel) * y_rel * f[i][j + 1] +
                   x_rel * y_rel * f[i + 1][j + 1];

    return result;
}

float interpolate_3D(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, vector<vector<vector<float>>> &f) {
    // Interpolates a 3D function f(x, y, z) given on a grid between x_min and x_max
    // and y_min and y_max and z_min and z_max at the point (x_val, y_val, z_val) using trilinear interpolation
    // x_val, y_val, z_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // z_min, z_max: minimum and maximum values of z
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val, z_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    z_val = sanitize_within_bounds(z_val, z_min, z_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (z_val < z_min || z_val > z_max) throw out_of_range("z_val out of bounds");
    if (f.size() < 2) throw invalid_argument("f size too small");
    if (f[0].size() < 2) throw invalid_argument("f[0] size too small");
    if (f[0][0].size() < 2) throw invalid_argument("f[0][0] size too small");

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);
    float dz = (z_max - z_min) / (f[0][0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    if (i < 0 || i >= f.size()) throw out_of_range("i out of bounds");
    if (j < 0 || j >= f[0].size()) throw out_of_range("j out of bounds");
    if (k < 0 || k >= f[0][0].size()) throw out_of_range("k out of bounds");

    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;
    if (k == f[0][0].size() - 1) k--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");
    if (z_rel < 0 || z_rel > 1) throw out_of_range("z_rel out of bounds");

    float result = (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * f[i][j][k] +
                   x_rel * (1 - y_rel) * (1 - z_rel) * f[i + 1][j][k] +
                   (1 - x_rel) * y_rel * (1 - z_rel) * f[i][j + 1][k] +
                   x_rel * y_rel * (1 - z_rel) * f[i + 1][j + 1][k] +
                   (1 - x_rel) * (1 - y_rel) * z_rel * f[i][j][k + 1] +
                   x_rel * (1 - y_rel) * z_rel * f[i + 1][j][k + 1] +
                   (1 - x_rel) * y_rel * z_rel * f[i][j + 1][k + 1] +
                   x_rel * y_rel * z_rel * f[i + 1][j + 1][k + 1];

    return result;
}

float interpolate_4D(float x_val, float y_val, float z_val, float w_val, 
        float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, 
        float w_min, float w_max, vector<vector<vector<vector<float>>>> &f) {
    // Interpolates a 4D function f(x, y, z, w) given on a grid between x_min and x_max
    // and y_min and y_max and z_min and z_max and w_min and w_max at the point (x_val, y_val, z_val, w_val) using trilinear interpolation
    // x_val, y_val, z_val, w_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // z_min, z_max: minimum and maximum values of z
    // w_min, w_max: minimum and maximum values of w
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val, z_val, w_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    z_val = sanitize_within_bounds(z_val, z_min, z_max);
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (z_val < z_min || z_val > z_max) throw out_of_range("z_val out of bounds");
    if (w_val < w_min || w_val > w_max) throw out_of_range("w_val out of bounds");
    if (f.size() < 2) throw invalid_argument("f size too small");
    if (f[0].size() < 2) throw invalid_argument("f[0] size too small");
    if (f[0][0].size() < 2) throw invalid_argument("f[0][0] size too small");
    if (f[0][0][0].size() < 2) throw invalid_argument("f[0][0][0] size too small");

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);
    float dz = (z_max - z_min) / (f[0][0].size() - 1);
    float dw = (w_max - w_min) / (f[0][0][0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    int l = (w_val - w_min) / dw;
    if (i < 0 || i >= f.size()) throw out_of_range("i out of bounds");
    if (j < 0 || j >= f[0].size()) throw out_of_range("j out of bounds");
    if (k < 0 || k >= f[0][0].size()) throw out_of_range("k out of bounds");
    if (l < 0 || l >= f[0][0][0].size()) throw out_of_range("l out of bounds");

    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;
    if (k == f[0][0].size() - 1) k--;
    if (l == f[0][0][0].size() - 1) l--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    float w_rel = (w_val - w_min) / dw - l;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");
    if (z_rel < 0 || z_rel > 1) throw out_of_range("z_rel out of bounds");
    if (w_rel < 0 || w_rel > 1) throw out_of_range("w_rel out of bounds");

float result =
        (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * (1 - w_rel) * f[i][j][k][l] +
        x_rel * (1 - y_rel) * (1 - z_rel) * (1 - w_rel) * f[i + 1][j][k][l] +
        (1 - x_rel) * y_rel * (1 - z_rel) * (1 - w_rel) * f[i][j + 1][k][l] +
        x_rel * y_rel * (1 - z_rel) * (1 - w_rel) * f[i + 1][j + 1][k][l] +
        (1 - x_rel) * (1 - y_rel) * z_rel * (1 - w_rel) * f[i][j][k + 1][l] +
        x_rel * (1 - y_rel) * z_rel * (1 - w_rel) * f[i + 1][j][k + 1][l] +
        (1 - x_rel) * y_rel * z_rel * (1 - w_rel) * f[i][j + 1][k + 1][l] +
        x_rel * y_rel * z_rel * (1 - w_rel) * f[i + 1][j + 1][k + 1][l] +
        (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * w_rel * f[i][j][k][l + 1] +
        x_rel * (1 - y_rel) * (1 - z_rel) * w_rel * f[i + 1][j][k][l + 1] +
        (1 - x_rel) * y_rel * (1 - z_rel) * w_rel * f[i][j + 1][k][l + 1] +
        x_rel * y_rel * (1 - z_rel) * w_rel * f[i + 1][j + 1][k][l + 1] +
        (1 - x_rel) * (1 - y_rel) * z_rel * w_rel * f[i][j][k + 1][l + 1] +
        x_rel * (1 - y_rel) * z_rel * w_rel * f[i + 1][j][k + 1][l + 1] +
        (1 - x_rel) * y_rel * z_rel * w_rel * f[i][j + 1][k + 1][l + 1] +
        x_rel * y_rel * z_rel * w_rel * f[i + 1][j + 1][k + 1][l + 1];

    return result;
}

complex<float> interpolate_1D_complex(float x_val, float x_min, float x_max, vector<complex<float>> &f) {
    // Interpolates a 1D complex function f(x) given on a grid between x_min and x_max
    // at the point x_val using linear interpolation
    // x_val: value at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (f.size() < 2) throw invalid_argument("f size too small");

    float dx = (x_max - x_min) / (f.size() - 1);
    int i = (x_val - x_min) / dx;
    if (i < 0 || i >= f.size()) throw out_of_range("i out of bounds");
    if (i == f.size() - 1) i--;

    float x_rel = (x_val - x_min) / dx - i;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    complex<float> result = f[i] + x_rel * (f[i + 1] - f[i]);

    return result;
}

complex<float> interpolate_2D_complex(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, vector<vector<complex<float>>> &f) {
    // Interpolates a 2D complex function f(x, y) given on a grid between x_min and x_max
    // and y_min and y_max at the point (x_val, y_val) using bilinear interpolation
    // x_val, y_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (f.size() < 2) throw invalid_argument("f size too small");
    if (f[0].size() < 2) throw invalid_argument("f[0] size too small");

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    if (i < 0 || i >= f.size()) throw out_of_range("i out of bounds");
    if (j < 0 || j >= f[0].size()) throw out_of_range("j out of bounds");

    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");

    complex<float> fx1 = (1 - x_rel) * f[i][j] + x_rel * f[i + 1][j];
    complex<float> fx2 = (1 - x_rel) * f[i][j + 1] + x_rel * f[i + 1][j + 1];

    complex<float> result = (1 - y_rel) * fx1 + y_rel * fx2;
    return result;

}

complex<float> interpolate_3D_complex(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, vector<vector<vector<complex<float>>>> &f) {
    // Interpolates a 3D complex function f(x, y, z) given on a grid between x_min and x_max
    // and y_min and y_max and z_min and z_max at the point (x_val, y_val, z_val) using trilinear interpolation
    // x_val, y_val, z_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // z_min, z_max: minimum and maximum values of z
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val, z_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    z_val = sanitize_within_bounds(z_val, z_min, z_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (z_val < z_min || z_val > z_max) throw out_of_range("z_val out of bounds");
    if (f.size() < 2) throw invalid_argument("f size too small");
    if (f[0].size() < 2) throw invalid_argument("f[0] size too small");
    if (f[0][0].size() < 2) throw invalid_argument("f[0][0] size too small");

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);
    float dz = (z_max - z_min) / (f[0][0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    if (i < 0 || i >= f.size()) throw out_of_range("i out of bounds");
    if (j < 0 || j >= f[0].size()) throw out_of_range("j out of bounds");
    if (k < 0 || k >= f[0][0].size()) throw out_of_range("k out of bounds");

    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;
    if (k == f[0][0].size() - 1) k--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");
    if (z_rel < 0 || z_rel > 1) throw out_of_range("z_rel out of bounds");

    complex<float> fx1 = (1 - x_rel) * f[i][j][k] + x_rel * f[i + 1][j][k];
    complex<float> fx2 = (1 - x_rel) * f[i][j+1][k] + x_rel * f[i + 1][j+1][k];
    complex<float> fx3 = (1 - x_rel) * f[i][j][k+1] + x_rel * f[i + 1][j][k+1];
    complex<float> fx4 = (1 - x_rel) * f[i][j+1][k+1] + x_rel * f[i + 1][j+1][k+1];

    complex<float> fy1 = (1 - y_rel) * fx1 + y_rel * fx2;
    complex<float> fy2 = (1 - y_rel) * fx3 + y_rel * fx4;

    complex<float> result = (1 - z_rel) * fy1 + z_rel * fy2;

    return result;
}

complex<float> interpolate_4D_complex(float x_val, float y_val, float z_val, float w_val, 
        float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, 
        float w_min, float w_max, vector<vector<vector<vector<complex<float>>>>> &f) {
    // Interpolates a 4D complex function f(x, y, z, w) given on a grid between x_min and x_max
    // and y_min and y_max and z_min and z_max and w_min and w_max at the point (x_val, y_val, z_val, w_val) using trilinear interpolation
    // x_val, y_val, z_val, w_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // z_min, z_max: minimum and maximum values of z
    // w_min, w_max: minimum and maximum values of w
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val, z_val, w_val)

    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (z_val < z_min || z_val > z_max) throw out_of_range("z_val out of bounds");
    if (w_val < w_min || w_val > w_max) throw out_of_range("w_val out of bounds");
    if (f.size() < 2) throw invalid_argument("f size too small");
    if (f[0].size() < 2) throw invalid_argument("f[0] size too small");
    if (f[0][0].size() < 2) throw invalid_argument("f[0][0] size too small");
    if (f[0][0][0].size() < 2) throw invalid_argument("f[0][0][0] size too small");

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);
    float dz = (z_max - z_min) / (f[0][0].size() - 1);
    float dw = (w_max - w_min) / (f[0][0][0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    int l = (w_val - w_min) / dw;
    if (i < 0 || i >= f.size()) throw out_of_range("i out of bounds");
    if (j < 0 || j >= f[0].size()) throw out_of_range("j out of bounds");
    if (k < 0 || k >= f[0][0].size()) throw out_of_range("k out of bounds");
    if (l < 0 || l >= f[0][0][0].size()) throw out_of_range("l out of bounds");
    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;
    if (k == f[0][0].size() - 1) k--;
    if (l == f[0][0][0].size() - 1) l--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    float w_rel = (w_val - w_min) / dw - l;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");
    if (z_rel < 0 || z_rel > 1) throw out_of_range("z_rel out of bounds");
    if (w_rel < 0 || w_rel > 1) throw out_of_range("w_rel out of bounds");

    complex<float> fx1 = f[i][j][k][l] * (1 - x_rel) + f[i + 1][j][k][l] * x_rel;
    complex<float> fx2 = f[i][j+1][k][l] * (1 - x_rel) + f[i + 1][j+1][k][l] * x_rel;
    complex<float> fx3 = f[i][j][k+1][l] * (1 - x_rel) + f[i+1][j][k+1][l] * x_rel;
    complex<float> fx4 = f[i][j+1][k+1][l] * (1 - x_rel) + f[i + 1][j + 1][k+1][l] * x_rel;
    complex<float> fx5 = f[i][j][k][l+1] * (1 - x_rel) + f[i + 1][j][k][l+1] * x_rel;
    complex<float> fx6 = f[i][j+1][k][l+1] * (1 - x_rel) + f[i + 1][j+1][k][l+1] * x_rel;
    complex<float> fx7 = f[i][j][k+1][l+1] * (1 - x_rel) + f[i+1][j][k+1][l+1] * x_rel;
    complex<float> fx8 = f[i][j+1][k+1][l+1] * (1 - x_rel) + f[i + 1][j + 1][k+1][l+1] * x_rel;

    complex<float> fy1 = (1 - y_rel) * fx1 + y_rel * fx2;
    complex<float> fy2 = (1 - y_rel) * fx3 + y_rel * fx4;
    complex<float> fy3 = (1 - y_rel) * fx5 + y_rel * fx6;
    complex<float> fy4 = (1 - y_rel) * fx7 + y_rel * fx8;

    complex<float> fz1 = (1 - z_rel) * fy1 + z_rel * fy2;
    complex<float> fz2 = (1 - z_rel) * fy3 + z_rel * fy4;

    complex<float> result = (1 - w_rel) * fz1 + w_rel * fz2;
    return result;
}


float interpolate_2D(float x_val, float y_val, float x_min, float x_max, float y_min, float y_max, int nx, int ny, vector<float> &f) {
    // Interpolates a 2D function f(x, y) given on a grid between x_min and x_max
    // and y_min and y_max at the point (x_val, y_val) using bilinear interpolation
    // x_val, y_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (f.size() != nx * ny) throw invalid_argument("f size does not match nx and ny");

    float dx = (x_max - x_min) / (nx - 1);
    float dy = (y_max - y_min) / (ny - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    if (i < 0 || i >= nx) throw out_of_range("i out of bounds");
    if (j < 0 || j >= ny) throw out_of_range("j out of bounds");
    if (i == nx - 1) i--;
    if (j == ny - 1) j--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");

    float result = (1 - x_rel) * (1 - y_rel) * f[i*nx+j] +
                   x_rel * (1 - y_rel) * f[(i + 1)*nx+j] +
                   (1 - x_rel) * y_rel * f[i*nx + j + 1] +
                   x_rel * y_rel * f[(i + 1)*nx+j + 1];

    return result;
}

float interpolate_3D(float x_val, float y_val, float z_val, float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, int nx, int ny, int nz, vector<float> &f) {
    // Interpolates a 3D function f(x, y, z) given on a grid between x_min and x_max
    // and y_min and y_max and z_min and z_max at the point (x_val, y_val, z_val) using trilinear interpolation
    // x_val, y_val, z_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // z_min, z_max: minimum and maximum values of z
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val, z_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    z_val = sanitize_within_bounds(z_val, z_min, z_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (z_val < z_min || z_val > z_max) throw out_of_range("z_val out of bounds");

    float dx = (x_max - x_min) / (nx - 1);
    float dy = (y_max - y_min) / (ny - 1);
    float dz = (z_max - z_min) / (nz - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    if (i < 0 || i >= nx) throw out_of_range("i out of bounds");
    if (j < 0 || j >= ny) throw out_of_range("j out of bounds");
    if (k < 0 || k >= nz) throw out_of_range("k out of bounds");

    if (i == nx - 1) i--;
    if (j == ny - 1) j--;
    if (k == nz - 1) k--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");
    if (z_rel < 0 || z_rel > 1) throw out_of_range("z_rel out of bounds");

    float result = (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * f[i*nx*nx+j*ny+k] +
                   x_rel * (1 - y_rel) * (1 - z_rel) * f[(i + 1)*nx*nx+j*ny+k] +
                   (1 - x_rel) * y_rel * (1 - z_rel) * f[i*nx*nx + (j + 1)*ny+k] +
                   x_rel * y_rel * (1 - z_rel) * f[(i + 1)*nx*nx+(j + 1)*ny+k] +
                   (1 - x_rel) * (1 - y_rel) * z_rel * f[i*nx*nx+j*ny+k + 1] +
                   x_rel * (1 - y_rel) * z_rel * f[(i + 1)*nx*nx+j*ny+k + 1] +
                   (1 - x_rel) * y_rel * z_rel * f[i*nx*nx + (j + 1)*ny+k + 1] +
                   x_rel * y_rel * z_rel * f[(i + 1)*nx*nx+(j + 1)*ny+k + 1];

    return result;
}

float interpolate_4D(float x_val, float y_val, float z_val, float w_val, 
        float x_min, float x_max, float y_min, float y_max, float z_min, float z_max, 
        float w_min, float w_max, int nx, int ny, int nz, int nw, vector<float> &f) {
    // Interpolates a 4D function f(x, y, z, w) given on a grid between x_min and x_max
    // and y_min and y_max and z_min and z_max and w_min and w_max at the point (x_val, y_val, z_val, w_val) using trilinear interpolation
    // x_val, y_val, z_val, w_val: values at which to interpolate
    // x_min, x_max: minimum and maximum values of x
    // y_min, y_max: minimum and maximum values of y
    // z_min, z_max: minimum and maximum values of z
    // w_min, w_max: minimum and maximum values of w
    // f: vector of function values at the grid points
    // returns: interpolated value of f(x_val, y_val, z_val, w_val)

    x_val = sanitize_within_bounds(x_val, x_min, x_max);
    y_val = sanitize_within_bounds(y_val, y_min, y_max);
    z_val = sanitize_within_bounds(z_val, z_min, z_max);
    w_val = sanitize_within_bounds(w_val, w_min, w_max);
    if (x_val < x_min || x_val > x_max) throw out_of_range("x_val out of bounds");
    if (y_val < y_min || y_val > y_max) throw out_of_range("y_val out of bounds");
    if (z_val < z_min || z_val > z_max) throw out_of_range("z_val out of bounds");
    if (w_val < w_min || w_val > w_max) throw out_of_range("w_val out of bounds");
    int n = pow(f.size(), 0.25);

    float dx = (x_max - x_min) / (nx - 1);
    float dy = (y_max - y_min) / (ny - 1);
    float dz = (z_max - z_min) / (nz - 1);
    float dw = (w_max - w_min) / (nw - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    int l = (w_val - w_min) / dw;
    if (i < 0 || i >= nx) throw out_of_range("i out of bounds");
    if (j < 0 || j >= ny) throw out_of_range("j out of bounds");
    if (k < 0 || k >= nz) throw out_of_range("k out of bounds");
    if (l < 0 || l >= nw) throw out_of_range("l out of bounds");
    if (i == nx - 1) i--;
    if (j == ny - 1) j--;
    if (k == nz - 1) k--;
    if (l == nw - 1) l--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    float w_rel = (w_val - w_min) / dw - l;
    if (x_rel < 0 || x_rel > 1) throw out_of_range("x_rel out of bounds");
    if (y_rel < 0 || y_rel > 1) throw out_of_range("y_rel out of bounds");
    if (z_rel < 0 || z_rel > 1) throw out_of_range("z_rel out of bounds");
    if (w_rel < 0 || w_rel > 1) throw out_of_range("w_rel out of bounds");

float result = 
        (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * (1 - w_rel) * f[i*nx*nx*nx+j*ny*ny+k*nz+l] + 
        x_rel * (1 - y_rel) * (1 - z_rel) * (1 - w_rel) * f[(i + 1)*nx*nx*nx+j*ny*ny+k*nz+l] + 
        (1 - x_rel) * y_rel * (1 - z_rel) * (1 - w_rel) * f[i*nx*nx*nx + (j + 1)*ny*ny+k*nz+l] +
        x_rel * y_rel * (1 - z_rel) * (1 - w_rel) * f[(i + 1)*nx*nx*nx+(j + 1)*ny*ny+k*nz+l] +
        (1 - x_rel) * (1 - y_rel) * z_rel * (1 - w_rel) * f[i*nx*nx*nx+j*ny*ny+k*nz+l + 1] +
        x_rel * (1 - y_rel) * z_rel * (1 - w_rel) * f[(i + 1)*nx*nx*nx+j*ny*ny+k*nz+l + 1] +
        (1 - x_rel) * y_rel * z_rel * (1 - w_rel) * f[i*nx*nx*nx + (j + 1)*ny*ny+k*nz+l + 1] +
        x_rel * y_rel * z_rel * (1 - w_rel) * f[(i + 1)*nx*nx*nx+(j + 1)*ny*ny+k*nz+l + 1] +
        (1 - x_rel) * (1 - y_rel) * (1 - z_rel) * w_rel * f[i*nx*nx*nx+j*ny*ny+k*nz+l + 1] +
        x_rel * (1 - y_rel) * (1 - z_rel) * w_rel * f[(i + 1)*nx*nx*nx+j*ny*ny+k*nz+l + 1] +
        (1 - x_rel) * y_rel * (1 - z_rel) * w_rel * f[i*nx*nx*nx + (j + 1)*ny*ny+k*nz+l + 1] +
        x_rel * y_rel * (1 - z_rel) * w_rel * f[(i + 1)*nx*nx*nx+(j + 1)*ny*ny+k*nz+l + 1] +
        (1 - x_rel) * (1 - y_rel) * z_rel * w_rel * f[i*nx*nx*nx+j*ny*ny+k*nz+l + 1] +
        x_rel * (1 - y_rel) * z_rel * w_rel * f[(i + 1)*nx*nx*nx+j*ny*ny+k*nz+l + 1] +
        (1 - x_rel) * y_rel * z_rel * w_rel * f[i*nx*nx*nx + (j + 1)*ny*ny+k*nz+l + 1] +
        x_rel * y_rel * z_rel * w_rel * f[(i + 1)*nx*nx*nx+(j + 1)*ny*ny+k*nz+l + 1];

    return result;
}
