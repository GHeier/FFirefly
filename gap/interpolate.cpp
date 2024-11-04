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
    assert(x_val >= x_min && x_val <= x_max);
    assert(f.size() > 1);

    float dx = (x_max - x_min) / (f.size() - 1);
    int i = (x_val - x_min) / dx;
    assert(i >= 0 && i <= f.size() - 1);
    if (i == f.size() - 1) i--;

    float x_rel = (x_val - x_min) / dx - i;
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
    assert(x_val >= x_min && x_val <= x_max);
    assert(y_val >= y_min && y_val <= y_max);
    assert(f.size() > 1);
    assert(f[0].size() > 1);

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    assert(i >= 0 && i <= f.size() - 1);
    assert(j >= 0 && j <= f[0].size() - 1);
    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;


    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;

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
    assert(x_val >= x_min && x_val <= x_max);
    assert(y_val >= y_min && y_val <= y_max);
    assert(z_val >= z_min && z_val <= z_max);
    assert(f.size() > 1);
    assert(f[0].size() > 1);
    assert(f[0][0].size() > 1);

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);
    float dz = (z_max - z_min) / (f[0][0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    assert(i >= 0 && i <= f.size() - 1);
    assert(j >= 0 && j <= f[0].size() - 1);
    assert(k >= 0 && k <= f[0][0].size() - 1);
    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;
    if (k == f[0][0].size() - 1) k--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;

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
    assert(x_val >= x_min && x_val <= x_max);
    assert(y_val >= y_min && y_val <= y_max);
    assert(z_val >= z_min && z_val <= z_max);
    assert(w_val >= w_min && w_val <= w_max);
    assert(f.size() > 1);
    assert(f[0].size() > 1);
    assert(f[0][0].size() > 1);
    assert(f[0][0][0].size() > 1);

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);
    float dz = (z_max - z_min) / (f[0][0].size() - 1);
    float dw = (w_max - w_min) / (f[0][0][0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    int l = (w_val - w_min) / dw;
    assert(i >= 0 && i <= f.size() - 1);
    assert(j >= 0 && j <= f[0].size() - 1);
    assert(k >= 0 && k <= f[0][0].size() - 1);
    assert(l >= 0 && l <= f[0][0][0].size() - 1);
    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;
    if (k == f[0][0].size() - 1) k--;
    if (l == f[0][0][0].size() - 1) l--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    float w_rel = (w_val - w_min) / dw - l;

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
    assert(x_val >= x_min && x_val <= x_max);
    assert(f.size() > 1);

    float dx = (x_max - x_min) / (f.size() - 1);
    int i = (x_val - x_min) / dx;
    assert(i >= 0 && i <= f.size() - 1);
    if (i == f.size() - 1) i--;

    float x_rel = (x_val - x_min) / dx - i;
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
    assert(x_val >= x_min && x_val <= x_max);
    assert(y_val >= y_min && y_val <= y_max);
    assert(f.size() > 1);
    assert(f[0].size() > 1);

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    assert(i >= 0 && i <= f.size() - 1);
    assert(j >= 0 && j <= f[0].size() - 1);
    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;

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
    assert(x_val >= x_min && x_val <= x_max);
    assert(y_val >= y_min && y_val <= y_max);
    assert(z_val >= z_min && z_val <= z_max);
    assert(f.size() > 1);
    assert(f[0].size() > 1);
    assert(f[0][0].size() > 1);

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);
    float dz = (z_max - z_min) / (f[0][0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    assert(i >= 0 && i <= f.size() - 1);
    assert(j >= 0 && j <= f[0].size() - 1);
    assert(k >= 0 && k <= f[0][0].size() - 1);
    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;
    if (k == f[0][0].size() - 1) k--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;

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

    assert(x_val >= x_min && x_val <= x_max);
    assert(y_val >= y_min && y_val <= y_max);
    assert(z_val >= z_min && z_val <= z_max);
    assert(w_val >= w_min && w_val <= w_max);
    assert(f.size() > 1);
    assert(f[0].size() > 1);
    assert(f[0][0].size() > 1);
    assert(f[0][0][0].size() > 1);

    float dx = (x_max - x_min) / (f.size() - 1);
    float dy = (y_max - y_min) / (f[0].size() - 1);
    float dz = (z_max - z_min) / (f[0][0].size() - 1);
    float dw = (w_max - w_min) / (f[0][0][0].size() - 1);

    int i = (x_val - x_min) / dx;
    int j = (y_val - y_min) / dy;
    int k = (z_val - z_min) / dz;
    int l = (w_val - w_min) / dw;
    assert(i >= 0 && i <= f.size() - 1);
    assert(j >= 0 && j <= f[0].size() - 1);
    assert(k >= 0 && k <= f[0][0].size() - 1);
    assert(l >= 0 && l <= f[0][0][0].size() - 1);
    if (i == f.size() - 1) i--;
    if (j == f[0].size() - 1) j--;
    if (k == f[0][0].size() - 1) k--;
    if (l == f[0][0][0].size() - 1) l--;

    float x_rel = (x_val - x_min) / dx - i;
    float y_rel = (y_val - y_min) / dy - j;
    float z_rel = (z_val - z_min) / dz - k;
    float w_rel = (w_val - w_min) / dw - l;

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


