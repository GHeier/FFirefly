/*
 * This file contains the implementation of the classes ScalarField and VectorField.
 * They interpolate scalar and vector fields, respectively.
 * They have been implemented to allow for input via a file or a list of points and values.
 * They are constructed to exist over a mesh, so it is ideal for Brillouin Zone calculations.
 *
 * Author: Griffin Heier
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <complex>

#include "vec.hpp"
#include "../algorithms/interpolate.hpp"
#include "field.hpp"

using namespace std;

ScalarField::ScalarField() {
    points = vector<Vec>();
    values = vector<float>();
    ivalues = vector<float>();
    domain = vector<Vec>();
    inv_domain = vector<Vec>();
    xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0, wmin = 0, wmax = 0;
    nx = 0, ny = 0, nz = 0, nw = 0;
    dimension = 0;
    is_complex = false;
}

VectorField::VectorField() {
    points = vector<Vec>();
    values = vector<Vec>();
    ivalues = vector<Vec>();
    domain = vector<Vec>();
    inv_domain = vector<Vec>();
    xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0, wmin = 0, wmax = 0;
    nx = 0, ny = 0, nz = 0, nw = 0;
    dimension = 0;
    is_complex = false;
}

bool domain_vec_found(Vec a, float dx, float dy, float dz, float dw) {
    return a(0) == 0 and a(1) == 0 and a(2) == 0 and a(3) == 0 
        or fabs(a(0) - dx) < 1e-4 
        and fabs(a(1) - dy) < 1e-4 
        and fabs(a(2) - dz) < 1e-4 
        and fabs(a(3) - dw) < 1e-4;
}
// Function to invert a matrix represented as vector<Vec>
vector<Vec> invertMatrix(vector<Vec>& matrix, int n) {
    // Create augmented matrix
    vector<Vec> augmented(n, Vec(2 * n));
    for (int i = 0; i < n; ++i) {
        for (int j = 1; j <= n; ++j) {
            augmented[i](j) = matrix[i](j);         // Original matrix
            augmented[i](j + n) = (i == j - 1) ? 1 : 0; // Identity matrix
        }
    }

    // Perform Gaussian elimination
    for (int i = 0; i < n; ++i) {
        // Find pivot
        double pivot = augmented[i](i + 1);
        if (pivot == 0) {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }

        // Normalize pivot row
        for (int j = 1; j <= 2 * n; ++j) {
            augmented[i](j) /= pivot;
        }

        // Eliminate other rows
        for (int k = 0; k < n; ++k) {
            if (k == i) continue;
            double factor = augmented[k](i + 1);
            for (int j = 1; j <= 2 * n; ++j) {
                augmented[k](j) -= factor * augmented[i](j);
            }
        }
    }

    // Extract the inverse matrix
    std::vector<Vec> inverse(n, Vec(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 1; j <= n; ++j) {
            inverse[i](j) = augmented[i](j + n);
        }
    }

    return inverse;
}

void ScalarField::get_values_for_interpolation() {
    set<float> x_values;
    set<float> y_values;
    set<float> z_values;
    set<float> w_values;

    float dx, dy, dz, dw;
    Vec a(0,0,0,0,0,dimension);
    Vec b(0,0,0,0,0,dimension);
    Vec c(0,0,0,0,0,dimension);
    Vec d(0,0,0,0,0,dimension);
    int index = 0;
    for (Vec point : points) {
        if (index > 0 and index < points.size()) {
            dx = point.x - points[index-1].x;
            dy = point.y - points[index-1].y;
            dz = point.z - points[index-1].z;
            dw = point.w - points[index-1].w;
            if (domain_vec_found(a, dx, dy, dz, dw)) {
                a(0) += dx; a(1) += dy; a(2) += dz; a(3) += dw;
            }
            else if (domain_vec_found(b, dx, dy, dz, dw)) {
                b(0) += dx; b(1) += dy; b(2) += dz; b(3) += dw;
            }
            else if (domain_vec_found(c, dx, dy, dz, dw)) {
                c(0) += dx; c(1) += dy; c(2) += dz; c(3) += dw;
            }
            else if (domain_vec_found(d, dx, dy, dz, dw)) {
                d(0) += dx; d(1) += dy; d(2) += dz; d(3) += dw;
            }
        }
        if (point.x < xmin) xmin = point.x;
        if (point.x > xmax) xmax = point.x;
        if (point.y < ymin) ymin = point.y;
        if (point.y > ymax) ymax = point.y;
        if (point.z < zmin) zmin = point.z;
        if (point.z > zmax) zmax = point.z;
        if (point.w < wmin) wmin = point.w;
        if (point.w > wmax) wmax = point.w;
        x_values.insert(point.x);
        y_values.insert(point.y);
        z_values.insert(point.z);
        w_values.insert(point.w);
        index++;
    }
    nx = x_values.size();
    ny = y_values.size();
    nz = z_values.size();
    nw = w_values.size();
    if (dimension > 0) domain.push_back(a);
    if (dimension > 1) domain.push_back(b);
    if (dimension > 2) domain.push_back(c);
    if (dimension > 3) domain.push_back(d);
    inv_domain = invertMatrix(domain, dimension);
}

void VectorField::get_values_for_interpolation() {
    set<float> x_values;
    set<float> y_values;
    set<float> z_values;
    set<float> w_values;

    float dx, dy, dz, dw;
    Vec a(0,0,0,0,0,dimension);
    Vec b(0,0,0,0,0,dimension);
    Vec c(0,0,0,0,0,dimension);
    Vec d(0,0,0,0,0,dimension);
    int index = 0;
    for (Vec point : points) {
        if (index > 0 and index < points.size()) {
            dx = point.x - points[index-1].x;
            dy = point.y - points[index-1].y;
            dz = point.z - points[index-1].z;
            dw = point.w - points[index-1].w;
            if (domain_vec_found(a, dx, dy, dz, dw)) {
                a(0) += dx; a(1) += dy; a(2) += dz; a(3) += dw;
            }
            else if (domain_vec_found(b, dx, dy, dz, dw)) {
                b(0) += dx; b(1) += dy; b(2) += dz; b(3) += dw;
            }
            else if (domain_vec_found(c, dx, dy, dz, dw)) {
                c(0) += dx; c(1) += dy; c(2) += dz; c(3) += dw;
            }
            else if (domain_vec_found(d, dx, dy, dz, dw)) {
                d(0) += dx; d(1) += dy; d(2) += dz; d(3) += dw;
            }
        }
        if (point.x < xmin) xmin = point.x;
        if (point.x > xmax) xmax = point.x;
        if (point.y < ymin) ymin = point.y;
        if (point.y > ymax) ymax = point.y;
        if (point.z < zmin) zmin = point.z;
        if (point.z > zmax) zmax = point.z;
        if (point.w < wmin) wmin = point.w;
        if (point.w > wmax) wmax = point.w;
        x_values.insert(point.x);
        y_values.insert(point.y);
        z_values.insert(point.z);
        w_values.insert(point.w);
        index++;
    }
    nx = x_values.size();
    ny = y_values.size();
    nz = z_values.size();
    nw = w_values.size();
    if (dimension > 0) domain.push_back(a);
    if (dimension > 1) domain.push_back(b);
    if (dimension > 2) domain.push_back(c);
    if (dimension > 3) domain.push_back(d);
    inv_domain = invertMatrix(domain, dimension);
}


ScalarField::ScalarField(vector<Vec> points, vector<float> values, int dimension, bool is_complex) {
    this->points = points;
    this->values = values;
    this->dimension = dimension;
    this->is_complex = is_complex;
    get_values_for_interpolation();
}

ScalarField::ScalarField(vector<Vec> points, vector<complex<float>> values, int dimension, bool is_complex) {
    this->points = points;
    this->dimension = dimension;
    this->is_complex = is_complex;
    transform(values.begin(), values.end(), values.begin(),
            [](complex<float> c) { return c.real(); });
    transform(values.begin(), values.end(), ivalues.begin(),
            [](complex<float> c) { return c.imag(); });
    get_values_for_interpolation();
    is_complex = true;
}

VectorField::VectorField(vector<Vec> points, vector<complex<float>> values, int dimension, bool is_complex) {
    this->points = points;;
    this->dimension = dimension;
    this->is_complex = is_complex;
    transform(values.begin(), values.end(), values.begin(),
            [](complex<float> c) { return c.real(); });
    transform(values.begin(), values.end(), ivalues.begin(),
            [](complex<float> c) { return c.imag(); });
    get_values_for_interpolation();
    is_complex = true;
}

VectorField::VectorField(vector<Vec> points, vector<Vec> values, int dimension, bool is_complex) {
    this->points = points;
    this->values = values;
    this->dimension = dimension;
    this->is_complex = is_complex;
    get_values_for_interpolation();
}

float ScalarField::operator() (Vec point) {
    if (dimension == 1) return interpolate_1D(point.x, xmin, xmax, values);
    if (dimension == 2) return interpolate_2D(point.x, point.y, xmin, xmax, ymin, ymax, nx, ny, values);
    if (dimension == 3) return interpolate_3D(point.x, point.y, point.z, xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, values);
    if (dimension == 4) return interpolate_4D(point.x, point.y, point.z, point.w, xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax, nx, ny, nz, nw, values);
    else {
        cerr << "Invalid dimension." << endl;
        exit(1);
    }
}

ScalarField::ScalarField(string filename, int dimension, bool is_complex) {
    this->dimension = dimension;
    this->is_complex = is_complex;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open the file." << endl;
        exit(1);
    }
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        Vec point;
        float value;
        float ivalue;
        for (int i = 0; i < dimension; i++) {
            float temp;
            iss >> temp;
            point(i) = temp;
        }
        iss >> value;
        points.push_back(point);
        values.push_back(value);
        if (is_complex) {
            iss >> ivalue;
            ivalues.push_back(ivalue);
        }
    }
    get_values_for_interpolation();
}

VectorField::VectorField(string filename, int dimension, bool is_complex) {
    this->dimension = dimension;
    this->is_complex = is_complex;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open the file." << endl;
        exit(1);
    }
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        Vec point;
        Vec value;
        Vec ivalue;
        for (int i = 0; i < dimension; i++) {
            float temp;
            iss >> temp;
            point(i) = temp;
        }
        if (is_complex) {
            for (int i = 0; i < 2*dimension; i++) {
                float temp;
                iss >> temp;
                if (i % 2 == 0) value(i) = temp;
                else ivalue(i) = temp;
            }
        }
        else {
            for (int i = 0; i < dimension; i++) {
                float temp;
                iss >> temp;
                value(i) = temp;
            }
        }
        points.push_back(point);
        values.push_back(value);
        if (is_complex) 
            ivalues.push_back(ivalue);
    }
    get_values_for_interpolation();
}

