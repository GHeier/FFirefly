/*
 * This file contains the implementation of the classes FastScalarField and FastVectorField.
 * They interpolate scalar and vector fields, respectively.
 * They have been implemented to allow for input via a file or a list of points and values.
 * They are constructed to exist over a mesh, so it is ideal for Brillouin Zone calculations.
 *
 * Author: Griffin Heier
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
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
#include "fastfield.hpp"

using namespace std;

FastScalarField::FastScalarField() {
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

FastVectorField::FastVectorField() {
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

// Function to invert a matrix represented as vector<Vec>
vector<Vec> invertMatrix(vector<Vec>& matrix, int n) {
    //print incoming matrix
    cout << "Incoming matrix: " << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i](j) << " ";
        }
        cout << endl;
    }

    // Create augmented matrix [A|I]
    vector<vector<float>> augmented(n, std::vector<float>(2 * n, 0.0f));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            augmented[i][j] = matrix[i](j);
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
            throw std::runtime_error("vector<vector<float>> is singular and cannot be inverted.");
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

vector<float> matrix_multiplication(vector<Vec>& matrix, Vec& vec, int n) {
    vector<float> result(n, 0.0f);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += matrix[i](j) * vec(j);
        }
    }
    return result;
}

Vec vec_matrix_multiplication(vector<Vec>& matrix, Vec& vec, int n) {
    Vec result;
    for (int i = 0; i < n; i++) {
        result(i) = 0;
        for (int j = 0; j < n; j++) {
            result(i) += matrix[i](j) * vec(j);
        }
    }
    return result;
}

void FastScalarField::get_values_for_interpolation() {
    int section = 0;
    vector<int> section_sizes = {1, 1, 1, 1};
    vector<Vec> domain;
    Vec first = points[0];
    Vec lattice_vec = (points[1] - first).round();
    int jump = 1;
    for (int i = 1; i < points.size(); i += jump) {
        Vec current_vec = (points[i] / section_sizes[section] - first).round();
        if (current_vec == lattice_vec) section_sizes[section]++;
        else {
            domain.push_back((points[i-1] - first).round());
            jump *= section_sizes[section];
            section++;
            section_sizes[section]++;
            if (section >= dimension) break;
            lattice_vec = (points[i] - first).round();
        }
    }
    domain.push_back((points[points.size()-1] - first).round());
    printf("Domain: \n");
    cout << "point: " << points[points.size()-1] << endl;
    cout << "first: " << first << endl << endl;
    for (int i = 0; i < domain.size(); i++) {
        for (int j = 0; j < i; j++) 
            domain[i] = domain[i] - domain[j];
    }
    inv_domain = invertMatrix(domain, dimension);
    nx = section_sizes[0];
    if (dimension > 1) ny = section_sizes[1];
    if (dimension > 2) nz = section_sizes[2];
    if (dimension > 3) nw = section_sizes[3];
}

void FastVectorField::get_values_for_interpolation() {
    int section = 0;
    vector<int> section_sizes = {1, 1, 1, 1};
    vector<Vec> domain;
    Vec first = points[0];
    Vec lattice_vec = (points[1] - first).round();
    int jump = 1;
    for (int i = 1; i < points.size(); i += jump) {
        Vec current_vec = (points[i] / section_sizes[section] - first).round();
        if (current_vec == lattice_vec) section_sizes[section]++;
        else {
            domain.push_back((points[i-1] - first).round());
            jump *= section_sizes[section];
            section++;
            section_sizes[section]++;
            if (section >= dimension) break;
            lattice_vec = (points[i] - first).round();
        }
    }
    domain.push_back((points[points.size()-1] - first).round());
    for (int i = 0; i < domain.size(); i++) {
        for (int j = 0; j < i; j++) 
            domain[i] = domain[i] - domain[j];
    }
    inv_domain = invertMatrix(domain, dimension);
    Vec min(0, 0, 0, 0);
    vector<float> min_values = matrix_multiplication(inv_domain, min, dimension);
    Vec max(1, 1, 1, 1);
    vector<float> max_values = matrix_multiplication(inv_domain, max, dimension);
    xmin = min_values[0]; xmax = max_values[0];
    if (dimension > 1) {
        ymin = min_values[1]; ymax = max_values[1];
    }
    if (dimension > 2) {
        zmin = min_values[2]; zmax = max_values[2];
    }
    if (dimension > 3) {
        wmin = min_values[3]; wmax = max_values[3];
    }
    nx = section_sizes[0];
    if (dimension > 1) ny = section_sizes[1];
    if (dimension > 2) nz = section_sizes[2];
    if (dimension > 3) nw = section_sizes[3];
}


FastScalarField::FastScalarField(vector<Vec> points, vector<float> values, int dimension, bool is_complex) {
    this->points = points;
    this->values = values;
    this->dimension = dimension;
    this->is_complex = is_complex;
    get_values_for_interpolation();
}

FastScalarField::FastScalarField(vector<Vec> points, vector<complex<float>> values, int dimension, bool is_complex) {
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

FastVectorField::FastVectorField(vector<Vec> points, vector<complex<float>> values, int dimension, bool is_complex) {
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

FastVectorField::FastVectorField(vector<Vec> points, vector<Vec> values, int dimension, bool is_complex) {
    this->points = points;
    this->values = values;
    this->dimension = dimension;
    this->is_complex = is_complex;
    get_values_for_interpolation();
}

float FastScalarField::operator() (Vec point) {
    Vec p = vec_matrix_multiplication(inv_domain, point, dimension);
    if (dimension == 1) return interpolate_1D(p.x, 0, 1, values);
    if (dimension == 2) return interpolate_2D(p.x, p.y, 0, 1, 0, 1, nx, ny, values);
    if (dimension == 3) return interpolate_3D(p.x, p.y, p.z, 0, 1, 0, 1, 0, 1, nx, ny, nz, values);
    if (dimension == 4) return interpolate_4D(p.x, p.y, p.z, p.w, 0, 1, 0, 1, 0, 1, 0, 1, nx, ny, nz, nw, values);
    else {
        cerr << "Invalid dimension." << endl;
        exit(1);
    }
}

float FastScalarField::operator()(py::array_t<double> array) {
        float dim = array.size();
        auto buf = array.unchecked<1>();
        Vec point; point.dimension = dim;
        for (int i = 0; i < dim; i++) {
            point(i) = buf(i);
        }
        return (*this)(point);
}

FastScalarField::FastScalarField(string filename, int dimension, bool is_complex) {
    this->dimension = dimension;
    this->is_complex = is_complex;
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Failed to open the file.");
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
        point.dimension = dimension;
        points.push_back(point);
        values.push_back(value);
        if (is_complex) {
            iss >> ivalue;
            ivalues.push_back(ivalue);
        }
    }
    get_values_for_interpolation();
}

FastVectorField::FastVectorField(string filename, int dimension, bool is_complex) {
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
        point.dimension = dimension;
        points.push_back(point);
        values.push_back(value);
        if (is_complex) 
            ivalues.push_back(ivalue);
    }
    get_values_for_interpolation();
}

