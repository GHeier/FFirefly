#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <complex>

#include "vec.h"
#include "interpolate.h"
#include "field.h"

using namespace std;

ScalarField::ScalarField() {
    points = vector<Vec>();
    values = vector<float>();
    ivalues = vector<float>();
    xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0, wmin = 0, wmax = 0;
    nx = 0, ny = 0, nz = 0, nw = 0;
    dimension = 0;
    is_complex = false;
}

VectorField::VectorField() {
    points = vector<Vec>();
    values = vector<Vec>();
    ivalues = vector<Vec>();
    xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0, wmin = 0, wmax = 0;
    nx = 0, ny = 0, nz = 0, nw = 0;
    dimension = 0;
}

void ScalarField::get_values_for_interpolation() {
    set<float> x_values;
    set<float> y_values;
    set<float> z_values;
    set<float> w_values;
    for (Vec point : points) {
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
    }
    nx = x_values.size();
    ny = y_values.size();
    nz = z_values.size();
    nw = w_values.size();
    dimension = 3;
}

void VectorField::get_values_for_interpolation() {
    set<float> x_values;
    set<float> y_values;
    set<float> z_values;
    set<float> w_values;
    for (Vec point : points) {
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
        dimension = point.dimension;
    }
    nx = x_values.size();
    ny = y_values.size();
    nz = z_values.size();
    nw = w_values.size();
}


ScalarField::ScalarField(vector<Vec> points, vector<float> values) {
    this->points = points;
    this->values = values;
    get_values_for_interpolation();
}

ScalarField::ScalarField(vector<Vec> points, vector<complex<float>> values) {
    this->points = points;
    transform(values.begin(), values.end(), values.begin(),
            [](complex<float> c) { return c.real(); });
    transform(values.begin(), values.end(), ivalues.begin(),
            [](complex<float> c) { return c.imag(); });
    get_values_for_interpolation();
    is_complex = true;
}

VectorField::VectorField(vector<Vec> points, vector<complex<float>> values) {
    this->points = points;
    transform(values.begin(), values.end(), values.begin(),
            [](complex<float> c) { return c.real(); });
    transform(values.begin(), values.end(), ivalues.begin(),
            [](complex<float> c) { return c.imag(); });
    get_values_for_interpolation();
    is_complex = true;
}

VectorField::VectorField(vector<Vec> points, vector<Vec> values) {
    this->points = points;
    this->values = values;
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

