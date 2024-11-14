#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <set>

#include "vec.h"
#include "interpolate.h"
#include "field.h"

using namespace std;

ScalarField::ScalarField() {
    points = vector<Vec>();
    values = vector<float>();
    xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0, wmin = 0, wmax = 0;
    nx = 0, ny = 0, nz = 0, nw = 0;
    int dimension = 0;
}

VectorField::VectorField() {
    points = vector<Vec>();
    values = vector<Vec>();
    xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0, wmin = 0, wmax = 0;
    nx = 0, ny = 0, nz = 0, nw = 0;
    int dimension = 0;
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

VectorField::VectorField(vector<Vec> points, vector<Vec> values) {
    this->points = points;
    this->values = values;
    get_values_for_interpolation();
}

float ScalarField::operator() (Vec point) {
    if (dimension == 1) return interpolate_1D(point.x, xmin, xmax, values);
    if (dimension == 2) return interpolate_2D(point.x, point.y, xmin, xmax, ymin, ymax, nx, ny, values);
    if (dimension == 3) return interpolate_3D(point.x, point.y, point.z, xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, values);
    return interpolate_4D(point.x, point.y, point.z, point.w, xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax, nx, ny, nz, nw, values);
}

ScalarField::ScalarField(string filename, int dimension) {
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
        for (int i = 0; i < dimension; i++) {
            float temp;
            iss >> temp;

        }
        iss >> value;
        points.push_back(point);
        values.push_back(value);
    }
    get_values_for_interpolation();
}

VectorField::VectorField(string filename, int dimension) {
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
        for (int i = 0; i < dimension; i++) {
            float temp;
            iss >> temp;
            point(i) = temp;
        }
        for (int i = 0; i < dimension; i++) {
            float temp;
            iss >> temp;
            value(i) = temp;
        }
        points.push_back(point);
        values.push_back(value);
    }
    get_values_for_interpolation();
}

