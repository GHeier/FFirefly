#pragma once

#include "vec.hpp"

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace std;

class __attribute__((visibility("default"))) Field {
public:
    py::object field_obj;
    int dimension;
    bool is_complex;
    bool is_vector;

    Field();
    Field(int dimension, bool is_complex, bool is_vector);
    py::module init();
    vector<double> Vec_to_vector(Vec point);
};

class ScalarField : public Field {
public:
    ScalarField();
    ScalarField(vector<Vec> points, vector<float> values, int dimension);
    ScalarField(vector<vector<double>> points, vector<float> values, int dimension);
    ScalarField(string filename, int dimension);
    float operator()(vector<double> coords);
    float operator()(Vec coords);
};

class ComplexField : public Field {
public:
    ComplexField();
    ComplexField(vector<Vec> points, vector<complex<float>> values, int dimension);
    ComplexField(vector<vector<double>> points, vector<complex<float>> values, int dimension);
    ComplexField(string filename, int dimension);
    complex<float> operator()(vector<double> coords);
    complex<float> operator()(Vec coords);
};

class VectorField : public Field {
public:
    VectorField();
    VectorField(vector<Vec> points, vector<float> values, int dimension);
    VectorField(vector<vector<double>> points, vector<float> values, int dimension);
    VectorField(string filename, int dimension);
    vector<float> operator()(vector<double> coords);
    vector<float> operator()(Vec coords);
};

class ComplexVectorField : public Field {
public:
    ComplexVectorField();
    ComplexVectorField(vector<Vec> points, vector<complex<float>> values, int dimension);
    ComplexVectorField(vector<vector<double>> points, vector<complex<float>> values, int dimension);
    ComplexVectorField(string filename, int dimension);
    vector<complex<float>> operator()(vector<double> coords);
    vector<complex<float>> operator()(Vec coords);
};

vector<vector<double>> Vec_to_doubles(vector<Vec> points);
