/*
 * This file contains the definition of the ScalarField and VectorField classes.
 * These classes are used to represent scalar and vector fields, respectively.
 * They are implemented over a mesh of points and values, and can be used to interpolate
 * The mesh is interpolated automatically when the field is created. Ideal for BZ calculations
 *
 * Author: Griffin Heier
 */
#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

#include <string>
#include <complex>
#include "vec.hpp"

using namespace std;

class ScalarField {
    public:
        vector<Vec> points;
        vector<float> values;
        vector<float> ivalues;
        vector<Vec> domain;
        vector<Vec> inv_domain;
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx, ny, nz, nw;
        int dimension;
        bool is_complex;

        ScalarField();
        ScalarField(vector<Vec> points, vector<float> values, int dimension=3, bool is_complex=false);
        ScalarField(vector<Vec> points, vector<complex<float>> values, int dimension=3, bool is_complex=true);
        ScalarField(string filename, int dimension=3, bool is_complex=false);

        void get_values_for_interpolation();
        float operator() (Vec point);
        float operator()(py::array_t<double> array);
};

class VectorField {
    public:
        vector<Vec> points;
        vector<Vec> values;
        vector<Vec> ivalues;
        vector<Vec> domain;
        vector<Vec> inv_domain;
        float xmin, xmax, ymin, ymax, zmin, zmax, wmin, wmax;
        int nx, ny, nz, nw;
        int dimension;
        bool is_complex;

        VectorField();
        VectorField(vector<Vec> points, vector<Vec> values, int dimension=3, bool is_complex=false);
        VectorField(vector<Vec> points, vector<complex<float>> values, int dimension=3, bool is_complex=true);
        VectorField(string filename, int dimension=3, bool is_complex=false);

        void get_values_for_interpolation();
        Vec operator() (Vec point);
};

vector<float> matrix_multiplication(vector<Vec>& matrix, Vec& vec, int n);
Vec vec_matrix_multiplication(vector<float>& matrix, Vec& vec, int n);
bool domain_vec_found(float a[4], float dx, float dy, float dz, float dw);
vector<Vec> invertMatrix(vector<Vec>& matrix, int n);

