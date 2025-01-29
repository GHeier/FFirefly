/*
 * This file contains the definition of the ScalarField and VectorField classes.
 * These classes are used to represent scalar and vector fields, respectively.
 * They are implemented over a mesh of points and values, and can be used to interpolate
 * The mesh is interpolated automatically when the field is created. Ideal for BZ calculations
 *
 * Author: Griffin Heier
 */
#pragma once

#include <memory>
#include <string>
#include <complex>
#include "vec.hpp"
#include "vec.hpp"

using namespace std;

class FastScalarField {
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

        FastScalarField();
        ~FastScalarField();
        FastScalarField(vector<Vec> points, vector<float> values, int dimension=3, bool is_complex=false);
        FastScalarField(vector<Vec> points, vector<complex<float>> values, int dimension=3, bool is_complex=true);
        FastScalarField(string filename, int dimension=3, bool is_complex=false);

        void get_values_for_interpolation();
        void spline_axis4();
        float operator() (Vec point);
private:
    struct Fit;
    unique_ptr<Fit> w_func;
};

class FastVectorField {
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

        FastVectorField();
        FastVectorField(vector<Vec> points, vector<Vec> values, int dimension=3, bool is_complex=false);
        FastVectorField(vector<Vec> points, vector<complex<float>> values, int dimension=3, bool is_complex=true);
        FastVectorField(string filename, int dimension=3, bool is_complex=false);

        void get_values_for_interpolation();
        Vec operator() (Vec point);
};

vector<float> matrix_multiplication(vector<Vec>& matrix, Vec& vec, int n);
Vec vec_matrix_multiplication(vector<float>& matrix, Vec& vec, int n);
bool domain_vec_found(float a[4], float dx, float dy, float dz, float dw);
vector<Vec> invertMatrix(vector<Vec>& matrix, int n);

